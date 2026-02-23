#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Create mask for protein coding genes
"""
Merge VEP annotations with Illumina precomputed SpliceAI (SNV + INDEL),
then generate functional categories and a SAIGE group file.
"""

from __future__ import annotations

import argparse
import re
import shlex
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import DefaultDict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


# ==============================================================================
# Logging + small utilities
# ==============================================================================

def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def eprint(*args, **kwargs) -> None:
    print(f"[{ts()}]", *args, file=sys.stderr, **kwargs)


def die(msg: str, code: int = 1) -> None:
    eprint(f"[ERROR] {msg}")
    raise SystemExit(code)


def run_and_capture(cmd: Sequence[str]) -> str:
    """Run a command and return stdout (stripped)."""
    try:
        return subprocess.check_output(list(cmd), text=True).strip()
    except subprocess.CalledProcessError as ex:
        die(f"Command failed: {' '.join(cmd)}\n{ex}")


def check_exists(path: str, kind: str = "file") -> None:
    p = Path(path)
    if kind == "file" and not p.is_file():
        die(f"Missing file: {path}")
    if kind == "dir" and not p.is_dir():
        die(f"Missing directory: {path}")


def extract_tail_snpid(s: str) -> str:
    """
    From 'DRAGEN:chr20:60001:T:A' -> 'chr20:60001:T:A'
    Also handles 'DRAGEN:20:...' by enforcing 'chr' prefix in the returned ID.
    """
    if not isinstance(s, str):
        return s  # type: ignore[return-value]
    parts = s.split(":")
    if len(parts) >= 4:
        tail = ":".join(parts[-4:])
        chrom = tail.split(":")[0]
        if not chrom.startswith("chr"):
            tail_parts = tail.split(":")
            tail_parts[0] = "chr" + tail_parts[0]
            tail = ":".join(tail_parts)
        return tail
    return s


# ==============================================================================
# Tabix helpers
# ==============================================================================

def list_contigs(vcf_gz: str) -> List[str]:
    """
    List contigs from a bgzipped + tabix-indexed VCF.
    Requires: tabix on PATH.
    """
    try:
        out = run_and_capture(["tabix", "-l", vcf_gz])
    except Exception as ex:
        die(f"tabix failed on {vcf_gz}. Is it bgzipped & indexed (.tbi)?\n{ex}")
    contigs = [c for c in out.splitlines() if c]
    if not contigs:
        die(f"No contigs reported by `tabix -l` for {vcf_gz}.")
    return contigs


def normalize_region_arg(available: List[str], contig: str) -> str:
    """
    Return a contig name that exists in the file, accepting '22' or 'chr22'.
    """
    if contig in available:
        return contig
    if contig.startswith("chr") and contig[3:] in available:
        return contig[3:]
    if ("chr" + contig) in available:
        return "chr" + contig
    die(f"Contig '{contig}' not found in VCF index. Example contigs: {', '.join(available[:10])} ...")
    return contig  # unreachable


# ==============================================================================
# SpliceAI parsing
# ==============================================================================

_RE_SPLICE_PAYLOAD = re.compile(r"SpliceAI=([^;\t]+)")


def parse_spliceai_info(payload: str) -> List[Tuple[str, float]]:
    """
    Parse a SpliceAI INFO payload and return list of (gene, DSmax).

    The function:
    - Supports payloads with or without leading "SpliceAI=".
    - Each record is split by comma.
    - Each record is split by '|', and DSmax is max of the 4 DS fields
      (commonly DS_AG, DS_AL, DS_DG, DS_DL) in positions [2:6].

    Output is collapsed to max DS per gene across records.
    """
    if not isinstance(payload, str) or payload.strip() in {"", "."}:
        return []

    if payload.startswith("SpliceAI="):
        payload = payload.split("=", 1)[1]

    per_gene: DefaultDict[str, float] = defaultdict(float)

    for rec in payload.split(","):
        parts = rec.split("|")
        if len(parts) < 6:
            continue

        ds_vals: List[float] = []
        for x in parts[2:6]:
            try:
                ds_vals.append(float(x))
            except Exception:
                pass
        if not ds_vals:
            continue

        ds_max = float(np.max(ds_vals))

        gene_field = parts[1]
        # Some layouts have "SYMBOL---ENSG..." or other joined pieces
        if "---" in gene_field:
            gene_id: Optional[str] = None
            for chunk in gene_field.split("---"):
                if chunk.startswith("ENSG"):
                    gene_id = chunk.split(".")[0]
                    break
            gene = gene_id if gene_id else gene_field.split("---")[0]
        else:
            gene = gene_field

        if ds_max > per_gene[gene]:
            per_gene[gene] = ds_max

    return list(per_gene.items())


def read_spliceai_region(vcf_gz: str, region: str) -> pd.DataFrame:
    """
    Stream a single region/contig from a tabix-indexed SpliceAI VCF.

    Returns a DataFrame with columns:
      - ID (chr:pos:ref:alt)
      - SP_GENE
      - SpliceAI_max (per (ID, gene) collapsed)

    Note: This function expects the SpliceAI payload to be in INFO as "SpliceAI=...".
    """
    cmd = f"tabix {shlex.quote(vcf_gz)} {shlex.quote(region)}"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True, errors="replace")
    if proc.stdout is None:
        die(f"Failed to read {vcf_gz} {region} with tabix.")

    ids: List[str] = []
    genes: List[str] = []
    dsvals: List[float] = []

    for line in proc.stdout:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            continue

        chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
        chrom = chrom if chrom.startswith("chr") else ("chr" + chrom)
        var_id = f"{chrom}:{pos}:{ref}:{alt}"

        info_field = fields[7]
        m = _RE_SPLICE_PAYLOAD.search(info_field)
        if not m:
            continue

        for g, ds in parse_spliceai_info(m.group(1)):
            if isinstance(g, str) and g.startswith("ENSG"):
                g = g.split(".")[0]
            ids.append(var_id)
            genes.append(g)
            dsvals.append(ds)

    proc.stdout.close()
    proc.wait()

    if not ids:
        return pd.DataFrame(columns=["ID", "SP_GENE", "SpliceAI_max"])

    sp = pd.DataFrame({"ID": ids, "SP_GENE": genes, "SpliceAI_max": dsvals})
    sp["SpliceAI_max"] = pd.to_numeric(sp["SpliceAI_max"], errors="coerce")
    sp = sp.groupby(["ID", "SP_GENE"], as_index=False, sort=False)["SpliceAI_max"].max()
    return sp


def read_spliceai_all_by_contig(vcf_gz: str, n_jobs: int = 1) -> pd.DataFrame:
    """
    Parse all contigs in a tabix-indexed SpliceAI VCF, optionally parallelized by contig.
    Returns DataFrame ['ID','SP_GENE','SpliceAI_max'] with max per (ID, gene).
    """
    contigs = list_contigs(vcf_gz)
    eprint(f"[INFO] SpliceAI contigs in {Path(vcf_gz).name}: {len(contigs)}")

    if n_jobs > 1:
        workers = min(max(1, n_jobs), cpu_count())
        eprint(f"[INFO] Parallel SpliceAI parsing: {workers} workers")
        with Pool(processes=workers) as pool:
            parts = pool.starmap(read_spliceai_region, [(vcf_gz, c) for c in contigs])
    else:
        parts = []
        for c in contigs:
            eprint(f"[INFO] Parsing contig {c}")
            parts.append(read_spliceai_region(vcf_gz, c))

    if not parts:
        return pd.DataFrame(columns=["ID", "SP_GENE", "SpliceAI_max"])

    sp = pd.concat(parts, ignore_index=True)
    sp = sp.groupby(["ID", "SP_GENE"], as_index=False, sort=False)["SpliceAI_max"].max()
    return sp


# ==============================================================================
# BRaVa-like categorization
# ==============================================================================

PLOF_CSQS = {
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
}
MISSENSE_CSQS = {
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
}
SYNONYMOUS_CSQS = {"stop_retained_variant", "synonymous_variant"}
OTHER_CSQS = {
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "feature_elongation",
    "regulatory_region_variant",
    "feature_truncation",
    "intergenic_variant",
}
INFRAME_CSQS = {"inframe_deletion", "inframe_insertion"}


def classify_variant(
    row: pd.Series,
    revel_cut: float,
    cadd_cut: float,
    spliceai_cut: float,
    col_lof: str,
    col_revel: str,
    col_cadd: str,
    col_csq: str,
) -> str:
    """
    Classify a variant into BRaVa-like categories based on VEP + LoFTEE + scores.
    """
    csq = str(row.get(col_csq, "")) if pd.notna(row.get(col_csq, np.nan)) else ""
    consequences = set(csq.split("&")) if csq else set()

    lof = str(row.get(col_lof, "")) if pd.notna(row.get(col_lof, np.nan)) else ""
    is_plof = lof == "HC"
    is_lof_low = lof == "LC"

    missense = bool(consequences & MISSENSE_CSQS)
    synonymous = bool(consequences & SYNONYMOUS_CSQS)
    other = bool(consequences & OTHER_CSQS)
    inframe = bool(consequences & INFRAME_CSQS)

    revel = row.get(col_revel, np.nan)
    cadd = row.get(col_cadd, np.nan)
    spmax = row.get("SpliceAI_max", np.nan)

    if is_plof:
        return "pLoF"

    if missense and (
        (pd.notna(revel) and float(revel) >= revel_cut) or (pd.notna(cadd) and float(cadd) >= cadd_cut)
    ):
        return "damaging_missense_or_protein_altering"

    if pd.notna(spmax) and float(spmax) >= spliceai_cut:
        return "damaging_missense_or_protein_altering"

    if is_lof_low:
        return "damaging_missense_or_protein_altering"

    if missense or inframe:
        return "other_missense_or_protein_altering"

    if other:
        return "non_coding"

    if synonymous:
        return "synonymous"

    return ""


# ==============================================================================
# CLI + main
# ==============================================================================

@dataclass(frozen=True)
class ColNames:
    snpid: str
    gene: str
    lof: str
    revel: str
    cadd: str
    csq: str
    canonical: str
    biotype: str
    mane: str


def parse_revel(cell: object) -> float:
    """
    REVEL may be '0.12&0.34' -> take first parseable numeric component.
    """
    if cell is None:
        return np.nan
    s = str(cell)
    if s in {"", "."}:
        return np.nan
    for p in s.split("&"):
        if p and p != ".":
            try:
                return float(p)
            except ValueError:
                continue
    return np.nan


def to_float(cell: object) -> float:
    try:
        s = str(cell)
        if s in {"", "."}:
            return np.nan
        return float(s)
    except Exception:
        return np.nan


def read_vep_table(vep_path: str, cn: ColNames) -> pd.DataFrame:
    """
    Read the VEP table (whitespace-delimited), keep core columns, and apply filters:
      - BIOTYPE == protein_coding
      - MANE_SELECT not-null OR (MANE_SELECT null AND CANONICAL == YES)
    """
    usecols = [cn.snpid, cn.gene, cn.lof, cn.revel, cn.cadd, cn.csq, cn.canonical, cn.biotype, cn.mane]

    eprint("[INFO] Reading VEP table ...")
    try:
        vep = pd.read_csv(
            vep_path,
            sep=r"\s+",
            usecols=usecols,
            na_values=".",
            dtype=str,
            engine="python",
            encoding="cp1252",
        )
    except Exception as ex:
        die(f"Failed reading VEP table '{vep_path}': {ex}")

    eprint(f"[INFO] VEP rows (raw): {len(vep):,}")

    vep["SNP_ID_final"] = vep[cn.snpid].astype(str)
    vep["SNP_ID_merge"] = vep["SNP_ID_final"].map(extract_tail_snpid)

    keep = (vep[cn.biotype] == "protein_coding") & (
        vep[cn.mane].notna() | (vep[cn.mane].isna() & (vep[cn.canonical] == "YES"))
    )
    before = len(vep)
    vep = vep.loc[keep].copy()
    eprint(f"[INFO] Filter protein_coding + MANE/CANONICAL: kept {len(vep):,} / {before:,}")

    eprint("[INFO] Parsing REVEL + CADD ...")
    vep[cn.revel] = vep[cn.revel].map(parse_revel)
    vep[cn.cadd] = vep[cn.cadd].map(to_float)

    return vep


def read_cadd_indels(cadd_path: str) -> pd.DataFrame:
    """
    Read a combined CADD indel file with columns:
      Chrom  Pos  Ref  Alt  RawScore  PHRED
    Returns DataFrame: snpid (chr:pos:ref:alt), PHRED_indel (float).
    """
    eprint(f"[INFO] Reading CADD indels: {cadd_path}")
    try:
        cadd = pd.read_csv(
            cadd_path,
            sep=r"\t",
            comment="#",
            header=None,
            names=["Chrom", "Pos", "Ref", "Alt", "RawScore", "PHRED"],
            dtype=str,
            engine="python",
        )
    except Exception as ex:
        die(f"Failed reading CADD indels '{cadd_path}': {ex}")

    if cadd.empty:
        return pd.DataFrame(columns=["snpid", "PHRED_indel"])

    cadd["Chrom"] = cadd["Chrom"].astype(str).str.replace(r"^chr", "", regex=True)
    cadd["snpid"] = "chr" + cadd["Chrom"] + ":" + cadd["Pos"].astype(str) + ":" + cadd["Ref"] + ":" + cadd["Alt"]
    cadd["PHRED_indel"] = pd.to_numeric(cadd["PHRED"], errors="coerce")
    cadd = cadd[["snpid", "PHRED_indel"]].dropna().drop_duplicates()
    return cadd


def build_spliceai_table(
    snv_vcf: str,
    indel_vcf: str,
    contig: Optional[str],
    n_jobs: int,
) -> pd.DataFrame:
    """
    Build a combined SpliceAI table from SNV + INDEL SpliceAI VCFs.
    Output columns: ID, SP_GENE, SpliceAI_max
    """
    eprint("[INFO] Reading SpliceAI (tabix streaming) ...")

    def load_one(vcf_gz: str, label: str) -> pd.DataFrame:
        if contig:
            contigs = list_contigs(vcf_gz)
            region = normalize_region_arg(contigs, contig)
            eprint(f"[INFO] {label}: using contig {region}")
            return read_spliceai_region(vcf_gz, region)
        return read_spliceai_all_by_contig(vcf_gz, n_jobs=n_jobs)

    sp_snv = load_one(snv_vcf, "SNV")
    sp_indel = load_one(indel_vcf, "INDEL")

    sp = pd.concat([sp_snv, sp_indel], ignore_index=True)
    if sp.empty:
        die("Parsed 0 SpliceAI rows from SNV+INDEL. Check tabix index and INFO parsing.")
    sp = sp.groupby(["ID", "SP_GENE"], as_index=False, sort=False)["SpliceAI_max"].max()
    eprint(f"[INFO] SpliceAI rows (per variant x gene): {len(sp):,}")
    return sp


def write_saige_group_file(df: pd.DataFrame, gene_col: str, out_path: str) -> None:
    """
    Write SAIGE group file in the "gene var" + "gene anno" two-line style.
    This script uses SNP_ID_final (original VEP ID) to preserve DRAGEN-style IDs if present.
    """
    eprint(f"[INFO] Writing SAIGE group file: {out_path}")
    with open(out_path, "w") as fout:
        for gene, sub in df.groupby(gene_col, sort=False):
            snps = " ".join(sub["SNP_ID_final"].astype(str).tolist())
            annos = " ".join(sub["annotation"].astype(str).tolist())
            fout.write(f"{gene} var {snps}\n")
            fout.write(f"{gene} anno {annos}\n")


def main() -> None:
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Merge VEP + SpliceAI (Illumina precomputed SNV+INDEL) -> long table + SAIGE group file",
    )

    ap.add_argument("--vep", "-v", required=True, help="VEP table (whitespace-delimited).")
    ap.add_argument("--spliceai-snv", required=True, help="Illumina precomputed SNV SpliceAI VCF (.vcf.gz)")
    ap.add_argument("--spliceai-indel", required=True, help="Illumina precomputed INDEL SpliceAI VCF (.vcf.gz)")
    ap.add_argument("--out", "-o", required=True, help="Output prefix. Writes <out>.long.csv.gz and <out> (group file).")
    ap.add_argument("--cadd-indels", default=None, help="Combined CADD indel file (.tsv or .tsv.gz).")

    # Performance / region control
    ap.add_argument("--contig", default=None, help="Contig to parse only (e.g. 22 or chr22).")
    ap.add_argument("--n-jobs", type=int, default=1, help="Parallel workers for scanning all contigs.")

    # VEP columns
    ap.add_argument("--col-snpid", default="SNP_ID")
    ap.add_argument("--col-gene", default="GENE")
    ap.add_argument("--col-lof", default="LOF")
    ap.add_argument("--col-revel", default="REVEL_SCORE")
    ap.add_argument("--col-cadd", default="CADD_PHRED")
    ap.add_argument("--col-csq", default="CSQ")
    ap.add_argument("--col-canonical", default="CANONICAL")
    ap.add_argument("--col-biotype", default="BIOTYPE")
    ap.add_argument("--col-mane", default="MANE_SELECT")

    # Cutoffs
    ap.add_argument("--revel-cut", type=float, default=0.773)
    ap.add_argument("--cadd-cut", type=float, default=28.1)
    ap.add_argument("--spliceai-cut", type=float, default=0.20)

    args = ap.parse_args()

    # Basic file checks
    check_exists(args.vep, "file")
    check_exists(args.spliceai_snv, "file")
    check_exists(args.spliceai_indel, "file")
    if args.cadd_indels:
        check_exists(args.cadd_indels, "file")

    # Column map
    cn = ColNames(
        snpid=args.col_snpid,
        gene=args.col_gene,
        lof=args.col_lof,
        revel=args.col_revel,
        cadd=args.col_cadd,
        csq=args.col_csq,
        canonical=args.col_canonical,
        biotype=args.col_biotype,
        mane=args.col_mane,
    )

    # 1) VEP
    vep = read_vep_table(args.vep, cn)

    # 2) SpliceAI
    sp = build_spliceai_table(
        snv_vcf=args.spliceai_snv,
        indel_vcf=args.spliceai_indel,
        contig=args.contig,
        n_jobs=max(1, args.n_jobs),
    )

    # 3) Merge: per-variant DSmax across genes
    eprint("[INFO] Merging VEP with SpliceAI on variant ID (per-variant DSmax across genes) ...")
    per_variant = sp.groupby("ID", as_index=False, sort=False)["SpliceAI_max"].max()
    merged = vep.merge(per_variant, left_on="SNP_ID_merge", right_on="ID", how="left")

    n_sp = int(merged["SpliceAI_max"].notna().sum())
    eprint(f"[CHECK] SpliceAI_max present: {n_sp:,} / {len(merged):,}")
    if n_sp == 0:
        die("SpliceAI coverage is 0 after merge. Check ID normalization and contig naming.")

    # 4) Optional: CADD indels patch-in where CADD is missing
    if args.cadd_indels:
        before = int(merged[cn.cadd].notna().sum())
        cadd = read_cadd_indels(args.cadd_indels)
        if not cadd.empty:
            merged = merged.merge(cadd, left_on="SNP_ID_merge", right_on="snpid", how="left")
            merged[cn.cadd] = merged[cn.cadd].where(
                merged[cn.cadd].notna(),
                merged["PHRED_indel"],
            )
            merged.drop(columns=["snpid", "PHRED_indel"], inplace=True, errors="ignore")

        merged[cn.cadd] = merged[cn.cadd].map(to_float)
        after = int(merged[cn.cadd].notna().sum())
        eprint(f"[INFO] CADD PHRED non-NA: before {before:,} -> after {after:,}")

    # 5) Classify
    eprint("[INFO] Classifying variants ...")
    merged["SpliceAI_max"] = pd.to_numeric(merged["SpliceAI_max"], errors="coerce")

    merged["annotation"] = merged.apply(
        lambda r: classify_variant(
            r,
            revel_cut=args.revel_cut,
            cadd_cut=args.cadd_cut,
            spliceai_cut=args.spliceai_cut,
            col_lof=cn.lof,
            col_revel=cn.revel,
            col_cadd=cn.cadd,
            col_csq=cn.csq,
        ),
        axis=1,
    )

    kept = merged.loc[merged["annotation"] != ""].copy()
    eprint(f"[INFO] Kept (non-empty annotation): {len(kept):,} / {len(merged):,}")
    eprint("[INFO] Annotation counts:\n" + kept["annotation"].value_counts().to_string())

    # 6) Outputs
    out_long = args.out + ".long.csv.gz"
    out_group = args.out  # keep your convention

    cols_out = [
        "SNP_ID_final",
        "SNP_ID_merge",
        cn.gene,
        cn.lof,
        cn.revel,
        cn.cadd,
        cn.csq,
        "SpliceAI_max",
        "annotation",
    ]
    cols_out = [c for c in cols_out if c in kept.columns]

    eprint(f"[INFO] Writing long table: {out_long}")
    kept[cols_out].to_csv(out_long, index=False, compression="gzip")

    write_saige_group_file(kept, gene_col=cn.gene, out_path=out_group)

    eprint("[INFO] Done.")


if __name__ == "__main__":
    main()