#!/usr/bin/env bash
# pipeline to run rare variant association analysis using SAIGE on UKB-RAP (DNAnexus)
set -euo pipefail

###############################################################################
# Step 1) Fit NULL model
###############################################################################

cmd="step1_fitNULLGLMM.R \
    --plinkFile=PLINK_for_vr_melanoma_whole_cohort \
    --phenoFile=Input_phenotype_for_nullmodel_all_melanoma_cohort.WGS.version.txt \
    --phenoCol=melanoma_status \
    --traitType=binary \
    --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,age_at_recruitment,genetic_sex,age_squared,sex_age,sex_age_squared,wgs_batch \
    --qCovarColList=genetic_sex,wgs_batch \
    --sampleIDColinphenoFile=eid \
    --nThreads=4 \
    --LOCO=FALSE \
    --outputPrefix=WGS_RV_melanoma_whole_cohort_null_model \
    --sparseGRMFile=saige_gene_step0_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMSampleIDFile=saige_gene_step0_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
    --isCateVarianceRatio=TRUE \
    --useSparseGRMtoFitNULL=TRUE \
    --useSparseGRMforVarRatio=TRUE \
    --SampleIDIncludeFile=WGS_shared_sample_ids.txt"

# ---- Inputs (hide project ID via env var) ----
# Set once:
#   export DX_PROJECT="project-xxxxxxxxxxxxxxxxxxxxxxxx"
: "${DX_PROJECT:?Set DX_PROJECT first (e.g., export DX_PROJECT='project-xxxx')}"

bedFile="${DX_PROJECT}:/2.RV_collapsing/Step1_fitnull_model/PLINK_for_vr_melanoma_whole_cohort.bed"
bimFile="${DX_PROJECT}:/2.RV_collapsing/Step1_fitnull_model/PLINK_for_vr_melanoma_whole_cohort.bim"
famFile="${DX_PROJECT}:/2.RV_collapsing/Step1_fitnull_model/PLINK_for_vr_melanoma_whole_cohort.fam"

subsampleFile="${DX_PROJECT}:/3.WGS_noncoding/WGS_shared_sample_ids.txt"

sparseGRM_matrix="${DX_PROJECT}:/1.GWAS-imputation/Step0.SparseGRM/sparseGRM/saige_gene_step0_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
sparseGRM_sample_file="${DX_PROJECT}:/1.GWAS-imputation/Step0.SparseGRM/sparseGRM/saige_gene_step0_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"

phenotype_file="${DX_PROJECT}:/3.WGS_noncoding/Input_phenotype_for_nullmodel_all_melanoma_cohort.WGS.version.txt"

dx run app-swiss-army-knife \
    -iin="${bedFile}" \
    -iin="${bimFile}" \
    -iin="${famFile}" \
    -iin="${sparseGRM_matrix}" \
    -iin="${sparseGRM_sample_file}" \
    -iin="${phenotype_file}" \
    -iin="${subsampleFile}" \
    -icmd="${cmd}" \
    -iimage_file="${DX_PROJECT}:/docker_images/saige-latest-new-1.4.0.tar.gz" \
    --folder="${DX_PROJECT}:/3.WGS_noncoding/Step1_fitnull_model/" \
    --instance-type=mem1_ssd2_v2_x4 \
    --name=WGS_step1_fitnull_model \
    --priority=high \
    -y

###############################################################################
# Step 2) Association tests (example: melanocyte cCREs)
###############################################################################

chr="${1:-}"
if [[ -z "${chr}" ]]; then
  echo "Usage: $0 <chr>" >&2
  echo "Example: $0 22" >&2
  exit 1
fi

cmd="step2_SPAtests.R \
    --bgenFile=/mnt/project/3.WGS_noncoding/WGS_bgen.1.2/WGS.chr${chr}.bgen.1.2version.QCed.bgen \
    --bgenFileIndex=/mnt/project/3.WGS_noncoding/WGS_bgen.1.2/WGS.chr${chr}.bgen.1.2version.QCed.bgen.bgi \
    --sampleFile=/mnt/project/3.WGS_noncoding/WGS_bgen.1.2/WGS.chr${chr}.bgen.1.2version.QCed.sample \
    --AlleleOrder=ref-first \
    --SAIGEOutputFile=noncoding.RV.melanocyte.cCRE.chr${chr}.txt \
    --chrom=${chr} \
    --LOCO=FALSE \
    --GMMATmodelFile=WGS_RV_melanoma_whole_cohort_null_model.rda \
    --varianceRatioFile=WGS_RV_melanoma_whole_cohort_null_model.varianceRatio.txt \
    --sparseGRMFile=saige_gene_step0_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMSampleIDFile=saige_gene_step0_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
    --groupFile=combined_mask.mel.cCRE.chr${chr}.txt \
    --annotation_in_groupTest=CADD:CADD_GERP:CADD_JARVIS:CADD_GERP_JARVIS,GERP:GERP_JARVIS:CADD_GERP:CADD_GERP_JARVIS,JARVIS:CADD_JARVIS:GERP_JARVIS:CADD_GERP_JARVIS,CADD:CADD_GERP:CADD_JARVIS:CADD_GERP_JARVIS:JARVIS:GERP_JARVIS,CADD:CADD_GERP:CADD_JARVIS:CADD_GERP_JARVIS:GERP:GERP_JARVIS,CADD:CADD_GERP:CADD_JARVIS:CADD_GERP_JARVIS:JARVIS:GERP_JARVIS:GERP,CADD:CADD_GERP:CADD_JARVIS:CADD_GERP_JARVIS:JARVIS:GERP_JARVIS:GERP:Null \
    --maxMAF_in_groupTest=0.0001,0.001,0.01 \
    --is_output_markerList_in_groupTest=TRUE \
    --subSampleFile=WGS_shared_sample_ids.txt \
    --is_single_in_groupTest=TRUE \
    --is_output_moreDetails=TRUE"

echo "${cmd}"

subsampleFile="${DX_PROJECT}:/3.WGS_noncoding/WGS_shared_sample_ids.txt"
step1_model_file="${DX_PROJECT}:/3.WGS_noncoding/Step1_fitnull_model/WGS_RV_melanoma_whole_cohort_null_model.rda"
step1_varianceratio="${DX_PROJECT}:/3.WGS_noncoding/Step1_fitnull_model/WGS_RV_melanoma_whole_cohort_null_model.varianceRatio.txt"
sparseGRM_matrix="${DX_PROJECT}:/1.GWAS-imputation/Step0.SparseGRM/sparseGRM/saige_gene_step0_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
sparseGRM_sample_file="${DX_PROJECT}:/1.GWAS-imputation/Step0.SparseGRM/sparseGRM/saige_gene_step0_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
target_groupfile="${DX_PROJECT}:/3.WGS_noncoding/melanocyte_cCRE_masks/combined_mask.mel.cCRE.chr${chr}.txt"
output_folder="${DX_PROJECT}:/3.WGS_noncoding/SAIGE.Results/melanocyte.cCRE/chr${chr}"

dx run app-swiss-army-knife \
    -iin="${subsampleFile}" \
    -iin="${step1_model_file}" \
    -iin="${step1_varianceratio}" \
    -iin="${sparseGRM_matrix}" \
    -iin="${sparseGRM_sample_file}" \
    -iin="${target_groupfile}" \
    -icmd="${cmd}" \
    -iimage_file="${DX_PROJECT}:/docker_images/saige.1.4.0.tar.gz" \
    --folder="${output_folder}" \
    --instance-type=mem2_ssd2_v2_x4 \
    --name=melanoma.cCRE.chr${chr} \
    --priority=high \
    -y \
    --allow-ssh

# usage : bash WGS.noncoding.SAIGE.cCRE.pipeline.sh <chr>