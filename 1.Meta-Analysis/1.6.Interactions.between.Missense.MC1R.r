############################################################
## MC1R missense variants × circulating KIT (interaction)
## - Additive + interaction regression (covariate-adjusted)
## - Descriptive 4-group carrier mean comparison + one-way ANOVA
## Author: Shelley
## Date:   2025-10-31
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(broom)
  library(ggplot2)
})

## -----------------------------
## 0) Paths (edit these as needed)
## -----------------------------y
kit_file  <- "MC1R.KIT.merged.tsv"
meta_file <- "MC1R.meta.info.tsv"

stopifnot(file.exists(kit_file), file.exists(meta_file))

## -----------------------------
## 1) Read data
## -----------------------------
kit_dt  <- fread(kit_file)   # KIT + genotypes (one row per sample)
meta_dt <- fread(meta_file)  # covariates

cat(sprintf("[INFO] kit_dt:  %d rows x %d cols\n", nrow(kit_dt),  ncol(kit_dt)))
cat(sprintf("[INFO] meta_dt: %d rows x %d cols\n", nrow(meta_dt), ncol(meta_dt)))

stopifnot("eid" %in% names(kit_dt), "eid" %in% names(meta_dt))

## -----------------------------
## 2) Merge by eid (left join; keep KIT/genotype rows)
## -----------------------------
merged_dt <- merge(
  kit_dt,
  meta_dt,
  by = "eid",
  all.x = TRUE,
  sort = FALSE
)
cat(sprintf("[INFO] merged_dt: %d rows x %d cols\n", nrow(merged_dt), ncol(merged_dt)))

## -----------------------------
## 3) Normalize key covariate names
## -----------------------------
setnames(merged_dt, old = "p74_i0",  new = "fasting_time", skip_absent = TRUE)
setnames(merged_dt, old = "p21022",  new = "age",          skip_absent = TRUE)
setnames(merged_dt, old = "p31",     new = "sex",          skip_absent = TRUE)

# PCs: p22009_a1 ... p22009_a40 -> PC1 ... PC40
pc_old <- sprintf("p22009_a%d", 1:40)
pc_new <- sprintf("PC%d",        1:40)
present <- pc_old %in% names(merged_dt)
setnames(merged_dt, old = pc_old[present], new = pc_new[present], skip_absent = TRUE)

## -----------------------------
## 4) Define genotype columns (EDIT if column names differ)
## -----------------------------
stopifnot("kit_irnt" %in% names(merged_dt))  # requires precomputed IRNT KIT

v1_raw <- "DRAGEN:chr16:89919510:C:A"  # v1
v2_raw <- "DRAGEN:chr16:89919709:C:T"  # v2
stopifnot(v1_raw %in% names(merged_dt), v2_raw %in% names(merged_dt))

merged_dt[, mc1r_v1 := as.numeric(get(v1_raw))]
merged_dt[, mc1r_v2 := as.numeric(get(v2_raw))]

## -----------------------------
## 5) Build covariates + analysis dataset (complete cases)
## -----------------------------
covars <- c(
  paste0("PC", 1:20),
  "age",
  "sex",
  "I(age^2)",
  "sex:age",
  "sex:I(age^2)",
  "fasting_time"
)
covar_str <- paste(covars, collapse = " + ")

# make sex explicit (avoid 1/2 numeric pitfalls)
merged_dt[, sex := as.factor(sex)]

need_cols <- c(
  "kit_irnt", "mc1r_v1", "mc1r_v2",
  paste0("PC", 1:20),
  "age", "sex", "fasting_time"
)

ana_dt <- merged_dt[complete.cases(merged_dt[, ..need_cols])]
cat(sprintf("[INFO] analysis N (complete cases): %d\n", nrow(ana_dt)))

## -----------------------------
## 6) Linear models (covariate-adjusted)
## -----------------------------
f1 <- as.formula(paste0("kit_irnt ~ mc1r_v1 + ", covar_str))
f2 <- as.formula(paste0("kit_irnt ~ mc1r_v2 + ", covar_str))
f3 <- as.formula(paste0("kit_irnt ~ mc1r_v1 + mc1r_v2 + ", covar_str))
f4 <- as.formula(paste0("kit_irnt ~ mc1r_v1 * mc1r_v2 + ", covar_str))

m1 <- lm(f1, data = ana_dt)
m2 <- lm(f2, data = ana_dt)
m3 <- lm(f3, data = ana_dt)
m_int <- lm(f4, data = ana_dt)

## -----------------------------
## 7) Summarize MC1R terms across models
## -----------------------------
mods <- list(
  model1_v1 = m1,
  model2_v2 = m2,
  model3_both_add = m3,
  model4_interaction = m_int
)

tidy_dt <- rbindlist(lapply(names(mods), function(nm) {
  tt <- broom::tidy(mods[[nm]])
  tt$model <- nm
  tt
}), use.names = TRUE, fill = TRUE)

keep_terms <- c("mc1r_v1", "mc1r_v2", "mc1r_v1:mc1r_v2")
mc1r_dt <- tidy_dt[term %in% keep_terms]

mc1r_dt[, term := factor(term, levels = keep_terms)]
mc1r_dt[, model := factor(model,
                          levels = c("model1_v1","model2_v2","model3_both_add","model4_interaction"))]
setorder(mc1r_dt, model, term)

summary_tab <- mc1r_dt[, .(model, term, estimate, std.error, statistic, p.value)]
print(summary_tab)

fwrite(summary_tab,
       file = "MC1R_KIT.interaction.models.summary.txt",
       sep = "\t", quote = FALSE, na = "NA")
cat("[WRITE] MC1R_KIT.interaction.models.summary.txt\n")

## -----------------------------
## 8) Descriptive 4-group carrier analysis (UNADJUSTED)
## -----------------------------
dt <- copy(ana_dt)

dt[, v1_car := fifelse(mc1r_v1 > 0, 1L, 0L)]
dt[, v2_car := fifelse(mc1r_v2 > 0, 1L, 0L)]

# Explicit mapping avoids v1/v2 swap bugs
dt[, geno4 := fcase(
  v1_car == 0L & v2_car == 0L, "none",
  v1_car == 0L & v2_car == 1L, "v2 only",
  v1_car == 1L & v2_car == 0L, "v1 only",
  v1_car == 1L & v2_car == 1L, "both",
  default = NA_character_
)]
dt[, geno4 := factor(geno4, levels = c("none", "v2 only", "v1 only", "both"))]

cat("[CHECK] group counts:\n")
print(dt[!is.na(geno4), .N, by = geno4])

dt_plot <- dt[!is.na(geno4) & !is.na(kit_irnt)]

# One-way ANOVA (unadjusted; descriptive)
fit_aov  <- aov(kit_irnt ~ geno4, data = dt_plot)
p_global <- summary(fit_aov)[[1]][["Pr(>F)"]][1]
p_lab    <- paste0("One-way ANOVA (unadjusted) P = ", formatC(p_global, format = "e", digits = 2))

# group summary
sum_dt <- dt_plot[, .(
  mean_kit = mean(kit_irnt, na.rm = TRUE),
  sd_kit   = sd(kit_irnt,   na.rm = TRUE),
  n        = .N
), by = geno4]

sum_dt[, se   := sd_kit / sqrt(n)]
sum_dt[, ci_l := mean_kit - 1.96 * se]
sum_dt[, ci_u := mean_kit + 1.96 * se]

# expected additive for "both"
mean_none <- sum_dt[geno4 == "none",    mean_kit]
mean_v2   <- sum_dt[geno4 == "v2 only", mean_kit]
mean_v1   <- sum_dt[geno4 == "v1 only", mean_kit]
expected_both <- mean_v1 + mean_v2 - mean_none

expected_df <- data.table(
  geno4    = factor("both", levels = levels(sum_dt$geno4)),
  expected = expected_both
)

## -----------------------------
## 9) Simple plot (GitHub-friendly; you will make Nature theme later)
## -----------------------------
col_map <- c(
  "none"    = "#2F3B4C",
  "v2 only" = "#4A7DB2",
  "v1 only" = "#4A7DB2",
  "both"    = "#5B4B8A"
)

p <- ggplot(sum_dt, aes(x = geno4, y = mean_kit)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey80") +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u, color = geno4),
                width = 0.16, linewidth = 0.75, show.legend = FALSE) +
  geom_point(aes(color = geno4), size = 1.7, show.legend = FALSE) +
  geom_point(data = expected_df, aes(x = geno4, y = expected),
             inherit.aes = FALSE, shape = 4, stroke = 0.8, size = 2.0, color = "grey35") +
  geom_segment(inherit.aes = FALSE,
               aes(x = 3.7, xend = 4.3, y = expected_both, yend = expected_both),
               color = "grey35", linetype = "dotted", linewidth = 0.55) +
  scale_color_manual(values = col_map) +
  labs(x = NULL, y = "KIT (IRNT)",
       title = "MC1R allelic combinations and circulating KIT",
       subtitle = p_lab) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(size = 9))

ggsave("MC1R_KIT_interaction_plot.github.pdf",
       p, width = 3.8, height = 3.0, units = "in",
       device = cairo_pdf, dpi = 300, useDingbats = FALSE)

fwrite(sum_dt,
       file = "MC1R_KIT_genotype_group_summary.tsv",
       sep = "\t", quote = FALSE, na = "NA")

cat("[WRITE] MC1R_KIT_interaction_plot.github.pdf\n")
cat("[WRITE] MC1R_KIT_genotype_group_summary.tsv\n")