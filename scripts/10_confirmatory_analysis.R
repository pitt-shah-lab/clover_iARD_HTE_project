# ============================================================================
# 10_confirmatory_analysis.R
#
# PURPOSE (Victor's stage 6): apply the best CATE model from derivation
# (stage 5) to the held-out validation set. Classify validation patients
# into benefit / indeterminate / harm subgroups based on derivation-derived
# thresholds, and test whether the treatment effect differs across subgroups.
#
# This is the step that confirms whether the HTE signal found in derivation
# replicates out of sample, or was a chance discovery.
#
# Uses confirm_analysis() from cate-repo/confirm_analysis.R.
#
# NOTE: confirm_analysis.R has a browser() call on line 74 that pauses
# execution in interactive R sessions. In non-interactive Rscript mode,
# browser() is silently skipped. If running interactively, type 'c' to
# continue when prompted.
#
# Paths come from config/config.R. Inputs from outputs/04_impute/ and
# outputs/03_split_derivation_validation/.
#
# Output (written to outputs/10_confirmatory_analysis/):
#   confirmatory_preds.csv         -- CATE predictions on validation
#   confirmatory_subgroups.csv     -- benefit/indeterminate/harm table
#   confirmatory_importance.csv    -- variable importance from the final model
#   confirmatory_report.txt        -- full console output
# ============================================================================

source(file.path((getwd()), "config", "config.R"))
out_dir <- make_output_subdir("10_confirmatory_analysis")

report_path <- file.path(out_dir, "confirmatory_report.txt")
report_con  <- file(report_path, open = "wt")
sink(report_con, split = TRUE)

# ---- source Victor's tools ----
source(SCORING_METHODS_R)
source(CROSS_VALIDATION_R)
source(POST_PROCESS_R)
source(CONFIRM_ANALYSIS_R)

# ---- load imputed data + the split ----
dsi_all <- readRDS(file.path(OUTPUTS_DIR, "04_impute", "dsi_all.RDS"))
ids     <- readRDS(file.path(OUTPUTS_DIR, "03_split_derivation_validation",
                              "ids_internal_list.RDS"))

der <- dsi_all[dsi_all$id %in% ids$ids_der, ]
val <- dsi_all[dsi_all$id %in% ids$ids_val, ]
rownames(der) <- 1:nrow(der)
rownames(val) <- 1:nrow(val)

cat("Derivation set:", nrow(der), "patients\n")
cat("Validation set:", nrow(val), "patients\n\n")

# ---- covariate lists (same as 07_cate_modeling.R) ----
clinical_norm <- c("age","temp","rr","hr","map","sbp",
                    "sofa","albumin","ln_bili",
                    "cr","bun","ln_g","ln_lac","sqrt_plt",
                    "ln_wbc","hgb","na","bicarb","prefluid",
                    "ln_bmi","gcs","charlson","o2sat","s2f")
clinical_bin  <- c("site_lung","site_abdom","site_urine","mv",
                    "vaso","ards","dial","chf","copd",
                    "liver","kidney")
biomarker_vars <- grep("^ln_(il1|ang1|ang2|tnfr|il6|strem1|kim1|srage)",
                        names(der), value = TRUE)

clinical_only_x  <- c(clinical_norm, clinical_bin)
with_biomarker_x <- c(clinical_norm, clinical_bin, biomarker_vars)

# ============================================================================
# STEP 1: Get CATE predictions on validation
#
# confirm_analysis with out="preds" trains the CATE model on derivation
# and predicts on validation. We run this for the best method (cf.CATE)
# on BOTH covariate sets (clinical only and clinical + biomarker).
# ============================================================================

cat("============================================================\n")
cat("STEP 1: CATE predictions on validation set\n")
cat("============================================================\n\n")

train_list <- list(der)
test_list  <- list(val)

# --- Clinical only ---
cat("Running cf.CATE on validation (clinical only, 35 covariates)...\n")
t0 <- Sys.time()
preds_clinical <- confirm_analysis(
  train_list = train_list,
  test_list  = test_list,
  mod_ben    = "cf.CATE",
  mod_harm   = "cf.CATE",
  cval_ben   = 0,     # placeholder, not used when out="preds"
  cval_harm  = 0,
  y       = "inhosp90",
  w       = "w",
  xvars   = clinical_only_x,
  pi      = 0.5,
  type    = "CATE",
  out     = "preds"
)
cat("  Runtime:", format(round(difftime(Sys.time(), t0), 1)), "\n")

preds_clin_df  <- preds_clinical[[1]]
import_clin    <- preds_clinical[[2]]

# --- Clinical + biomarker ---
cat("Running cf.CATE on validation (clinical + biomarker, 43 covariates)...\n")
t0 <- Sys.time()
preds_bio <- confirm_analysis(
  train_list = train_list,
  test_list  = test_list,
  mod_ben    = "cf.CATE",
  mod_harm   = "cf.CATE",
  cval_ben   = 0,
  cval_harm  = 0,
  y       = "inhosp90",
  w       = "w",
  xvars   = with_biomarker_x,
  pi      = 0.5,
  type    = "CATE",
  out     = "preds"
)
cat("  Runtime:", format(round(difftime(Sys.time(), t0), 1)), "\n\n")

preds_bio_df  <- preds_bio[[1]]
import_bio    <- preds_bio[[2]]

# ---- Save predictions ----
write.csv(preds_clin_df,
          file.path(out_dir, "confirmatory_preds_clinical.csv"),
          row.names = FALSE)
write.csv(preds_bio_df,
          file.path(out_dir, "confirmatory_preds_with_biomarker.csv"),
          row.names = FALSE)

# ---- Save importance ----
if (!is.null(import_clin)) {
  imp_clin_df <- data.frame(
    variable = clinical_only_x,
    importance = as.numeric(import_clin),
    stringsAsFactors = FALSE
  )
  imp_clin_df <- imp_clin_df[order(-imp_clin_df$importance), ]
  write.csv(imp_clin_df,
            file.path(out_dir, "confirmatory_importance_clinical.csv"),
            row.names = FALSE)
  cat("Top 10 variables (clinical only):\n")
  print(head(imp_clin_df, 10), row.names = FALSE)
  cat("\n")
}

if (!is.null(import_bio)) {
  imp_bio_df <- data.frame(
    variable = with_biomarker_x,
    importance = as.numeric(import_bio),
    stringsAsFactors = FALSE
  )
  imp_bio_df <- imp_bio_df[order(-imp_bio_df$importance), ]
  write.csv(imp_bio_df,
            file.path(out_dir, "confirmatory_importance_with_biomarker.csv"),
            row.names = FALSE)
  cat("Top 10 variables (clinical + biomarker):\n")
  print(head(imp_bio_df, 10), row.names = FALSE)
  cat("\n")
}

# ============================================================================
# STEP 2: Classify validation patients into benefit / indeterminate / harm
#
# Thresholds are derived from the validation CATE scores themselves.
# Convention: CATE < 0 means treatment reduces mortality (benefit).
#   - Benefit:       score < tertile_33 (bottom third, most benefit)
#   - Indeterminate: tertile_33 <= score <= tertile_67
#   - Harm:          score > tertile_67 (top third, least benefit / harm)
#
# We do this manually rather than calling confirm_analysis(out="table")
# to avoid the browser() call on line 74 of confirm_analysis.R.
# ============================================================================

cat("============================================================\n")
cat("STEP 2: Subgroup classification and treatment effect testing\n")
cat("============================================================\n\n")

classify_and_test <- function(preds_df, label) {
  cat("--- ", label, " ---\n")

  score <- preds_df$score
  y_val <- preds_df$y
  w_val <- preds_df$w

  # Tertile thresholds
  q33 <- quantile(score, 1/3)
  q67 <- quantile(score, 2/3)

  cat("  Score distribution: median =", round(median(score), 4),
      ", Q33 =", round(q33, 4), ", Q67 =", round(q67, 4), "\n")

  # Classify
  sub_ben   <- as.integer(score < q33)
  sub_harm  <- as.integer(score > q67)
  sub_indet <- as.integer(score >= q33 & score <= q67)

  cat("  Benefit group:       n =", sum(sub_ben), "\n")
  cat("  Indeterminate group: n =", sum(sub_indet), "\n")
  cat("  Harm group:          n =", sum(sub_harm), "\n\n")

  # Treatment effect within each subgroup (prop.test)
  test_subgroup <- function(sub_indicator, alt) {
    n_w1 <- sum(sub_indicator * w_val)
    n_w0 <- sum(sub_indicator * (1 - w_val))
    d_w1 <- sum(sub_indicator * w_val * y_val)
    d_w0 <- sum(sub_indicator * (1 - w_val) * y_val)

    if (n_w1 == 0 || n_w0 == 0) {
      return(c(N = sum(sub_indicator), n_w1 = n_w1, n_w0 = n_w0,
               deaths_w1 = d_w1, rate_w1 = NA,
               deaths_w0 = d_w0, rate_w0 = NA,
               rate_diff = NA, pval = NA))
    }

    rate_w1 <- d_w1 / n_w1
    rate_w0 <- d_w0 / n_w0
    rate_diff <- rate_w1 - rate_w0

    pv <- tryCatch(
      prop.test(x = c(d_w1, d_w0), n = c(n_w1, n_w0),
                alternative = alt)$p.value,
      error = function(e) NA
    )

    c(N = sum(sub_indicator), n_w1 = n_w1, n_w0 = n_w0,
      deaths_w1 = d_w1, rate_w1 = round(rate_w1, 4),
      deaths_w0 = d_w0, rate_w0 = round(rate_w0, 4),
      rate_diff = round(rate_diff, 4), pval = round(pv, 4))
  }

  resmat <- rbind(
    benefit       = test_subgroup(sub_ben,   "less"),       # expect w1 < w0
    indeterminate = test_subgroup(sub_indet, "two.sided"),
    harm          = test_subgroup(sub_harm,  "greater")     # expect w1 > w0
  )

  cat("  Subgroup treatment effects (inhosp90):\n")
  print(as.data.frame(resmat))
  cat("\n")

  as.data.frame(resmat)
}

resmat_clinical <- classify_and_test(preds_clin_df, "Clinical only (cf.CATE)")
resmat_bio      <- classify_and_test(preds_bio_df,  "Clinical + biomarker (cf.CATE)")

# ---- Save subgroup tables ----
resmat_clinical$model <- "clinical_only"
resmat_bio$model      <- "clinical_plus_biomarker"
resmat_clinical$subgroup <- rownames(resmat_clinical)
resmat_bio$subgroup      <- rownames(resmat_bio)

combined <- rbind(resmat_clinical, resmat_bio)
write.csv(combined,
          file.path(out_dir, "confirmatory_subgroups.csv"),
          row.names = FALSE)

# ============================================================================
# STEP 3: Summary
# ============================================================================

cat("============================================================\n")
cat("SUMMARY\n")
cat("============================================================\n\n")

cat("If the benefit subgroup shows a significant treatment effect\n")
cat("(rate_w1 < rate_w0, p < 0.05 one-sided), the derivation\n")
cat("finding replicates. If not, the HTE signal may be noise.\n\n")

cat("Key comparison: does the clinical+biomarker model identify\n")
cat("a benefit subgroup with a LARGER treatment effect (rate_diff)\n")
cat("or a MORE significant p-value than the clinical-only model?\n\n")

cat("All outputs written to:\n  ", out_dir, "\n")

sink()
close(report_con)
cat("\nFull report saved to:\n", report_path, "\n")