# ============================================================================
# 08_gates_inference.R
#
# This script will run the analysis for GATES (Group Average Treatment Effects) and BLP (Best Linear
# Predictor). However, the calibration test on the validation-set predictions is made by
# 10_confirmatory_analysis.R.
#
# The BLP p-value is a test of whether CATE scores capture real
# treatment effect heterogeneity (primary validation metric)
#
# Inputs:
#   outputs/10_confirmatory_analysis/confirmatory_preds_clinical.csv
#   outputs/10_confirmatory_analysis/confirmatory_preds_with_biomarker.csv
#
# Outputs:
#   outputs/08_gates_inference/gates_results_clinical.csv
#   outputs/08_gates_inference/gates_results_with_biomarker.csv
#   outputs/08_gates_inference/gates_summary.csv
#   outputs/08_gates_inference/gates_report.txt
#
# Sources:
#   cate-repo/gates_inference.R  (gates.fun)
# ============================================================================

source(file.path(getwd(), "config", "config.R"))
out_dir <- make_output_subdir("08_gates_inference")

# ---- Source shared CATE repo functions ----
source(GATES_INFERENCE_R)

# ---- Load validation predictions ----
preds_clin <- read.csv(file.path(OUTPUTS_DIR,
  "10_confirmatory_analysis", "confirmatory_preds_clinical.csv"))
preds_bio  <- read.csv(file.path(OUTPUTS_DIR,
  "10_confirmatory_analysis", "confirmatory_preds_with_biomarker.csv"))

cat("Loaded clinical predictions:", nrow(preds_clin), "rows\n")
cat("Loaded biomarker predictions:", nrow(preds_bio), "rows\n")

# ---- Verify required columns ----
required_cols <- c("y", "w", "score", "b", "s")
stopifnot(all(required_cols %in% names(preds_clin)))
stopifnot(all(required_cols %in% names(preds_bio)))

# ---- Settings ----
pi    <- 0.5    # known randomization probability (1:1 trial)
num.q <- 5      # number of quantile groups (quintiles)

# ---- Helper: run gates.fun on a single prediction data frame ----
# gates.fun expects master = list of data frames (one per CV fold).
# For validation (single held-out set, no repeated splits), we pass a
# list of length 1, with loo=TRUE to use standard z = 1.96 CIs and
# to NOT double the p-values (the doubling is for internal CV medians).
run_gates <- function(preds_df, label) {
  cat("\n============================================================\n")
  cat("Running GATES on:", label, "\n")
  cat("============================================================\n")

  master <- list(preds_df)

  t0 <- Sys.time()
  gates_out <- gates.fun(master    = master,
                         pi        = pi,
                         num.q     = num.q,
                         loo       = TRUE,
                         irev      = NULL,
                         lintest   = TRUE,
                         unknown.pi = FALSE,
                         alpha     = 0.05,
                         vcov.type = "HC3")
  elapsed <- round(difftime(Sys.time(), t0, units = "secs"), 1)
  cat("  Runtime:", elapsed, "secs\n")

  # ---- Print key results ----
  cat("\n  Per-quintile group ATEs (median theta):\n")
  for (q in 1:num.q) {
    cat(sprintf("    Q%d: %+.4f  [%.4f, %.4f]\n",
        q, gates_out$`Theta Median`[q],
        gates_out$`Theta Lower Median`[q],
        gates_out$`Theta Upper Median`[q]))
  }
  cat(sprintf("\n  Monotonicity p-value:  %.4f\n", gates_out$`Pval Monoton`))
  cat(sprintf("  Extreme-bin p-value:   %.4f\n", gates_out$`Pval Extreme`))
  cat(sprintf("  BLP slope p-value:     %.4f  <--- primary validation test\n",
              gates_out$`Pval BLP test`))
  cat(sprintf("  BLP intercept coef:    %.4f\n", gates_out$`Pval BLP coef intercept`))
  cat(sprintf("  BLP slope coef:        %.4f\n", gates_out$`Pval BLP coef slope`))

  gates_out
}

# ---- Run both passes ----
sink_file <- file.path(out_dir, "gates_report.txt")
sink(sink_file, split = TRUE)

cat("GATES Inference on Validation Set\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("Validation n:", nrow(preds_clin), "\n")
cat("Quantile groups:", num.q, "\n")
cat("Randomization probability (pi):", pi, "\n")

gates_clin <- run_gates(preds_clin, "Clinical only (35 covariates)")
gates_bio  <- run_gates(preds_bio,  "Clinical + biomarker (43 covariates)")

# ---- Build per-quintile results tables ----
build_quintile_table <- function(gates_out, label) {
  data.frame(
    model       = label,
    quintile    = 1:num.q,
    theta       = gates_out$`Theta Median`,
    theta_lower = gates_out$`Theta Lower Median`,
    theta_upper = gates_out$`Theta Upper Median`,
    mean_score  = gates_out$`Mean prediction in quantiles`
  )
}

qtab_clin <- build_quintile_table(gates_clin, "clinical_only")
qtab_bio  <- build_quintile_table(gates_bio,  "clinical_plus_biomarker")

write.csv(qtab_clin, file.path(out_dir, "gates_results_clinical.csv"),
          row.names = FALSE)
write.csv(qtab_bio, file.path(out_dir, "gates_results_with_biomarker.csv"),
          row.names = FALSE)

# ---- Build summary table ----
summary_df <- data.frame(
  model               = c("clinical_only", "clinical_plus_biomarker"),
  pval_monotonicity    = c(gates_clin$`Pval Monoton`,  gates_bio$`Pval Monoton`),
  pval_extreme_bin     = c(gates_clin$`Pval Extreme`,  gates_bio$`Pval Extreme`),
  pval_blp_slope       = c(gates_clin$`Pval BLP test`, gates_bio$`Pval BLP test`),
  blp_intercept_coef   = c(gates_clin$`Pval BLP coef intercept`,
                           gates_bio$`Pval BLP coef intercept`),
  blp_slope_coef       = c(gates_clin$`Pval BLP coef slope`,
                           gates_bio$`Pval BLP coef slope`)
)

write.csv(summary_df, file.path(out_dir, "gates_summary.csv"),
          row.names = FALSE)

cat("\n\n============================================================\n")
cat("SUMMARY\n")
cat("============================================================\n\n")
cat("                         Clinical only    +Biomarker\n")
cat(sprintf("  Monotonicity p:        %.4f           %.4f\n",
    summary_df$pval_monotonicity[1], summary_df$pval_monotonicity[2]))
cat(sprintf("  Extreme-bin p:         %.4f           %.4f\n",
    summary_df$pval_extreme_bin[1], summary_df$pval_extreme_bin[2]))
cat(sprintf("  BLP slope p:           %.4f           %.4f\n",
    summary_df$pval_blp_slope[1], summary_df$pval_blp_slope[2]))
cat(sprintf("  BLP intercept coef:    %.4f           %.4f\n",
    summary_df$blp_intercept_coef[1], summary_df$blp_intercept_coef[2]))
cat(sprintf("  BLP slope coef:        %.4f           %.4f\n",
    summary_df$blp_slope_coef[1], summary_df$blp_slope_coef[2]))

cat("\nIf BLP slope p < 0.05, the CATE scores capture real HTE.\n")
cat("If monotonicity p < 0.05, treatment effects increase across quintiles.\n")
cat("If extreme-bin p < 0.05, top and bottom quintiles differ significantly.\n")

cat("\nAll outputs written to:\n  ", out_dir, "\n")
sink()

cat("\nGATES inference complete.\n")
cat("Report saved to:", sink_file, "\n")