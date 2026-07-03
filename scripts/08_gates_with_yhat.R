# ============================================================================
# 08b_gates_with_yhat.R
#
# Runs GATES / BLP on the Y.hat-corrected validation predictions from
# 10b_confirmatory_analysis_with_yhat.R.
#
# Original results (without Y.hat) remain in outputs/08_gates_inference/.
# These results go to outputs/08b_gates_with_yhat/.
# ============================================================================

this_script <- sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])
source(file.path(dirname(dirname(this_script)), "config", "config.R"))
out_dir <- make_output_subdir("08b_gates_with_yhat")

source(GATES_INFERENCE_R)

# ---- Load Y.hat predictions ----
preds_clin <- read.csv(file.path(OUTPUTS_DIR,
  "10b_confirmatory_with_yhat", "confirmatory_preds_clinical.csv"))
preds_bio  <- read.csv(file.path(OUTPUTS_DIR,
  "10b_confirmatory_with_yhat", "confirmatory_preds_with_biomarker.csv"))

cat("Loaded clinical predictions:", nrow(preds_clin), "rows\n")
cat("Loaded biomarker predictions:", nrow(preds_bio), "rows\n")

required_cols <- c("y", "w", "score", "b", "s")
stopifnot(all(required_cols %in% names(preds_clin)))
stopifnot(all(required_cols %in% names(preds_bio)))

pi    <- 0.5
num.q <- 5

run_gates <- function(preds_df, label) {
  cat("\n============================================================\n")
  cat("Running GATES on:", label, "\n")
  cat("============================================================\n")

  master <- list(preds_df)
  t0 <- Sys.time()
  gates_out <- gates.fun(master = master, pi = pi, num.q = num.q,
                         loo = TRUE, irev = NULL, lintest = TRUE,
                         unknown.pi = FALSE, alpha = 0.05, vcov.type = "HC3")
  cat("  Runtime:", round(difftime(Sys.time(), t0, units = "secs"), 1), "secs\n")

  cat("\n  Per-quintile group ATEs:\n")
  for (q in 1:num.q) {
    cat(sprintf("    Q%d: %+.4f  [%.4f, %.4f]\n",
        q, gates_out$`Theta Median`[q],
        gates_out$`Theta Lower Median`[q],
        gates_out$`Theta Upper Median`[q]))
  }
  cat(sprintf("\n  Monotonicity p:  %.4f\n", gates_out$`Pval Monoton`))
  cat(sprintf("  Extreme-bin p:   %.4f\n", gates_out$`Pval Extreme`))
  cat(sprintf("  BLP slope p:     %.4f  <--- primary validation test\n",
              gates_out$`Pval BLP test`))
  cat(sprintf("  BLP intercept:   %.4f\n", gates_out$`Pval BLP coef intercept`))
  cat(sprintf("  BLP slope coef:  %.4f\n", gates_out$`Pval BLP coef slope`))

  gates_out
}

sink_file <- file.path(out_dir, "gates_report.txt")
sink(sink_file, split = TRUE)

cat("GATES Inference — WITH Y.hat (risk-model predictions)\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("Validation n:", nrow(preds_clin), "\n\n")

gates_clin <- run_gates(preds_clin, "Clinical only + Y.hat")
gates_bio  <- run_gates(preds_bio,  "Clinical + biomarker + Y.hat")

# ---- Save results ----
summary_df <- data.frame(
  model             = c("clinical_with_yhat", "biomarker_with_yhat"),
  pval_monotonicity  = c(gates_clin$`Pval Monoton`,  gates_bio$`Pval Monoton`),
  pval_extreme_bin   = c(gates_clin$`Pval Extreme`,  gates_bio$`Pval Extreme`),
  pval_blp_slope     = c(gates_clin$`Pval BLP test`, gates_bio$`Pval BLP test`),
  blp_intercept_coef = c(gates_clin$`Pval BLP coef intercept`,
                         gates_bio$`Pval BLP coef intercept`),
  blp_slope_coef     = c(gates_clin$`Pval BLP coef slope`,
                         gates_bio$`Pval BLP coef slope`)
)

write.csv(summary_df, file.path(out_dir, "gates_summary.csv"), row.names = FALSE)

# ---- Side-by-side with original ----
cat("\n\n============================================================\n")
cat("COMPARISON: without Y.hat vs with Y.hat\n")
cat("============================================================\n\n")

old_path <- file.path(OUTPUTS_DIR, "08_gates_inference", "gates_summary.csv")
if (file.exists(old_path)) {
  old <- read.csv(old_path)
  cat(sprintf("%-30s %10s %10s\n", "", "Without Y.hat", "With Y.hat"))
  cat(sprintf("%-30s %10.4f %10.4f\n", "Clinical BLP slope p:",
      old$pval_blp_slope[1], summary_df$pval_blp_slope[1]))
  cat(sprintf("%-30s %10.4f %10.4f\n", "+Biomarker BLP slope p:",
      old$pval_blp_slope[2], summary_df$pval_blp_slope[2]))
} else {
  cat("Original GATES results not found for comparison.\n")
}

cat("\n\nAll outputs written to:\n  ", out_dir, "\n")
sink()
cat("\nDone. Report:", sink_file, "\n")