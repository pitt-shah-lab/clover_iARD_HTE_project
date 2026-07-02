# ============================================================================
# 02c_check_missingness_and_qc.R
#
# PURPOSE: Run this AFTER 02_build_flat_file.R (and, optionally, after
# 02b_build_variable_source_table.R -- order between 02b and 02c doesn't
# matter, both just read 02's output).
# Loads that script's output and reports missingness, sanity-checks the
# log-transformed biomarkers, and gives a clear pass/fail on imputation
# readiness.
#
# Paths come from config/config.R -- edit PROJECT_ROOT there, not here.
#
# Does not modify anything. Safe to re-run.
# ============================================================================

source(file.path(dirname(getwd()), "config", "config.R"))
out_dir <- make_output_subdir("02c_check_missingness_and_qc")

# ---- capture everything below to a single text file, AND keep printing
# to the console at the same time (split = TRUE) ----
report_path <- file.path(out_dir, "qc_report.txt")
report_con <- file(report_path, open = "wt")
sink(report_con, split = TRUE)

in_path <- file.path(OUTPUTS_DIR, "02_build_flat_file",
                      "clovers_flat_file.csv")
dat <- read.csv(in_path)

cat("============================================================\n")
cat("Loaded:", in_path, "--", nrow(dat), "rows,", ncol(dat), "columns\n")
cat("============================================================\n\n")

# ---- 1. missingness per column ----
cat("Missingness by column (only columns with >0 missing shown):\n")
miss_counts <- sapply(dat, function(col) sum(is.na(col)))
miss_counts <- miss_counts[miss_counts > 0]
miss_pct <- round(100 * miss_counts / nrow(dat), 1)
miss_tbl <- data.frame(column = names(miss_counts), n_missing = miss_counts, pct_missing = miss_pct)
miss_tbl <- miss_tbl[order(-miss_tbl$n_missing), ]
print(miss_tbl, row.names = FALSE)

miss_tbl_path <- file.path(out_dir, "missingness_by_column.csv")
write.csv(miss_tbl, miss_tbl_path, row.names = FALSE)
cat("\nWrote:", miss_tbl_path, "\n\n")

# ---- 2. log-transformed biomarker sanity check ----
cat("============================================================\n")
cat("Log-transformed biomarker ranges (should all be finite, no -Inf/NaN)\n")
cat("============================================================\n")

ln_bio_cols <- grep("^ln_(il1|ang1|ang2|tnfr|il6|strem1|kim1|srage)", names(dat), value = TRUE)
for (col in ln_bio_cols) {
  v <- dat[[col]]
  n_inf <- sum(is.infinite(v), na.rm = TRUE)
  n_nan <- sum(is.nan(v))
  cat(sprintf("%-26s  min=%8.2f  max=%8.2f  n_missing=%4d  n_-Inf=%2d  n_NaN=%2d\n",
              col, suppressWarnings(min(v, na.rm = TRUE)), suppressWarnings(max(v, na.rm = TRUE)),
              sum(is.na(v)), n_inf, n_nan))
}
cat("\n")

# ---- 3. treatment/outcome sanity check ----
cat("============================================================\n")
cat("Treatment and outcome distribution in this cohort\n")
cat("============================================================\n")
cat("Treatment (w):\n");      print(table(dat$w))
cat("Outcome (inhosp90):\n"); print(table(dat$inhosp90))
cat("sofa_diff missing:", sum(is.na(dat$sofa_diff)), "of", nrow(dat), "\n\n")

# ---- 4. final readiness check ----
cat("============================================================\n")
cat("READINESS CHECK\n")
cat("============================================================\n")
covariate_cols <- setdiff(names(dat), c("id","inhosp90","sofa_diff","w"))
any_missing <- sum(miss_counts[names(miss_counts) %in% covariate_cols])
cat("Total missing covariate cells across the dataset:", any_missing, "\n")

if (any_missing > 0) {
  cat("--> NOT yet ready for causal_forest()/bart(). Run imputation\n")
  cat("    (missforest_train.R, sourced via config.R) on the covariate\n")
  cat("    columns first.\n")
} else {
  cat("--> No missing covariate values. Ready for the CATE scoring functions.\n")
}

sink()
close(report_con)
cat("\nFull QC report saved to:\n", report_path, "\n")