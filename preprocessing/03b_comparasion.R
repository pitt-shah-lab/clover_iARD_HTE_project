# ============================================================================
# 03b_compare_der_val.R
#
# Objective (Victor's stage 2b): table comparing every covariate,
# outcome, and treatment variable between the derivation and validation
# halves, to confirm the split produced two equivalent groups.
#
# For continuous variables: mean, SD, and a Welch two-sample t-test p-value.
# For binary variables: proportion in each half, and a chi-squared test
# p-value (Fisher's exact if any expected cell count is small).
#
# Large p-values indicate the two halves are statistically equivalent, as expected
# from a proper random split. A few values of < ,.05 is okay. 
#
# Paths come from config/config.R
#
# Input:
#   outputs/03_split_derivation_validation/flat_file_der.csv
#   outputs/03_split_derivation_validation/flat_file_val.csv
# Output (written to outputs/03b_compare_der_val/):
#   der_val_comparison_table.csv
#   der_val_comparison_report.txt
# ============================================================================

source(file.path("config", "config.R"))

out_dir <- make_output_subdir("03b_compare_der_val")

split_dir <- file.path(OUTPUTS_DIR, "03_split_derivation_validation")
der_path <- file.path(split_dir, "flat_file_der.csv")
val_path <- file.path(split_dir, "flat_file_val.csv")

if (!file.exists(der_path) || !file.exists(val_path)) {
  stop("Derivation/validation files not found in ", split_dir,
       " -- run 03_split_derivation_validation.R first.")
}

der <- read.csv(der_path, stringsAsFactors = FALSE)
val <- read.csv(val_path, stringsAsFactors = FALSE)

report_path <- file.path(out_dir, "der_val_comparison_report.txt")
report_con <- file(report_path, open = "wt")
sink(report_con, split = TRUE)

cat("============================================================\n")
cat("DERIVATION vs VALIDATION COMPARISON\n")
cat("============================================================\n")
cat("Derivation n:", nrow(der), " | Validation n:", nrow(val), "\n\n")

# ---- variable classification ----
# binary covariates from the flat file, used as-is (0/1)
binary_vars <- c("w","inhosp90","site_lung","site_abdom","site_urine","mv",
                  "vaso","ards","dial","chf","copd","liver","kidney")

# everything else numeric in the flat file is treated as continuous
exclude_from_continuous <- c("id", binary_vars)
all_numeric_cols <- names(der)[sapply(der, is.numeric)]
continuous_vars <- setdiff(all_numeric_cols, exclude_from_continuous)
binary_vars <- intersect(binary_vars, names(der))   # only ones actually present

cat("Continuous variables compared:", length(continuous_vars), "\n")
cat("Binary variables compared:    ", length(binary_vars), "\n\n")

results <- data.frame(
  variable = character(0), type = character(0),
  der_summary = character(0), val_summary = character(0),
  p_value = numeric(0), stringsAsFactors = FALSE
)

# ---- continuous variables: mean (SD), Welch t-test ----
for (v in continuous_vars) {
  d <- der[[v]]; va <- val[[v]]
  d_ok <- d[!is.na(d)]; va_ok <- va[!is.na(va)]

  if (length(d_ok) < 2 || length(va_ok) < 2) {
    p <- NA
  } else {
    p <- tryCatch(t.test(d_ok, va_ok)$p.value, error = function(e) NA)
  }

  results <- rbind(results, data.frame(
    variable = v, type = "continuous",
    der_summary = sprintf("%.2f (%.2f)", mean(d_ok), sd(d_ok)),
    val_summary = sprintf("%.2f (%.2f)", mean(va_ok), sd(va_ok)),
    p_value = p, stringsAsFactors = FALSE
  ))
}

# ---- binary variables: proportion, chi-squared (or Fisher's if needed) ----
for (v in binary_vars) {
  d <- der[[v]]; va <- val[[v]]
  d_ok <- d[!is.na(d)]; va_ok <- va[!is.na(va)]

  tab <- matrix(c(sum(d_ok == 1), sum(d_ok == 0),
                  sum(va_ok == 1), sum(va_ok == 0)), nrow = 2)

  p <- tryCatch({
    ch <- suppressWarnings(chisq.test(tab))
    if (any(ch$expected < 5)) {
      fisher.test(tab)$p.value
    } else {
      ch$p.value
    }
  }, error = function(e) NA)

  results <- rbind(results, data.frame(
    variable = v, type = "binary",
    der_summary = sprintf("%d/%d (%.1f%%)", sum(d_ok == 1), length(d_ok), 100 * mean(d_ok == 1)),
    val_summary = sprintf("%d/%d (%.1f%%)", sum(va_ok == 1), length(va_ok), 100 * mean(va_ok == 1)),
    p_value = p, stringsAsFactors = FALSE
  ))
}

results <- results[order(results$p_value), ]

cat("Comparison table (sorted by p-value, smallest first):\n")
cat(sprintf("%-26s%-12s%-18s%-18s%s\n", "variable", "type", "derivation", "validation", "p_value"))
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  p_str <- if (is.na(r$p_value)) "NA" else sprintf("%.3f", r$p_value)
  cat(sprintf("%-26s%-12s%-18s%-18s%s\n", r$variable, r$type, r$der_summary, r$val_summary, p_str))
}

n_sig <- sum(results$p_value < 0.05, na.rm = TRUE)
n_tested <- sum(!is.na(results$p_value))
cat("\n", n_sig, "of", n_tested, "variables show p < 0.05.\n")
cat("With", n_tested, "variables tested at alpha=0.05, roughly",
    round(n_tested * 0.05, 1), "would be expected to cross that threshold\n")
cat("by chance alone, even with two truly equivalent groups. A small number\n")
cat("of borderline p-values is expected and not evidence of a bad split.\n")
cat("Investigate further only if a p-value is very small or the variable\n")
cat("involved is a key confounder/effect-modifier you plan to rely on.\n")

# ---- Average treatment effects per half (Victor request) ----
cat("\n============================================================\n")
cat("AVERAGE TREATMENT EFFECTS BY HALF\n")
cat("============================================================\n")

ate_per_half <- function(df, label) {
  y1 <- df$inhosp90[df$w == 1]
  y0 <- df$inhosp90[df$w == 0]
  n1 <- length(y1); n0 <- length(y0)
  r1 <- mean(y1);   r0 <- mean(y0)
  ate <- r1 - r0
  pt <- prop.test(c(sum(y1), sum(y0)), c(n1, n0))
  cat(sprintf("  %s (n=%d):\n", label, nrow(df)))
  cat(sprintf("    Treated (w=1):   %d/%d = %.1f%%\n", sum(y1), n1, 100*r1))
  cat(sprintf("    Control (w=0):   %d/%d = %.1f%%\n", sum(y0), n0, 100*r0))
  cat(sprintf("    ATE (risk diff): %+.1f%%  (p = %.3f)\n\n", 100*ate, pt$p.value))
  data.frame(half = label, n = nrow(df),
             n_w1 = n1, n_w0 = n0,
             rate_w1 = round(r1, 4), rate_w0 = round(r0, 4),
             ATE = round(ate, 4), p_value = round(pt$p.value, 4))
}

ate_der <- ate_per_half(der, "Derivation")
ate_val <- ate_per_half(val, "Validation")
ate_table <- rbind(ate_der, ate_val)

write.csv(ate_table, file.path(out_dir, "ate_by_half.csv"), row.names = FALSE)
cat("Wrote: ate_by_half.csv\n")

write.csv(results, file.path(out_dir, "der_val_comparison_table.csv"), row.names = FALSE)

sink()
close(report_con)
cat("\nWrote comparison table and full report to:\n  ", out_dir, "\n")