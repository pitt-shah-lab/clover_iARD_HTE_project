# ============================================================================
# build_interaction_table.R
#
# Builds a summary table from michaelq()'s output (michael_q.R): the
# 2xK interaction test used e.g. across quantile bins or trial sites.
# ============================================================================

source(file.path(dirname(getwd()), "config", "config.R"))
out_dir <- make_output_subdir("09_build_cate_outputs")

build_interaction_table <- function(mq_out, label,
                                    path_csv = file.path(out_dir, "interaction_test_summary.csv")) {
  tab <- data.frame(comparison = label,
                    difference = paste(round(mq_out$Difference, 3), collapse = "; "),
                    Q = mq_out$Q, p_value = mq_out$pval)
  write.csv(tab, path_csv, row.names = FALSE, append = file.exists(path_csv))
  cat("Wrote/updated:", path_csv, "\n")
  tab
}

# ---------------------------------------------------------------------------
# Example usage:
# mq_out <- michaelq(ps, ns, alpha = 0.05)
# build_interaction_table(mq_out, "site")
# ---------------------------------------------------------------------------