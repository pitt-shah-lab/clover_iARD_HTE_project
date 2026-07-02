# ============================================================================
# build_gates_forest.R
#
# Renders the GATES forest plot from gates.fun()'s output (gates_inference.R):
# per-quantile group ATE with CIs, plus the monotonicity and extreme-bin
# test p-values in the subtitle.
# ============================================================================

source(file.path(dirname(getwd()), "config", "config.R"))
out_dir <- make_output_subdir("09_build_cate_outputs")

library(ggplot2)

plot_gates_forest <- function(gates_out, method_name,
                              path_png = file.path(out_dir, paste0("gates_forest_", method_name, ".png"))) {
  num.q <- length(gates_out$`Theta Median`)
  df <- data.frame(
    bin = factor(1:num.q),
    est = gates_out$`Theta Median`,
    lo  = gates_out$`Theta Lower Median`,
    hi  = gates_out$`Theta Upper Median`
  )
  subtitle <- sprintf("Monotonicity p = %.3f | Extreme-bin p = %.3f",
                      gates_out$`Pval Monoton`, gates_out$`Pval Extreme`)
  p <- ggplot(df, aes(bin, est)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15) +
    geom_point(size = 2.5) +
    labs(title = paste("GATES --", method_name), subtitle = subtitle,
         x = "Score quantile (1 = predicted most benefit)", y = "Group ATE") +
    theme_minimal()
  ggsave(path_png, p, width = 6, height = 4.5, dpi = 300)
  cat("Wrote:", path_png, "\n")
  invisible(p)
}

# ---------------------------------------------------------------------------
# Example usage:
# gates_out <- gates.fun(master_list, pi = 0.5, num.q = 5)
# plot_gates_forest(gates_out, "cf.CATE")
# ---------------------------------------------------------------------------