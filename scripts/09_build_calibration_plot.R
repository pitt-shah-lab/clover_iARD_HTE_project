# ============================================================================
# build_calibration_plot.R
#
# Renders a binned calibration check: mean predicted CATE within each score
# quantile (x) vs. the realized group ATE within that quantile (y), from
# gates.fun()'s output (gates_inference.R). This is a coarser stand-in for
# the full best-linear-predictor (BLP) line -- gates.fun() doesn't currently
# return the fitted BLP line itself, only test statistics.
#
# NOTE on an upstream bug: in gates_inference.R, the returned
# "Pval BLP coef intercept"/"Pval BLP coef slope" reference `blp.summary`
# *after* the per-repetition loop closes, so they reflect only the LAST
# repetition's fit rather than a median across repetitions (unlike
# "Pval BLP test", which correctly aggregates via pval_lintest[i] across the
# whole loop). If you want the true BLP slope/intercept rendered here rather
# than this binned proxy, fix that in gates_inference.R first by collecting
# blp_intercept[i] <- coef(fit.blp)[1] and blp_slope[i] <- coef(fit.blp)[2]
# inside the loop and taking median() after it.
# ============================================================================

source(file.path(dirname(getwd()), "config", "config.R"))
out_dir <- make_output_subdir("09_build_cate_outputs")

library(ggplot2)

plot_calibration <- function(gates_out, method_name,
                             path_png = file.path(out_dir, paste0("calibration_", method_name, ".png"))) {
  df <- data.frame(
    mean_score = gates_out$`Mean prediction in quantiles`,
    realized_ATE = gates_out$`Theta Median`
  )
  p <- ggplot(df, aes(mean_score, realized_ATE)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(size = 2.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
    labs(title = paste("Calibration --", method_name),
         x = "Mean predicted CATE within bin", y = "Realized group ATE within bin") +
    theme_minimal()
  ggsave(path_png, p, width = 5.5, height = 4.5, dpi = 300)
  cat("Wrote:", path_png, "\n")
  invisible(p)
}

# ---------------------------------------------------------------------------
# Example usage:
# gates_out <- gates.fun(master_list, pi = 0.5, num.q = 5)
# plot_calibration(gates_out, "cf.CATE")
# ---------------------------------------------------------------------------