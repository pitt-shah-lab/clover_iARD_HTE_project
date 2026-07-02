# ============================================================================
# build_importance_plot.R
#
# Renders the variable-importance bar chart from the table built by
# build_importance_table.R (build_importance_table()).
# ============================================================================

source(file.path(dirname(getwd()), "config", "config.R"))
out_dir <- make_output_subdir("09_build_cate_outputs")

library(ggplot2)

plot_importance_bar <- function(importance_table, method_name, top_n = 15,
                                path_png = file.path(out_dir, paste0("importance_", method_name, ".png"))) {
  df <- head(importance_table[order(-importance_table$mean_importance), ], top_n)
  df$variable <- factor(df$variable, levels = rev(df$variable))
  p <- ggplot(df, aes(variable, mean_importance)) +
    geom_col() +
    coord_flip() +
    labs(title = paste("Variable importance --", method_name),
         x = NULL, y = "Mean importance across CV folds") +
    theme_minimal()
  ggsave(path_png, p, width = 6, height = 5, dpi = 300)
  cat("Wrote:", path_png, "\n")
  invisible(p)
}

# ---------------------------------------------------------------------------
# Example usage:
# imp_tab <- build_importance_table(imp_list, xvars, "cf.CATE")  # from build_importance_table.R
# plot_importance_bar(imp_tab, "cf.CATE")
# ---------------------------------------------------------------------------