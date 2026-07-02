# ============================================================================
# 09_build_cate_outputs_FIXED.R
#
# Corrected version of 09_build_cate_outputs.R, updated to source() the
# ACTUAL current filenames after the recent renames:
#
#   build_gates_forest.R      -> 09_build_forest.R
#   build_confirm_table.R     -> 09_build_confirm_grp_table.R
#   build_interaction_table.R -> 09_build_interaction_plot.R
#
# UNRESOLVED before this will fully work:
#   - build_toc_curve.R currently contains build_rate_table.R's code
#     (confirmed by diff -- identical body, different header comment only).
#     plot_toc_curve() does not exist under any filename right now. Restore
#     the real build_toc_curve.R content (see the version from earlier in
#     this conversation, or rewrite it) before this loader will expose
#     plot_toc_curve().
#   - build_gates_table.R was not part of the last upload, so I can't
#     confirm whether it's still named that or was renamed too. Check
#     before relying on build_gates_table().
# ============================================================================

this_file   <- sys.frame(1)$ofile
script_dir  <- if (!is.null(this_file)) dirname(this_file) else getwd()

for (f in c("build_toc_curve.R",            # still broken -- see note above
            "build_rate_table.R",
            "09_build_forest.R",            # was build_gates_forest.R
            "build_gates_table.R",          # unconfirmed current name
            "09_build_calibration_plot.R",  # was build_calibration_plot.R
            "09_build_confirm_grp_table.R", # was build_confirm_table.R
            "build_importance_table.R",
            "09_build_importance_plot.R",   # was build_importance_plot.R
            "09_build_interaction_plot.R")) { # was build_interaction_table.R
  path <- file.path(script_dir, f)
  if (file.exists(path)) {
    source(path)
  } else {
    warning("Could not find ", f, " in ", script_dir,
            " -- functions from that file will not be available.")
  }
}

cat("Done loading build_*/09_*.R functions (see warnings above for any misses).\n")