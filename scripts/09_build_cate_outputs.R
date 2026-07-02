# ============================================================================
# 09_build_cate_outputs.R
#
# Sources all the individual build_*.R scripts below so their functions are
# available in one call. Each script is standalone and can also be sourced
# on its own -- see the header comment in each for what it needs and an
# example call.
#
#   build_toc_curve.R          -> plot_toc_curve()
#   build_rate_table.R         -> build_rate_table()          (eTable 6)
#   build_gates_forest.R       -> plot_gates_forest()
#   build_gates_table.R        -> build_gates_table()
#   build_calibration_plot.R   -> plot_calibration()
#   build_confirm_table.R      -> build_confirm_table()
#   build_importance_table.R   -> build_importance_table()    (eTable 7)
#   build_importance_plot.R    -> plot_importance_bar()
#   build_interaction_table.R  -> build_interaction_table()
#
# Assumes all build_*.R files live in the same directory as this script.
# ============================================================================

this_file   <- sys.frame(1)$ofile
script_dir  <- if (!is.null(this_file)) dirname(this_file) else getwd()

for (f in c("build_toc_curve.R",
            "build_rate_table.R",
            "build_gates_forest.R",
            "build_gates_table.R",
            "build_calibration_plot.R",
            "build_confirm_table.R",
            "build_importance_table.R",
            "build_importance_plot.R",
            "build_interaction_table.R")) {
  source(file.path(script_dir, f))
}

cat("All build_*.R functions loaded.\n")