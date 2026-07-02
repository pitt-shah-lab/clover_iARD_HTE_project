# ============================================================================
# build_confirm_table.R
#
# Builds the confirmatory subgroup table (benefit / indeterminate / harm)
# from confirm_analysis()'s resmat output (confirm_analysis.R), called with
# out != "preds". Analogous to the treatment-effect-by-subgroup tables in
# the CLOVERS papers, but for benefit/harm classes derived from your CATE
# model's pre-specified threshold rather than pre-specified clinical
# subgroups.
# ============================================================================

source(file.path(dirname(getwd()), "config", "config.R"))
out_dir <- make_output_subdir("09_build_cate_outputs")

build_confirm_table <- function(resmat,
                                path_csv = file.path(out_dir, "confirmatory_subgroup_table.csv")) {
  tab <- as.data.frame(resmat)
  tab <- cbind(subgroup = rownames(resmat), tab)
  write.csv(tab, path_csv, row.names = FALSE)
  cat("Wrote:", path_csv, "\n")
  tab
}

# ---------------------------------------------------------------------------
# Example usage:
# confirm_out <- confirm_analysis(train_list, test_list,
#                                 cval_ben = ..., cval_harm = ...,
#                                 xvars = xvars, out = "table")
# build_confirm_table(confirm_out)
# ---------------------------------------------------------------------------