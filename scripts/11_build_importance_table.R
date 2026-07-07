# ============================================================================
# build_importance_table.R
#
# Builds eTable 7 -- the biomarker/variable-importance contribution table --
# from the per-fold `importance` vectors that every function in
# scoring_methods_v4.R returns as its 2nd list element. Collect these across
# cv_k folds for one method (e.g. from res.allmods[["cf.CATE"]], each
# fold's [[2]]) before calling this.
# ============================================================================

this_script <- sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])
source(file.path(dirname(dirname(this_script)), "config", "config.R"))
out_dir <- make_output_subdir("09_build_cate_outputs")

build_importance_table <- function(importance_list, xvars, method_name,
                                    path_csv = file.path(out_dir, "eTable_7_biomarker_contribution.csv")) {
  # importance_list: list of numeric vectors (one per CV fold), each of
  # length(xvars), in the same order as xvars was passed to the scoring fn.
  imp_mat <- do.call(cbind, importance_list)
  mean_imp <- rowMeans(imp_mat, na.rm = TRUE)
  sd_imp   <- apply(imp_mat, 1, sd, na.rm = TRUE)
  tab <- data.frame(method = method_name, variable = xvars,
                    mean_importance = mean_imp, sd_importance = sd_imp)
  tab <- tab[order(-tab$mean_importance), ]
  write.csv(tab, path_csv, row.names = FALSE, append = file.exists(path_csv))
  cat("Wrote/updated:", path_csv, "\n")
  tab
}

# ---------------------------------------------------------------------------
# Example usage:
# imp_list <- lapply(res.allmods[["cf.CATE"]], function(f) f[[2]])
# imp_tab  <- build_importance_table(imp_list, xvars, "cf.CATE")
# ---------------------------------------------------------------------------