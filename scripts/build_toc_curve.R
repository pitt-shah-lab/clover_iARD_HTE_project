# ============================================================================
# build_rate_table.R
#
# Builds eTable 6 -- the CATE method-by-method comparison table -- from the
# $ratelist element of post.cv.fun()'s output (post_process_v6.R):
# AUQINI, AUTOC, Kendall Tau, Adj AUQINI (and SEs if you requested
# se.method != "none"), for every method you ran through the CV harness.
# ============================================================================

source(file.path(dirname(getwd()), "config", "config.R"))
out_dir <- make_output_subdir("09_build_cate_outputs")

build_rate_table <- function(post_out_list, method_names,
                             path_csv = file.path(out_dir, "eTable_6_cate_method_comparison.csv")) {
  # post_out_list: named list of post.cv.fun() outputs, one per CATE method
  # in scoring_methods_v4.R you're comparing.
  rows <- lapply(seq_along(post_out_list), function(i) {
    rl <- post_out_list[[i]]$ratelist
    if (is.null(rl)) return(NULL)
    do.call(rbind, lapply(names(rl), function(dir_name) {
      cbind(method = method_names[i], direction = dir_name, rl[[dir_name]])
    }))
  })
  tab <- do.call(rbind, rows)
  tab <- tab[order(-tab$AUTOC), ]  # best AUTOC first
  write.csv(tab, path_csv, row.names = FALSE)
  cat("Wrote:", path_csv, "\n")
  tab
}

# ---------------------------------------------------------------------------
# Example usage:
# post_out_cf   <- post.cv.fun(cv.obj_cf,   stats = "stats.bin", score = "cate")
# post_out_bcf  <- post.cv.fun(cv.obj_bcf,  stats = "stats.bin", score = "cate")
# build_rate_table(list(post_out_cf, post_out_bcf), c("cf.CATE", "bcf.CATE"))
# ---------------------------------------------------------------------------