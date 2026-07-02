# ============================================================================
# 02b_build_variable_source_table.R
#
# PURPOSE: produces the two-column lookup table -- variable name on the
# left, source CSV on the right -- documenting where every column in the
# flat file (built by 02_build_flat_file.R) actually came from.
#
# This is a standalone, read-only script: it does not read or depend on
# 02_build_flat_file.R's output. It just documents the same column list by
# hand, so it can be run independently, and cross-checks itself against
# what 02_build_flat_file.R actually produced (if that's already been run).
#
# Paths come from config/config.R -- edit PROJECT_ROOT there, not here.
#
# Output (written to outputs/02b_build_variable_source_table/):
#   clovers_variable_source_table.csv
# ============================================================================

source(file.path(dirname(getwd()), "config", "config.R"))
out_dir <- make_output_subdir("02b_build_variable_source_table")

outcome.vars <- c("inhosp90","sofa_diff")
treat.var    <- "w"

norm.vars <- c("age","temp","rr","hr","map","sbp",
               "sofa","albumin","ln_bili",
               "cr","bun","ln_g","ln_lac","sqrt_plt",
               "ln_wbc","hgb","na","bicarb","prefluid",
               "ln_bmi","gcs","charlson","o2sat","s2f")

bin.vars <- c("site_lung","site_abdom","site_urine","mv",
              "vaso","ards","dial","chf","copd",
              "liver","kidney")

bio_raw_cols <- c("il1_pg_ml","ang1_pg_ml","ang2_pg_ml","tnfr_calibrated_rplex",
                  "il6_pg_ml","strem1_pg_ml","kim1_pg_ml","srage_pg_ml")
bio.vars <- paste0("ln_", bio_raw_cols)

flat_cols <- c("id", outcome.vars, treat.var, norm.vars, bin.vars, bio.vars)

# ============================================================================
# Variable -> source table
# ============================================================================

var_source <- data.frame(variable = character(0), source = character(0),
                          stringsAsFactors = FALSE)

add_rows <- function(vars, source_label) {
  data.frame(variable = vars, source = source_label, stringsAsFactors = FALSE)
}

var_source <- rbind(
  var_source,
  add_rows("id", "yw.csv / xvars_egdt.csv / xvars_other.csv / Share.4.29.26.csv (merge key, present in all)"),
  add_rows(outcome.vars, "yw.csv"),
  add_rows(treat.var, "yw.csv"),
  add_rows(c("age","temp","rr","hr","map","sbp","sofa","albumin",
             "cr","bun","hgb","gcs","charlson","o2sat",
             "site_lung","site_abdom","site_urine","mv","vaso"), "xvars_egdt.csv"),
  add_rows("ln_bili",   "xvars_egdt.csv (raw column: bili, log-transformed)"),
  add_rows("ln_g",      "xvars_egdt.csv (raw column: g, log-transformed)"),
  add_rows("ln_lac",    "xvars_egdt.csv (raw column: lac, log-transformed)"),
  add_rows("sqrt_plt",  "xvars_egdt.csv (raw column: plt, sqrt-transformed)"),
  add_rows("ln_wbc",    "xvars_egdt.csv (raw column: wbc, log-transformed)"),
  add_rows(c("na","bicarb","prefluid","s2f","ards","dial","chf","copd",
             "liver","kidney"), "xvars_other.csv"),
  add_rows("ln_bmi",    "xvars_other.csv (raw column: bmi, log-transformed)"),
  add_rows(c("ln_il1_pg_ml","ln_ang1_pg_ml","ln_ang2_pg_ml",
             "ln_tnfr_calibrated_rplex","ln_il6_pg_ml","ln_strem1_pg_ml",
             "ln_kim1_pg_ml","ln_srage_pg_ml"),
           "Share.4.29.26.csv, visit==V1 only (raw pg/mL columns, log-transformed)")
)

source_table_path <- file.path(out_dir, "clovers_variable_source_table.csv")
write.csv(var_source, source_table_path, row.names = FALSE)
cat("Wrote", source_table_path, "--", nrow(var_source), "rows\n")

# ---- self-check against the expected column list above ----
missing_from_table <- setdiff(flat_cols, var_source$variable)
extra_in_table      <- setdiff(var_source$variable, flat_cols)
if (length(missing_from_table) > 0) {
  cat("\nWARNING: these expected columns have NO entry in the source table:\n")
  print(missing_from_table)
}
if (length(extra_in_table) > 0) {
  cat("\nWARNING: these source-table entries are NOT in the expected column list:\n")
  print(extra_in_table)
}
if (length(missing_from_table) == 0 && length(extra_in_table) == 0) {
  cat("Source table matches the expected flat-file column list exactly.\n")
}

# ---- cross-check against the ACTUAL flat file, if it's already been built ----
flat_file_path <- file.path(OUTPUTS_DIR, "02_build_flat_file", "clovers_flat_file.csv")
if (file.exists(flat_file_path)) {
  actual_cols <- names(read.csv(flat_file_path, nrows = 1))
  mismatch_a <- setdiff(actual_cols, var_source$variable)
  mismatch_b <- setdiff(var_source$variable, actual_cols)
  cat("\n-- Cross-check against the actual flat file at:\n  ", flat_file_path, "\n")
  if (length(mismatch_a) == 0 && length(mismatch_b) == 0) {
    cat("Matches exactly -- every column in the real flat file has a source entry.\n")
  } else {
    if (length(mismatch_a) > 0) {
      cat("WARNING: flat file has columns not documented here:\n"); print(mismatch_a)
    }
    if (length(mismatch_b) > 0) {
      cat("WARNING: this table documents columns not in the flat file:\n"); print(mismatch_b)
    }
  }
} else {
  cat("\n(02_build_flat_file.R hasn't been run yet -- skipping cross-check against\n")
  cat("the actual flat file. The self-check above still applies.)\n")
}