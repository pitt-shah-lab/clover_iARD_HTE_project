# ============================================================================
# 02b_build_flat_file_and_variable_source_table.R
#
# PURPOSE (Victor's stage 1): build the analysis-ready flat file with every
# covariate, the outcome(s), treatment, and id -- NO imputation, NO modeling
# of any kind, including risk modeling. Leave missing data as NA.
#
# Also produces a second output: a two-column lookup table, variable name on
# the left and the source CSV it came from on the right, so anyone (Victor,
# Jaspreet, a reviewer) can see at a glance where every column originated.
#
# Outputs:
#   ../outputs/clovers_flat_file.csv             -- the analysis-ready data
#   ../outputs/clovers_variable_source_table.csv -- variable -> source table
# ============================================================================

data_dir <- "../data/Data"   # <- relative to preprocessing/, points at the unzipped Data.zip

# ---- load source files ----
yw    <- read.csv(file.path(data_dir, "Curated datasets/yw.csv"))
egdt  <- read.csv(file.path(data_dir, "Curated datasets/xvars_egdt.csv"))
other <- read.csv(file.path(data_dir, "Curated datasets/xvars_other.csv"))
bio   <- read.csv(file.path(data_dir, "Biomarker data/Share.4.29.26.csv"))

strip_index_col <- function(df) {
  if (names(df)[1] == "" || names(df)[1] == "X") df <- df[, -1]
  df
}
yw    <- strip_index_col(yw)
egdt  <- strip_index_col(egdt)
other <- strip_index_col(other)

bio_v1 <- bio[bio$visit == "V1", ]   # baseline draw only

bio_raw_cols <- c("il1_pg_ml","ang1_pg_ml","ang2_pg_ml","tnfr_calibrated_rplex",
                  "il6_pg_ml","strem1_pg_ml","kim1_pg_ml","srage_pg_ml")
# NOTE: tnfr_uplex (the alternate TNFR-1 platform) is excluded here in favor
# of tnfr_calibrated_rplex -- confirm this choice with Victor.
bio_v1 <- bio_v1[, c("id", bio_raw_cols)]

# ---- merge everything on id ----
dat <- merge(yw,    egdt,   by = "id")
dat <- merge(dat,   other,  by = "id")
dat <- merge(dat,   bio_v1, by = "id")   # inner join -> restricts to the
                                          # 1,340 patients with V1 biomarkers

cat("Flat file cohort:", nrow(dat), "patients\n")

# ============================================================================
# Transforms only -- same convention as the original clovers_forHTE_flat.R.
# These are unit/scale conversions, NOT imputation or modeling, so they
# belong in stage 1 per Victor's instructions.
# ============================================================================

dat$bili[dat$bili == 0] <- NA
dat$ln_bili <- log(dat$bili)
dat$ln_g    <- log(dat$g)
dat$ln_lac  <- log(dat$lac)
dat$sqrt_plt <- sqrt(dat$plt)
dat$wbc[dat$wbc == 0] <- NA
dat$ln_wbc  <- log(dat$wbc)
dat$ln_bmi  <- log(dat$bmi)

for (col in bio_raw_cols) {
  ln_col <- paste0("ln_", col)
  vals <- dat[[col]]
  vals[vals <= 0] <- NA
  dat[[ln_col]] <- log(vals)
}

# ============================================================================
# Final column selection for the flat file
# ============================================================================

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

bio.vars <- paste0("ln_", bio_raw_cols)

flat_cols <- c("id", outcome.vars, treat.var, norm.vars, bin.vars, bio.vars)
flat <- dat[, flat_cols]

write.csv(flat, "../outputs/clovers_flat_file.csv", row.names = FALSE)
cat("Wrote ../outputs/clovers_flat_file.csv --", nrow(flat), "rows,", ncol(flat), "columns\n")
cat("(no imputation or modeling applied -- missing values left as NA)\n\n")

# ============================================================================
# Variable -> source table
# ============================================================================
# Two columns: variable name (as it appears in the flat file above) and the
# source CSV it was pulled from. Transformed variables (ln_*, sqrt_*) list
# the source of their RAW input column, since that's where the data
# physically comes from.

var_source <- data.frame(
  variable = character(0),
  source   = character(0),
  stringsAsFactors = FALSE
)

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

write.csv(var_source, "../outputs/clovers_variable_source_table.csv", row.names = FALSE)
cat("Wrote ../outputs/clovers_variable_source_table.csv --", nrow(var_source), "rows\n")

# sanity check: every column in the flat file should have exactly one row
# in the source table, and vice versa
missing_from_table <- setdiff(flat_cols, var_source$variable)
extra_in_table      <- setdiff(var_source$variable, flat_cols)
if (length(missing_from_table) > 0) {
  cat("\nWARNING: these flat-file columns have NO entry in the source table:\n")
  print(missing_from_table)
}
if (length(extra_in_table) > 0) {
  cat("\nWARNING: these source-table entries do NOT appear in the flat file:\n")
  print(extra_in_table)
}
if (length(missing_from_table) == 0 && length(extra_in_table) == 0) {
  cat("Source table matches flat file exactly -- every column accounted for.\n")
}