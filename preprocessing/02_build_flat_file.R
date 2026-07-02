# ============================================================================
# 02_build_flat_file.R
#
# PURPOSE (Victor's stage 1): build the analysis-ready flat file with every
# covariate, the outcome(s), treatment, and id -- NO imputation, NO modeling
# of any kind, including risk modeling. Leave missing data as NA.
#
# Paths come from config/config.R -- edit PROJECT_ROOT there, not here.
#
# Output (written to outputs/02_build_flat_file/):
#   clovers_flat_file.csv
#
# See 02b_build_variable_source_table.R for the companion variable -> source
# lookup table (kept as a separate script on purpose).
# ============================================================================

source(file.path(dirname(getwd()), "config", "config.R"))
out_dir <- make_output_subdir("02_build_flat_file")

# ---- load source files (paths from config.R) ----
yw    <- read.csv(YW_CSV)
egdt  <- read.csv(EGDT_CSV)
other <- read.csv(OTHER_CSV)
bio   <- read.csv(BIOMARKER_CSV)   # Share.4.29.26.csv

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
                                          # cohort with V1 biomarkers

cat("Flat file cohort:", nrow(dat), "patients\n")

# ============================================================================
# Transforms only -- same convention as the original clovers_forHTE_flat.R.
# Unit/scale conversions, NOT imputation or modeling -- belongs in stage 1.
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

flat_path <- file.path(out_dir, "clovers_flat_file.csv")
write.csv(flat, flat_path, row.names = FALSE)

cat("Wrote", flat_path, "--", nrow(flat), "rows,", ncol(flat), "columns\n")
cat("(no imputation or modeling applied -- missing values left as NA)\n")
cat("\nNext: run 02b_build_variable_source_table.R for the variable/source\n")
cat("lookup table, then 03_check_missingness_and_qc.R to QC this file.\n")