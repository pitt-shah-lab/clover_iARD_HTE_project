# ============================================================================
# 11_table1.R
#
# Builds Table 1: baseline characteristics of the biomarker analytic cohort,
# stratified by treatment arm (w=0 restrictive, w=1 liberal).
#
# Uses the tableone package to produce a publication-ready table with
# means (SD) for continuous variables, n (%) for categorical, p-values
# from t-tests / chi-squared, and standardized mean differences (SMDs).
#
# Output (written to outputs/11_table1/):
#   table1.csv                -- machine-readable
#   table1_formatted.txt      -- print-formatted for pasting into manuscript
#   table1_by_cohort.csv      -- full cohort (1563) vs biomarker cohort (1340)
#add comment
# ============================================================================

source(file.path((getwd()), "config", "config.R"))
out_dir <- make_output_subdir("11_table1")

# ---- install tableone if needed ----
if (!requireNamespace("tableone", quietly = TRUE)) {
  install.packages("tableone", repos = "https://cloud.r-project.org")
}
library(tableone)

# ---- load the flat file (pre-imputation, so missingness is visible) ----
flat_path <- file.path(OUTPUTS_DIR, "02_build_flat_file", "clovers_flat_file.csv")
if (!file.exists(flat_path)) {
  stop("Flat file not found at ", flat_path, " -- run 02_build_flat_file.R first.")
}
ds <- read.csv(flat_path, stringsAsFactors = FALSE)
cat("Loaded flat file:", nrow(ds), "patients,", ncol(ds), "columns\n\n")

# ---- define variable lists ----
# Continuous variables (display as mean ± SD or median [IQR] if skewed)
cont_vars <- c("age", "temp", "rr", "hr", "map", "sbp",
               "sofa", "albumin", "cr", "bun", "hgb",
               "na", "bicarb", "prefluid", "gcs", "charlson",
               "o2sat", "s2f")

# Variables to display as median [IQR] (skewed)
nonnormal_vars <- c("sofa", "charlson", "prefluid", "cr", "bun")

# Log/sqrt transformed biomarkers — show on original scale
bio_orig <- c("ln_il1_pg_ml", "ln_ang1_pg_ml", "ln_ang2_pg_ml",
              "ln_tnfr_calibrated_rplex", "ln_il6_pg_ml",
              "ln_strem1_pg_ml", "ln_kim1_pg_ml", "ln_srage_pg_ml")

# Binary variables
bin_vars <- c("site_lung", "site_abdom", "site_urine",
              "mv", "vaso", "ards", "dial",
              "chf", "copd", "liver", "kidney")

# All variables for Table 1
all_vars <- c(cont_vars, bio_orig, bin_vars)

# Which are categorical (for tableone)
cat_vars <- bin_vars

# ---- treatment arm labels ----
ds$arm <- ifelse(ds$w == 0, "Restrictive", "Liberal")

# ============================================================================
# TABLE 1A: Biomarker cohort by treatment arm
# ============================================================================

cat("============================================================\n")
cat("TABLE 1A: Biomarker cohort (n=1,340) by treatment arm\n")
cat("============================================================\n\n")

tab1 <- CreateTableOne(
  vars      = all_vars,
  strata    = "arm",
  data      = ds,
  factorVars = cat_vars,
  addOverall = TRUE
)

# Print with SMD and p-values
tab1_print <- print(tab1,
                     nonnormal  = nonnormal_vars,
                     printToggle = FALSE,
                     smd = TRUE,
                     test = TRUE)

cat(capture.output(print(tab1, nonnormal = nonnormal_vars, smd = TRUE)),
    sep = "\n")

# ---- Save ----
write.csv(tab1_print,
          file.path(out_dir, "table1_by_arm.csv"))

# Also save a nicely formatted text version
sink(file.path(out_dir, "table1_formatted.txt"))
print(tab1, nonnormal = nonnormal_vars, smd = TRUE)
sink()

cat("\n\nWrote Table 1 to:\n")
cat("  ", file.path(out_dir, "table1_by_arm.csv"), "\n")
cat("  ", file.path(out_dir, "table1_formatted.txt"), "\n")

# ============================================================================
# TABLE 1B: Full CLOVERS cohort (1563) vs biomarker cohort (1340)
#
# Shows that excluded patients (no V1 biomarker) are similar to included
# ============================================================================

cat("\n============================================================\n")
cat("TABLE 1B: Included vs excluded patients\n")
cat("============================================================\n\n")

# Load the full yw.csv to get all 1563 patients
yw_full <- read.csv(YW_CSV, stringsAsFactors = FALSE)
egdt    <- read.csv(EGDT_CSV, stringsAsFactors = FALSE)
other   <- read.csv(OTHER_CSV, stringsAsFactors = FALSE)

# Merge clinical only (no biomarker requirement)
full_clinical <- merge(yw_full, egdt, by = "id", all.x = TRUE)
full_clinical <- merge(full_clinical, other, by = "id", all.x = TRUE)

# Flag who is in the biomarker cohort
full_clinical$in_biomarker_cohort <- ifelse(
  full_clinical$id %in% ds$id, "Included (n=1340)", "Excluded (n=223)"
)

# Use only the clinical variables (no biomarkers for excluded patients)
clinical_vars <- c(cont_vars, bin_vars)

tab1b <- CreateTableOne(
  vars       = clinical_vars[clinical_vars %in% names(full_clinical)],
  strata     = "in_biomarker_cohort",
  data       = full_clinical,
  factorVars = cat_vars[cat_vars %in% names(full_clinical)]
)

tab1b_print <- print(tab1b,
                      nonnormal  = nonnormal_vars,
                      printToggle = FALSE,
                      smd = TRUE,
                      test = TRUE)

cat(capture.output(print(tab1b, nonnormal = nonnormal_vars, smd = TRUE)),
    sep = "\n")

write.csv(tab1b_print,
          file.path(out_dir, "table1_included_vs_excluded.csv"))

cat("\n\nWrote:\n")
cat("  ", file.path(out_dir, "table1_included_vs_excluded.csv"), "\n")
cat("\nDone.\n")