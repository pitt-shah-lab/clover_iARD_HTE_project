# ============================================================================
# 03_split_derivation_validation.R
#
# PURPOSE (Victor's stage 2): split the flat file into derivation and
# validation halves.
#
# Method mirrors ADRENAL_internal_split.R exactly: a RANDOM split, done
# SEPARATELY within each treatment arm (w==0 and w==1), so the overall 50/50
# derivation/validation ratio is preserved within each arm. This is plain
# stratified-random -- not the covariate-balanced SIDES allocation procedure
# in cate-repo/balanced.split.R, which exists but was not what Victor used
# for ADRENAL. Confirmed against his own script before writing this.
#
# Paths come from config/config.R -- edit PROJECT_ROOT there, not here.
#
# Input:  outputs/02_build_flat_file/clovers_flat_file.csv
# Output (written to outputs/03_split_derivation_validation/):
#   ids_internal_list.RDS   -- list(ids_der = ..., ids_val = ...)
#   flat_file_der.csv       -- derivation half, full columns
#   flat_file_val.csv       -- validation half, full columns
# ============================================================================

source(file.path(dirname(getwd()), "config", "config.R"))
out_dir <- make_output_subdir("03_split_derivation_validation")

flat_path <- file.path(OUTPUTS_DIR, "02_build_flat_file", "clovers_flat_file.csv")
if (!file.exists(flat_path)) {
  stop("Flat file not found at ", flat_path, " -- run 02_build_flat_file.R first.")
}
ds0 <- read.csv(flat_path, stringsAsFactors = FALSE)

cat("Loaded flat file:", nrow(ds0), "patients,", ncol(ds0), "columns\n\n")

# ============================================================================
# Split randomly into derivation and validation, stratified by treatment (w)
# -- same approach as ADRENAL_internal_split.R, applied here to `id`/`w`.
# ============================================================================

SPLIT_SEED <- 03202026   # matching Victor's seed convention (date-coded);
                          # change if you want an independent split, but keep
                          # it FIXED and documented once chosen.
set.seed(SPLIT_SEED)

ids_w0_der <- sample(ds0[ds0$w == 0, "id"], nrow(ds0[ds0$w == 0, ]) / 2, replace = FALSE)
ids_w1_der <- sample(ds0[ds0$w == 1, "id"], nrow(ds0[ds0$w == 1, ]) / 2, replace = FALSE)
ids_w0_val <- ds0[ds0$w == 0 & !ds0$id %in% ids_w0_der, "id"]
ids_w1_val <- ds0[ds0$w == 1 & !ds0$id %in% ids_w1_der, "id"]

ids_internal_list <- list(
  "ids_der" = c(ids_w0_der, ids_w1_der),
  "ids_val" = c(ids_w0_val, ids_w1_val)
)

saveRDS(ids_internal_list, file.path(out_dir, "ids_internal_list.RDS"))

der <- ds0[ds0$id %in% ids_internal_list$ids_der, ]
val <- ds0[ds0$id %in% ids_internal_list$ids_val, ]

write.csv(der, file.path(out_dir, "flat_file_der.csv"), row.names = FALSE)
write.csv(val, file.path(out_dir, "flat_file_val.csv"), row.names = FALSE)

cat("Derivation set: ", nrow(der), "patients (w=0:", sum(der$w == 0), ", w=1:", sum(der$w == 1), ")\n")
cat("Validation set: ", nrow(val), "patients (w=0:", sum(val$w == 0), ", w=1:", sum(val$w == 1), ")\n")
cat("\nWrote:\n  ", file.path(out_dir, "ids_internal_list.RDS"), "\n")
cat("  ", file.path(out_dir, "flat_file_der.csv"), "\n")
cat("  ", file.path(out_dir, "flat_file_val.csv"), "\n")

cat("\nNext: run 03b_compare_der_val.R to check covariate balance between\n")
cat("the two halves.\n")