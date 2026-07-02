# ============================================================================
# 04_impute.R
#
# PURPOSE (Victor's stage 3): impute missing covariate values using the
# missForest-based tool in cate-repo/missforest_train.R.
#
# Method mirrors ADRENAL_internal_impute.R exactly: the imputation forests
# are TRAINED ONLY on the derivation set. That same trained forest is then
# applied to the validation set (impute_test_with_rfs) -- validation data
# never influences how its own missing values get filled in. This is the
# correct way to avoid leaking validation-set information into derivation,
# matching what Victor did for ADRENAL.
#
# Difference from his script: our flat file has TWO outcome columns
# (inhosp90, sofa_diff) instead of his single `y`, so the "first 3 columns"
# he excludes by position (id, y, w) becomes "these 4 named columns" for us
# (id, inhosp90, sofa_diff, w) -- selected by NAME here, not position, so a
# future column reorder won't silently break this.
#
# Paths come from config/config.R -- edit PROJECT_ROOT there, not here.
#
# Input:
#   outputs/03_split_derivation_validation/flat_file_der.csv
#   outputs/03_split_derivation_validation/flat_file_val.csv
# Output (written to outputs/04_impute/):
#   dsi_der.RDS, dsi_der.csv   -- imputed derivation set
#   dsi_val.RDS, dsi_val.csv   -- imputed validation set
#   dsi_all.RDS, dsi_all.csv   -- both combined, the "dsi0" needed downstream
# ============================================================================

source(file.path(dirname(getwd()), "config", "config.R"))
out_dir <- make_output_subdir("04_impute")

source(MISSFOREST_TRAIN_R)

split_dir <- file.path(OUTPUTS_DIR, "03_split_derivation_validation")
der_path <- file.path(split_dir, "flat_file_der.csv")
val_path <- file.path(split_dir, "flat_file_val.csv")

if (!file.exists(der_path) || !file.exists(val_path)) {
  stop("Derivation/validation files not found in ", split_dir,
       " -- run 03_split_derivation_validation.R first.")
}

der <- read.csv(der_path, stringsAsFactors = FALSE)
val <- read.csv(val_path, stringsAsFactors = FALSE)

cat("Derivation set:", nrow(der), "patients\n")
cat("Validation set:", nrow(val), "patients\n\n")

# ---- non-covariate columns, selected by NAME (not position) ----
id_outcome_treat_cols <- c("id", "inhosp90", "sofa_diff", "w")

missing_check <- setdiff(id_outcome_treat_cols, names(der))
if (length(missing_check) > 0) {
  stop("Expected columns not found in the flat file: ",
       paste(missing_check, collapse = ", "))
}

covariate_cols <- setdiff(names(der), id_outcome_treat_cols)
cat("Covariate columns to impute:", length(covariate_cols), "\n\n")

# ============================================================================
# Impute: fit on derivation only, apply that same fit to validation
# ============================================================================

set.seed(03202026)   # matching Victor's seed convention

cat("Training imputation forests on the derivation set...\n")
impute_train <- impute_train_and_save_rfs(der[, covariate_cols])

dsi_der <- cbind(der[, id_outcome_treat_cols], impute_train$imputed_data)

cat("Applying the trained forests to the validation set...\n")
xi_val <- impute_test_with_rfs(impute_train$final_rfs, val[, covariate_cols])

dsi_val <- cbind(val[, id_outcome_treat_cols], xi_val)

# ---- combined object, needed by the derivation/CATE script downstream ----
# (this is the "dsi0" gap identified earlier -- ADRENAL_internal_derivation_v2.R
# expects ONE combined file with both der and val rows, re-split internally
# using ids_internal_list.RDS. We build that combined object explicitly here
# so the gap doesn't propagate into our own pipeline.)
dsi_all <- rbind(dsi_der, dsi_val)

# ---- sanity check: no missing values remain in either imputed set ----
n_missing_der <- sum(is.na(dsi_der[, covariate_cols]))
n_missing_val <- sum(is.na(dsi_val[, covariate_cols]))
cat("\nRemaining missing covariate cells after imputation:\n")
cat("  Derivation:", n_missing_der, "\n")
cat("  Validation:", n_missing_val, "\n")
if (n_missing_der > 0 || n_missing_val > 0) {
  cat("WARNING: imputation did not resolve all missing values -- investigate\n")
  cat("before proceeding to CATE modeling.\n")
} else {
  cat("All covariate missingness resolved. Ready for CATE modeling.\n")
}

# ---- write outputs ----
saveRDS(dsi_der, file.path(out_dir, "dsi_der.RDS"))
saveRDS(dsi_val, file.path(out_dir, "dsi_val.RDS"))
saveRDS(dsi_all, file.path(out_dir, "dsi_all.RDS"))
write.csv(dsi_der, file.path(out_dir, "dsi_der.csv"), row.names = FALSE)
write.csv(dsi_val, file.path(out_dir, "dsi_val.csv"), row.names = FALSE)
write.csv(dsi_all, file.path(out_dir, "dsi_all.csv"), row.names = FALSE)

cat("\nWrote imputed datasets to:\n  ", out_dir, "\n")

# ---- flag a known quirk: missforest_train.R does not distinguish binary
# from continuous covariates, so imputed values for binary columns (e.g.
# kidney, mv, chf) may come back as fractional values (like 0.23) rather
# than clean 0/1, since randomForest() runs in regression mode for any
# numeric y. This is not necessarily wrong for causal_forest()/bart() --
# they can take a continuous value between 0 and 1 as a soft indicator --
# but it's worth knowing about and confirming with Victor whether his
# ADRENAL pipeline rounds these or leaves them fractional too. ----
binary_covariate_cols <- c("site_lung","site_abdom","site_urine","mv","vaso",
                            "ards","dial","chf","copd","liver","kidney")
binary_covariate_cols <- intersect(binary_covariate_cols, covariate_cols)

cat("\n============================================================\n")
cat("CHECK: fractional values in binary covariate columns after imputation\n")
cat("============================================================\n")
for (col in binary_covariate_cols) {
  vals_der <- dsi_der[[col]]
  n_frac_der <- sum(vals_der != 0 & vals_der != 1)
  vals_val <- dsi_val[[col]]
  n_frac_val <- sum(vals_val != 0 & vals_val != 1)
  if (n_frac_der > 0 || n_frac_val > 0) {
    cat(sprintf("  %-14s der: %d fractional values | val: %d fractional values\n",
                col, n_frac_der, n_frac_val))
  }
}
cat("(if any rows printed above, those are imputed values that are not\n")
cat(" clean 0/1 -- confirm with Victor whether to round these before modeling)\n")
cat("\nNext: build the CATE driver script (mirroring\n")
cat("ADRENAL_internal_derivation_v2.R), reading dsi_all.RDS and re-splitting\n")
cat("by ids_internal_list.RDS the same way his script does.\n")