# ============================================================================
# 04_audit_sofa_and_charlson.R
#
# PURPOSE: Independently reconstruct the `sofa`, `sofa_diff`, and `charlson`
# values found in xvars_egdt.csv, working backward from the raw CLOVERS data
# (DERIVED.csv and DATASET.csv), and report any discrepancies found.
#
# This is read-only and exploratory -- it does not modify any data or feed
# into the flat file. It exists to (a) document where these two variables
# actually come from, since neither is computed by any script we have, and
# (b) surface a known bug in `charlson` before it propagates further.
#
# Paths come from config/config.R -- edit PROJECT_ROOT there, not here.
#
# Output (written to outputs/04_audit_sofa_and_charlson/):
#   sofa_audit.csv           -- per-patient comparison, sofa + sofa_diff
#   charlson_audit.csv       -- per-patient comparison, charlson
#   audit_report.txt         -- full console output, saved
# ============================================================================

source(file.path(dirname(getwd()), "config", "config.R"))
out_dir <- make_output_subdir("04_audit_sofa_and_charlson")

report_path <- file.path(out_dir, "audit_report.txt")
report_con <- file(report_path, open = "wt")
sink(report_con, split = TRUE)

cat("============================================================\n")
cat("PART 1: SOFA audit\n")
cat("============================================================\n")
cat("Claim being tested: xvars_egdt.csv$sofa and the sofa_diff outcome in\n")
cat("yw.csv are taken directly from DERIVED.csv's official d_sofa_gcs and\n")
cat("d_sofa_gcs_change fields, with no separate computation.\n\n")

egdt    <- read.csv(EGDT_CSV, stringsAsFactors = FALSE)
yw      <- read.csv(YW_CSV, stringsAsFactors = FALSE)
derived <- read.csv(DERIVED_CSV, stringsAsFactors = FALSE)

derived_baseline <- derived[derived$intervalname == "Baseline", ]
derived_day3     <- derived[derived$intervalname == "Day 3", ]

sofa_check <- merge(
  egdt[, c("id", "sofa")],
  derived_baseline[, c("id", "d_sofa_gcs")],
  by = "id"
)
sofa_check$sofa_match <- sofa_check$sofa == sofa_check$d_sofa_gcs

n_sofa_total   <- nrow(sofa_check)
n_sofa_match   <- sum(sofa_check$sofa_match, na.rm = TRUE)
n_sofa_mismatch <- sum(!sofa_check$sofa_match, na.rm = TRUE)

cat("sofa (baseline covariate) vs DERIVED.csv$d_sofa_gcs:\n")
cat("  Patients compared:", n_sofa_total, "\n")
cat("  Exact matches:    ", n_sofa_match, "\n")
cat("  Mismatches:       ", n_sofa_mismatch, "\n\n")

if (n_sofa_mismatch > 0) {
  cat("MISMATCHES FOUND -- printing all of them:\n")
  print(sofa_check[!sofa_check$sofa_match & !is.na(sofa_check$sofa_match), ])
} else {
  cat("No mismatches. `sofa` matches DERIVED.csv exactly for every patient.\n")
}

# sofa_diff vs d_sofa_gcs_change
sofa_diff_check <- merge(
  yw[, c("id", "sofa_diff")],
  derived_day3[, c("id", "d_sofa_gcs_change")],
  by = "id", all.x = TRUE
)
sofa_diff_check$diff_match <- mapply(function(a, b) {
  if (is.na(a) && is.na(b)) return(TRUE)   # both NA = consistent (no Day 3 assessment)
  if (is.na(a) || is.na(b)) return(FALSE)  # one NA, one not = real discrepancy
  a == b
}, sofa_diff_check$sofa_diff, sofa_diff_check$d_sofa_gcs_change)

n_diff_total    <- nrow(sofa_diff_check)
n_diff_match    <- sum(sofa_diff_check$diff_match)
n_diff_mismatch <- sum(!sofa_diff_check$diff_match)

cat("\nsofa_diff (outcome) vs DERIVED.csv$d_sofa_gcs_change (Day 3 interval):\n")
cat("  Patients compared:", n_diff_total, "\n")
cat("  Consistent (match, or both NA): ", n_diff_match, "\n")
cat("  Real discrepancies:             ", n_diff_mismatch, "\n\n")

if (n_diff_mismatch > 0) {
  cat("DISCREPANCIES FOUND -- printing all of them:\n")
  print(sofa_diff_check[!sofa_diff_check$diff_match, ])
} else {
  cat("No real discrepancies. Every NA in sofa_diff corresponds to a patient\n")
  cat("with no Day 3 SOFA assessment recorded at all (death/discharge/missed\n")
  cat("visit) -- not a derivation error.\n")
}

write.csv(merge(sofa_check, sofa_diff_check, by = "id", all = TRUE),
          file.path(out_dir, "sofa_audit.csv"), row.names = FALSE)

cat("\n============================================================\n")
cat("PART 2: Charlson audit\n")
cat("============================================================\n")
cat("Claim being tested: xvars_egdt.csv$charlson = (sum of charl_* binary\n")
cat("flags from DATASET.csv, standard weights) + (age-bracket points).\n")
cat("Known issue going in: age exactly 80 may get 0 age-points instead of 4.\n\n")

dataset  <- read.csv(DATASET_CSV, stringsAsFactors = FALSE)
baseline <- dataset[dataset$intervalname == "Baseline", ]

age_pts <- function(age) {
  age <- as.numeric(age)
  ifelse(age < 50, 0,
  ifelse(age < 60, 1,
  ifelse(age < 70, 2,
  ifelse(age < 80, 3, 4))))
}

# weights confirmed empirically (see project history): tumor/liver/diabetes/
# kidney severity grades contribute 0 in this dataset's charlson, even though
# the CRF collects them -- confirm with Victor whether that's intentional.
w1 <- c("charl_myoinfarc","charl_congheart","charl_perivasc","charl_cerebvasc",
        "charl_dementia","charl_copd","charl_contis","charl_ulcer")
w2 <- c("charl_hemiplegia","charl_leuk","charl_lymph")

charl_cols_needed <- c(w1, w2)
missing_cols <- setdiff(charl_cols_needed, names(baseline))
if (length(missing_cols) > 0) {
  cat("ERROR: DATASET.csv is missing expected columns:\n")
  print(missing_cols)
  stop("Cannot continue Charlson audit -- check DATASET.csv structure.")
}

baseline$comorbid_pts <- rowSums(sapply(w1, function(c) baseline[[c]] == "Yes"), na.rm = FALSE) +
                          2 * rowSums(sapply(w2, function(c) baseline[[c]] == "Yes"), na.rm = FALSE)

charl_check <- merge(
  egdt[, c("id", "age", "charlson")],
  baseline[, c("id", "comorbid_pts")],
  by = "id"
)
charl_check <- charl_check[!is.na(charl_check$charlson) & charl_check$charlson != "", ]
charl_check$charlson <- as.numeric(charl_check$charlson)
charl_check$age_pts_expected <- age_pts(charl_check$age)
charl_check$computed_total   <- charl_check$comorbid_pts + charl_check$age_pts_expected
charl_check$diff             <- charl_check$computed_total - charl_check$charlson
charl_check$match            <- charl_check$diff == 0

n_charl_total    <- nrow(charl_check)
n_charl_match    <- sum(charl_check$match)
n_charl_mismatch <- sum(!charl_check$match)

cat("charlson vs (comorbidity points + age-bracket points):\n")
cat("  Patients compared:", n_charl_total, "\n")
cat("  Exact matches:    ", n_charl_match, "\n")
cat("  Mismatches:       ", n_charl_mismatch, "\n\n")

if (n_charl_mismatch > 0) {
  cat("MISMATCHES FOUND -- printing all of them:\n")
  mismatches <- charl_check[!charl_check$match, ]
  print(mismatches)

  ages_involved <- sort(unique(mismatches$age))
  cat("\nAges represented among mismatches:", paste(ages_involved, collapse = ", "), "\n")
  if (length(ages_involved) == 1 && ages_involved == 80) {
    cat("--> CONFIRMED: every mismatch is age exactly 80. This matches the\n")
    cat("    known bug where age==80 gets 0 age-points instead of 4.\n")
  } else {
    cat("--> NOTE: mismatches involve ages other than 80 -- this may be a\n")
    cat("    DIFFERENT or ADDITIONAL issue beyond the known age-80 bug.\n")
    cat("    Investigate before assuming the existing explanation covers it.\n")
  }
} else {
  cat("No mismatches. `charlson` matches the reconstructed formula exactly\n")
  cat("for every patient (including age==80 -- the known bug may already be\n")
  cat("fixed in this copy of xvars_egdt.csv).\n")
}

write.csv(charl_check, file.path(out_dir, "charlson_audit.csv"), row.names = FALSE)

cat("\n============================================================\n")
cat("END OF AUDIT\n")
cat("============================================================\n")

sink()
close(report_con)
cat("\nFull report saved to:\n", report_path, "\n")
cat("Per-patient detail saved to sofa_audit.csv and charlson_audit.csv in the\n")
cat("same folder.\n")