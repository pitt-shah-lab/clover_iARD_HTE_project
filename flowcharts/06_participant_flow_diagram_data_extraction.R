# ============================================================================
# 06_participant_flow_extract.R
#
# PURPOSE: extract the real patient counts at each stage of inclusion in
# this project's analytic cohort, and emit a Mermaid flowchart syntax block
# you can paste into mermaid.live, a GitHub markdown file, or any tool that
# renders Mermaid syntax.
#
# IMPORTANT HONESTY NOTE: this can only produce counts for the part of the
# flow that's downstream of randomization, because the files in this project
# represent the analytic cohort -- they don't contain the screening logs or
# excluded-pre-randomization patients. Upstream numbers (screened, refused
# consent, excluded for protocol reasons, etc.) must be filled in by hand
# from the published CLOVERS protocol or main-results paper. Placeholders
# are left in the Mermaid output for those numbers.
#
# Paths come from config/config.R.
#
# Output (written to outputs/06_participant_flow_extract/):
#   participant_flow_counts.csv  -- table of every count
#   participant_flow.mmd         -- Mermaid syntax, ready to render
#   participant_flow_report.txt  -- full console output captured
# ============================================================================

source(file.path(dirname(getwd()), "config", "config.R"))
out_dir <- make_output_subdir("06_participant_flow_extract")

report_path <- file.path(out_dir, "participant_flow_report.txt")
report_con  <- file(report_path, open = "wt")
sink(report_con, split = TRUE)

# ---- count each cohort, source by source ---------------------------------

yw    <- read.csv(YW_CSV)
egdt  <- read.csv(EGDT_CSV)
other <- read.csv(OTHER_CSV)
bio   <- read.csv(BIOMARKER_CSV)

strip_index_col <- function(df) {
  if (names(df)[1] == "" || names(df)[1] == "X") df <- df[, -1]
  df
}
yw    <- strip_index_col(yw)
egdt  <- strip_index_col(egdt)
other <- strip_index_col(other)

# Analytic cohort (post-randomization): everyone in yw.csv with both treatment
# arm and covariates -- this is what Victor calls the "curated cohort"
n_analytic <- nrow(yw)
n_analytic_w0 <- sum(yw$w == 0)
n_analytic_w1 <- sum(yw$w == 1)

# Have V1 biomarker draw
bio_v1_ids <- unique(bio$id[bio$visit == "V1"])
n_with_bio <- sum(yw$id %in% bio_v1_ids)
n_no_bio   <- n_analytic - n_with_bio

# After inner-join with biomarker file (the cohort we actually model)
biomarker_cohort <- yw[yw$id %in% bio_v1_ids, ]
n_biomarker_w0 <- sum(biomarker_cohort$w == 0)
n_biomarker_w1 <- sum(biomarker_cohort$w == 1)

# Derivation / validation split (from script 03's output, if available)
split_path <- file.path(OUTPUTS_DIR, "03_split_derivation_validation",
                         "ids_internal_list.RDS")
if (file.exists(split_path)) {
  ids <- readRDS(split_path)
  n_der <- length(ids$ids_der)
  n_val <- length(ids$ids_val)
  n_der_w0 <- sum(biomarker_cohort$id %in% ids$ids_der & biomarker_cohort$w == 0)
  n_der_w1 <- sum(biomarker_cohort$id %in% ids$ids_der & biomarker_cohort$w == 1)
  n_val_w0 <- sum(biomarker_cohort$id %in% ids$ids_val & biomarker_cohort$w == 0)
  n_val_w1 <- sum(biomarker_cohort$id %in% ids$ids_val & biomarker_cohort$w == 1)
} else {
  n_der <- n_val <- n_der_w0 <- n_der_w1 <- n_val_w0 <- n_val_w1 <- NA
}

# Outcome event counts
n_events_total      <- sum(yw$inhosp90 == 1, na.rm = TRUE)
n_events_biomarker  <- sum(biomarker_cohort$inhosp90 == 1, na.rm = TRUE)

cat("============================================================\n")
cat("PARTICIPANT FLOW COUNTS\n")
cat("============================================================\n\n")
cat("Analytic cohort (yw.csv):              ", n_analytic, "patients\n")
cat("  Arm w=0 (restrictive fluid):         ", n_analytic_w0, "\n")
cat("  Arm w=1 (liberal fluid):             ", n_analytic_w1, "\n")
cat("  Outcome events (inhosp90 == 1):      ", n_events_total,
    sprintf("(%.1f%%)\n", 100 * n_events_total / n_analytic))
cat("\n")
cat("With baseline (V1) biomarker draw:     ", n_with_bio, "patients\n")
cat("Excluded for no biomarker draw:        ", n_no_bio, "patients\n")
cat("  Arm w=0 in biomarker cohort:         ", n_biomarker_w0, "\n")
cat("  Arm w=1 in biomarker cohort:         ", n_biomarker_w1, "\n")
cat("  Outcome events in biomarker cohort:  ", n_events_biomarker,
    sprintf("(%.1f%%)\n", 100 * n_events_biomarker / n_with_bio))
cat("\n")
if (!is.na(n_der)) {
  cat("Derivation set (50%):                  ", n_der, "patients\n")
  cat("  w=0:", n_der_w0, " | w=1:", n_der_w1, "\n")
  cat("Validation set (50%):                  ", n_val, "patients\n")
  cat("  w=0:", n_val_w0, " | w=1:", n_val_w1, "\n")
} else {
  cat("Derivation/validation split:           (not yet built -- run 03_split first)\n")
}

# ---- write counts to CSV ---------------------------------------------------
counts_df <- data.frame(
  stage = c("analytic_cohort_total",
            "analytic_cohort_w0", "analytic_cohort_w1",
            "analytic_cohort_events",
            "with_v1_biomarker", "excluded_no_biomarker",
            "biomarker_cohort_w0", "biomarker_cohort_w1",
            "biomarker_cohort_events",
            "derivation_n", "derivation_w0", "derivation_w1",
            "validation_n", "validation_w0", "validation_w1"),
  count = c(n_analytic, n_analytic_w0, n_analytic_w1, n_events_total,
            n_with_bio, n_no_bio, n_biomarker_w0, n_biomarker_w1, n_events_biomarker,
            n_der, n_der_w0, n_der_w1, n_val, n_val_w0, n_val_w1),
  stringsAsFactors = FALSE
)
write.csv(counts_df, file.path(out_dir, "participant_flow_counts.csv"),
          row.names = FALSE)

# ---- emit Mermaid syntax --------------------------------------------------
# Upstream nodes are placeholders -- fill in from published CLOVERS trial
mermaid <- c(
  "flowchart TD",
  "    %% UPSTREAM NODES -- FILL IN FROM PUBLISHED CLOVERS PROTOCOL/RESULTS",
  "    A[Patients screened<br/>n = ???]",
  "    A --> B[Excluded pre-randomization<br/>n = ???<br/>Reasons: ???]",
  "    A --> C[Randomized<br/>n = ???]",
  "",
  "    %% DOWNSTREAM FROM THE ANALYTIC COHORT FILES -- REAL NUMBERS BELOW",
  sprintf("    C --> D[Analytic cohort<br/>n = %d<br/>w=0: %d &nbsp;w=1: %d<br/>Events: %d]",
          n_analytic, n_analytic_w0, n_analytic_w1, n_events_total),
  sprintf("    D --> E[Excluded for no baseline biomarker draw<br/>n = %d]",
          n_no_bio),
  sprintf("    D --> F[Biomarker analytic cohort<br/>n = %d<br/>w=0: %d &nbsp;w=1: %d<br/>Events: %d]",
          n_with_bio, n_biomarker_w0, n_biomarker_w1, n_events_biomarker)
)

if (!is.na(n_der)) {
  mermaid <- c(
    mermaid,
    sprintf("    F --> G[Derivation set<br/>n = %d<br/>w=0: %d &nbsp;w=1: %d]",
            n_der, n_der_w0, n_der_w1),
    sprintf("    F --> H[Validation set<br/>n = %d<br/>w=0: %d &nbsp;w=1: %d]",
            n_val, n_val_w0, n_val_w1),
    "    G --> I[Imputation fit on derivation]",
    "    H --> I",
    "    I --> J[CATE modeling: clinical-only vs clinical+biomarker]"
  )
}

mermaid <- c(mermaid,
  "",
  "    %% styling -- placeholder nodes shown in red so they're hard to miss",
  "    style A fill:#fee,stroke:#c00",
  "    style B fill:#fee,stroke:#c00",
  "    style C fill:#fee,stroke:#c00"
)

mermaid_text <- paste(mermaid, collapse = "\n")
mermaid_path <- file.path(out_dir, "participant_flow.mmd")
writeLines(mermaid, mermaid_path)

cat("\n============================================================\n")
cat("MERMAID FLOWCHART SYNTAX\n")
cat("============================================================\n")
cat("Saved to:", mermaid_path, "\n\n")
cat("Paste the contents below into https://mermaid.live to render it,\n")
cat("or embed in a markdown file inside a ```mermaid ... ``` block.\n")
cat("The red-shaded nodes at the top are placeholders -- replace the ???\n")
cat("values with numbers from the published CLOVERS trial paper.\n\n")
cat(mermaid_text)
cat("\n")

sink()
close(report_con)
cat("\nFull report saved to:\n  ", report_path, "\n")