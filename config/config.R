# ============================================================================
# config.R
#
# Central place for every path used across this project. Every script
# sources this file first and refers to these variables instead of hardcoding
# its own paths. If your folder structure changes, this is the ONLY file you
# need to edit.
#
# HOW TO USE: source this from any script with:
#   source(here::here("config", "config.R"))
# or, if not using the `here` package:
#   source("../config/config.R")   # adjust ../ depth based on script location
# ============================================================================

# ---- project root ----
# This should be the absolute path to the top-level project folder, i.e. the
# folder that directly contains data/, cate-repo/, preprocessing/, outputs/,
# config/. Set this once, here, and every other path below is built from it.

PROJECT_ROOT <- "/Users/jaspreetsingh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Clover_Dataset/clover_hte_biomarkers"
# ^ this is pre-filled based on your terminal path. If you move the project
#   folder, update this line to match the new location.

# ---- raw data locations (read-only -- scripts should never write here) ----
DATA_DIR        <- file.path(PROJECT_ROOT, "data", "Data")
CURATED_DIR     <- file.path(DATA_DIR, "Curated datasets")
BIOMARKER_DIR   <- file.path(DATA_DIR, "Biomarker data")
RAW_CSV_DIR     <- file.path(DATA_DIR, "data", "csv")
DOC_DIR         <- file.path(DATA_DIR, "documentation")

YW_CSV          <- file.path(CURATED_DIR, "yw.csv")
EGDT_CSV        <- file.path(CURATED_DIR, "xvars_egdt.csv")
OTHER_CSV       <- file.path(CURATED_DIR, "xvars_other.csv")
BIOMARKER_CSV   <- file.path(BIOMARKER_DIR, "Share.4.29.26.csv")
DATASET_CSV     <- file.path(RAW_CSV_DIR, "DATASET.csv")
DERIVED_CSV     <- file.path(RAW_CSV_DIR, "DERIVED.csv")

# ---- the original CATE repo (read-only -- scripts source() from here) ----
CATE_REPO_DIR        <- file.path(PROJECT_ROOT, "cate-repo")
SCORING_METHODS_R    <- file.path(CATE_REPO_DIR, "scoring.methods_v4.R")
CROSS_VALIDATION_R   <- file.path(CATE_REPO_DIR, "cross.validation.R")
GATES_INFERENCE_R    <- file.path(CATE_REPO_DIR, "gates.inference.R")
POST_PROCESS_R       <- file.path(CATE_REPO_DIR, "post.process_v6.R")
CONFIRM_ANALYSIS_R   <- file.path(CATE_REPO_DIR, "confirm.analysis.R")
MISSFOREST_TRAIN_R   <- file.path(CATE_REPO_DIR, "missforest_train.R")
BALANCED_SPLIT_R     <- file.path(CATE_REPO_DIR, "balanced.split.R")
UPLIFT_FAST_R        <- file.path(CATE_REPO_DIR, "uplift.fast.R")
MICHAEL_Q_R          <- file.path(CATE_REPO_DIR, "michael.q.R")

# ---- outputs ----
# Every script gets its OWN named subfolder under outputs/, named after the
# script itself. This keeps every script's output traceable to exactly what
# produced it, and avoids different scripts silently overwriting each
# other's files.

OUTPUTS_DIR <- file.path(PROJECT_ROOT, "outputs")

make_output_subdir <- function(script_name) {
  # Call this once near the top of any script, passing its own name
  # (without the .R extension), to get a dedicated, auto-created output
  # folder. Returns the path.
  #
  # Example, inside 02b_build_flat_file_and_variable_source_table.R:
  #   out_dir <- make_output_subdir("02b_build_flat_file_and_variable_source_table")
  #   write.csv(flat, file.path(out_dir, "clovers_flat_file.csv"), row.names = FALSE)
  out_dir <- file.path(OUTPUTS_DIR, script_name)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
    cat("Created output folder:", out_dir, "\n")
  }
  out_dir
}

# ---- convenience: a single check that the key raw files actually exist ----
check_config_paths <- function() {
  paths_to_check <- c(
    YW_CSV, EGDT_CSV, OTHER_CSV, BIOMARKER_CSV, DATASET_CSV, DERIVED_CSV,
    SCORING_METHODS_R, CROSS_VALIDATION_R, GATES_INFERENCE_R,
    POST_PROCESS_R, CONFIRM_ANALYSIS_R, MISSFOREST_TRAIN_R
  )
  names(paths_to_check) <- c(
    "YW_CSV","EGDT_CSV","OTHER_CSV","BIOMARKER_CSV","DATASET_CSV","DERIVED_CSV",
    "SCORING_METHODS_R","CROSS_VALIDATION_R","GATES_INFERENCE_R",
    "POST_PROCESS_R","CONFIRM_ANALYSIS_R","MISSFOREST_TRAIN_R"
  )
  ok <- TRUE
  for (i in seq_along(paths_to_check)) {
    exists <- file.exists(paths_to_check[i])
    status <- if (exists) "OK  " else "MISSING"
    cat(sprintf("[%-7s] %-20s %s\n", status, names(paths_to_check)[i], paths_to_check[i]))
    if (!exists) ok <- FALSE
  }
  if (ok) cat("\nAll config paths resolved successfully.\n")
  else cat("\nSome paths are missing -- fix PROJECT_ROOT or your folder layout before continuing.\n")
  invisible(ok)
}