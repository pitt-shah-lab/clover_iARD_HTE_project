# ============================================================================
# 01b_explore_all_datasets.R
#
# PURPOSE: Before deciding what goes in the flat file, see what's actually
# available. This scans every CSV and XLSX file under data/ (recursively --
# it doesn't matter exactly how deep your folders are nested, and it covers
# files sitting directly in data/ as well as anything nested under
# data/Data/) and prints, for each file:
#   - its path
#   - its shape (rows x columns)
#   - every column name
#   - a few sample values per column (first non-missing value found)
#
# This is read-only and exploratory -- it doesn't build anything. The goal
# is to give you and Victor a single printout you can scan to manually
# decide which columns belong in the flat file.
#
# Paths come from config/config.R -- edit PROJECT_ROOT there, not here.
#
# Output: printed to console AND saved as a text file in this script's own
# output subfolder, so you can scroll back through it later or send it to
# Victor without re-running anything.
# ============================================================================

source(file.path(dirname(getwd()), "config", "config.R"))
out_dir <- make_output_subdir("01b_explore_all_datasets")

# ---- find every CSV and XLSX under the whole data/ folder, however deep ----
# NOTE: scans the parent of DATA_DIR (i.e. data/, not just data/Data/) so it
# also picks up files like clovers_forHTE_flat.csv or any .xlsx that live
# directly in data/ rather than nested under data/Data/.
SCAN_ROOT <- dirname(DATA_DIR)
csv_files  <- list.files(SCAN_ROOT, pattern = "\\.csv$",  recursive = TRUE,
                          full.names = TRUE, ignore.case = TRUE)
xlsx_files <- list.files(SCAN_ROOT, pattern = "\\.xlsx$", recursive = TRUE,
                          full.names = TRUE, ignore.case = TRUE)

cat("Found", length(csv_files), "CSV files and", length(xlsx_files),
    "XLSX files under:", SCAN_ROOT, "\n\n")

if (length(csv_files) == 0 && length(xlsx_files) == 0) {
  stop("No CSV or XLSX files found. Check PROJECT_ROOT and DATA_DIR in config.R.")
}

# ---- helper: get a few real sample values for one column ----
sample_values <- function(col, n = 3) {
  non_na <- col[!is.na(col) & col != "" & col != "NA"]
  if (length(non_na) == 0) return("(all missing)")
  paste(utils::head(unique(non_na), n), collapse = ", ")
}

# ---- open a connection that writes to BOTH console and a text file ----
report_path <- file.path(out_dir, "dataset_exploration_report.txt")
report_con <- file(report_path, open = "wt")
sink(report_con, split = TRUE)   # split=TRUE -> still prints to console too

cat("============================================================\n")
cat("DATASET EXPLORATION REPORT\n")
cat("Generated from:", SCAN_ROOT, "\n")
cat("============================================================\n\n")

all_files <- c(csv_files, xlsx_files)

for (f in all_files) {
  cat("\n------------------------------------------------------------\n")
  cat("FILE:", f, "\n")

  is_xlsx <- grepl("\\.xlsx$", f, ignore.case = TRUE)

  df <- tryCatch(
    if (is_xlsx) {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("readxl package not installed -- run install.packages('readxl') first")
      }
      as.data.frame(readxl::read_excel(f, col_types = "text"))
    } else {
      read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    },
    error = function(e) {
      cat("  COULD NOT READ THIS FILE:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(df)) next

  cat("SHAPE:", nrow(df), "rows x", ncol(df), "columns\n")
  cat("------------------------------------------------------------\n")
  cat(sprintf("%-4s %-30s %-10s %-10s %s\n",
              "#", "column", "n_missing", "n_unique", "sample values"))

  for (i in seq_along(df)) {
    col <- df[[i]]
    n_missing <- sum(is.na(col) | col == "" | col == "NA")
    n_unique  <- length(unique(col[!is.na(col)]))
    cat(sprintf("%-4d %-30s %-10d %-10d %s\n",
                i, names(df)[i], n_missing, n_unique, sample_values(col)))
  }
}

cat("\n============================================================\n")
cat("END OF REPORT --", length(all_files), "files scanned\n")
cat("============================================================\n")

sink()
close(report_con)

cat("\n\nFull report also saved to:\n", report_path, "\n")
cat("Open that file (or scroll up) to manually decide what belongs in the flat file.\n")
