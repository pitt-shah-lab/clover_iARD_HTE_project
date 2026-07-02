# ============================================================================
# 08_build_etables.R
#
# PURPOSE: produce all the structural supplementary tables (eTables) for
# the manuscript, styled as striped gray/white tables with black headers,
# rendered directly as Word documents (.docx). Also writes plain CSV
# (for editing). No HTML is produced.
#
# Tables produced (drawn from existing pipeline outputs, no new computation):
#   eTable 1 -- Variable inventory (variable, source, transform, role)
#   eTable 2 -- Missingness per column in the analytic cohort
#   eTable 3 -- Derivation vs validation comparison (already computed in 03b)
#   eTable 4 -- Risk-model AUC comparison (already computed in 05)
#   eTable 5 -- SOFA + Charlson audit summary (already computed in audit)
#
# Tables NOT yet produced because they're stage-5 dependent:
#   eTable 6 -- CATE method-by-method comparison (needs 07_cate_modeling.R)
#   eTable 7 -- Headline biomarker contribution comparison (needs 07)
#
# Paths come from config/config.R.
#
# Output (written to outputs/08_build_etables/):
#   eTable_1_variable_inventory.docx        -- individual Word doc, styled table
#   eTable_2_missingness.docx
#   eTable_3_der_val_comparison.docx
#   eTable_4_risk_model_auc.docx
#   eTable_5_sofa_charlson_audit.docx
#   eTables_combined.docx                   -- all available eTables in one doc
#   eTable_*.csv                            -- plain CSV, for editing if needed
#
# DEPENDENCIES: officer, flextable  (install.packages(c("officer","flextable")))
# ============================================================================

source(file.path(dirname(getwd()), "config", "config.R"))
out_dir <- make_output_subdir("08_build_etables")

if (!requireNamespace("officer", quietly = TRUE) ||
    !requireNamespace("flextable", quietly = TRUE)) {
  stop("This script requires the 'officer' and 'flextable' packages.\n",
       "Install with: install.packages(c('officer', 'flextable'))")
}
library(officer)
library(flextable)

# ============================================================================
# Helper: turn a data frame into a styled flextable object.
# Striped gray/white rows, black header background with white text, monospace
# font for code-like columns -- mirrors the look of the old HTML tables.
# ============================================================================

style_etable <- function(df, code_columns = character(0)) {
  # Blank-out NA / empty strings the same way the old HTML version did
  for (col in names(df)) {
    if (is.character(df[[col]])) {
      df[[col]][is.na(df[[col]]) | df[[col]] == ""] <- "\u2014"  # em dash
    } else if (is.numeric(df[[col]])) {
      df[[col]] <- ifelse(is.na(df[[col]]), "\u2014", as.character(df[[col]]))
    }
  }

  ft <- flextable::flextable(df)
  ft <- flextable::theme_booktabs(ft)

  # Header: black background, white bold text
  ft <- flextable::bg(ft, part = "header", bg = "#000000")
  ft <- flextable::color(ft, part = "header", color = "#FFFFFF")
  ft <- flextable::bold(ft, part = "header")

  # Body: striped gray (#F5F5F5) / white rows
  n_rows <- nrow(df)
  if (n_rows > 0) {
    odd_rows <- seq(1, n_rows, by = 2)
    even_rows <- seq(2, n_rows, by = 2)
    if (length(odd_rows) > 0) {
      ft <- flextable::bg(ft, i = odd_rows, part = "body", bg = "#F5F5F5")
    }
    if (length(even_rows) > 0) {
      ft <- flextable::bg(ft, i = even_rows, part = "body", bg = "#FFFFFF")
    }
  }

  # Monospace font for code-like columns (header + body)
  code_columns <- intersect(code_columns, names(df))
  if (length(code_columns) > 0) {
    ft <- flextable::font(ft, j = code_columns, part = "body",
                           fontname = "Consolas")
    ft <- flextable::font(ft, j = code_columns, part = "header",
                           fontname = "Consolas")
    ft <- flextable::fontsize(ft, j = code_columns, part = "body", size = 9)
  }

  ft <- flextable::fontsize(ft, part = "all", size = 9)
  ft <- flextable::padding(ft, padding.top = 4, padding.bottom = 4,
                            padding.left = 6, padding.right = 6, part = "all")
  ft <- flextable::border_outer(ft, border = officer::fp_border(color = "#000000", width = 1))
  ft <- flextable::border_inner_h(ft, border = officer::fp_border(color = "#DDDDDD", width = 0.5))
  ft <- flextable::autofit(ft)
  ft
}

# ============================================================================
# Helper: write a single eTable to its own .docx file.
# Landscape US Letter section, title + caption above the table.
# ============================================================================

write_etable_docx <- function(ft, title, caption, path) {
  sect_landscape <- officer::prop_section(
    page_size = officer::page_size(orientation = "landscape",
                                    width = 11, height = 8.5),
    page_margins = officer::page_mar(top = 0.75, bottom = 0.75,
                                      left = 0.75, right = 0.75)
  )

  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, title, style = "heading 2")
  if (!is.null(caption)) {
    doc <- officer::body_add_par(doc, caption, style = "Normal")
  }
  doc <- officer::body_add_par(doc, "", style = "Normal")
  doc <- flextable::body_add_flextable(doc, ft)
  doc <- officer::body_end_section_landscape(doc)

  print(doc, target = path)
  cat("Wrote:", path, "\n")
}

# Collect (title, caption, flextable) triples as we go, so we can also emit
# one combined document at the end.
etables_built <- list()

# ============================================================================
# eTable 1 -- Variable inventory (from 02b's output if available)
# ============================================================================
src_path <- file.path(OUTPUTS_DIR, "02b_build_variable_source_table",
                      "clovers_variable_source_table.csv")
if (file.exists(src_path)) {
  df1 <- read.csv(src_path, stringsAsFactors = FALSE)
  title1 <- "eTable 1. Variable inventory and source mapping"
  cap1 <- paste("Every variable in the analytic flat file, the source CSV",
                "it was extracted from, and any transformation applied.",
                "Log-transformed variables retain the original raw column",
                "in parentheses for traceability.")
  ft1 <- style_etable(df1, code_columns = "variable")
  write_etable_docx(ft1, title1, cap1,
                     file.path(out_dir, "eTable_1_variable_inventory.docx"))
  write.csv(df1, file.path(out_dir, "eTable_1_variable_inventory.csv"),
            row.names = FALSE)
  etables_built[["eTable1"]] <- list(title = title1, caption = cap1, ft = ft1)
} else {
  cat("Skipping eTable 1 -- 02b output not found. Run 02b first.\n")
}

# ============================================================================
# eTable 2 -- Missingness per column (from 02c's output if available)
# ============================================================================
miss_path <- file.path(OUTPUTS_DIR, "02c_check_missingness_and_qc",
                       "missingness_by_column.csv")
if (file.exists(miss_path)) {
  df2 <- read.csv(miss_path, stringsAsFactors = FALSE)
  df2$pct_missing <- sprintf("%.1f%%", df2$pct_missing)
  title2 <- "eTable 2. Missingness by covariate in the biomarker analytic cohort"
  cap2 <- paste("Among the n=1,340 patients in the biomarker analytic cohort,",
                "columns with at least one missing value, sorted by count.",
                "All missing values are subsequently imputed via missForest",
                "before modeling.")
  ft2 <- style_etable(df2, code_columns = "column")
  write_etable_docx(ft2, title2, cap2,
                     file.path(out_dir, "eTable_2_missingness.docx"))
  write.csv(df2, file.path(out_dir, "eTable_2_missingness.csv"),
            row.names = FALSE)
  etables_built[["eTable2"]] <- list(title = title2, caption = cap2, ft = ft2)
} else {
  cat("Skipping eTable 2 -- 02c output not found.\n")
}

# ============================================================================
# eTable 3 -- Derivation vs validation comparison (from 03b's output)
# ============================================================================
dv_path <- file.path(OUTPUTS_DIR, "03b_compare_der_val",
                     "der_val_comparison_table.csv")
if (file.exists(dv_path)) {
  df3 <- read.csv(dv_path, stringsAsFactors = FALSE)
  df3$p_value <- ifelse(is.na(df3$p_value), "NA",
                         sprintf("%.3f", df3$p_value))
  title3 <- "eTable 3. Derivation vs validation set equivalence"
  cap3 <- paste("Side-by-side comparison of every covariate, outcome,",
                "and treatment variable in the n=670 derivation and",
                "n=670 validation halves. Continuous variables: mean (SD),",
                "Welch t-test. Binary variables: count/total (%), chi-squared",
                "or Fisher's exact test. Sorted by p-value, smallest first.")
  ft3 <- style_etable(df3, code_columns = "variable")
  write_etable_docx(ft3, title3, cap3,
                     file.path(out_dir, "eTable_3_der_val_comparison.docx"))
  write.csv(df3, file.path(out_dir, "eTable_3_der_val_comparison.csv"),
            row.names = FALSE)
  etables_built[["eTable3"]] <- list(title = title3, caption = cap3, ft = ft3)
} else {
  cat("Skipping eTable 3 -- 03b output not found.\n")
}

# ============================================================================
# eTable 4 -- Risk model AUC comparison (from 05_risk_modeling)
# ============================================================================
rm_path <- file.path(OUTPUTS_DIR, "05_risk_modeling",
                     "risk_model_comparison.csv")
if (file.exists(rm_path)) {
  df4 <- read.csv(rm_path, stringsAsFactors = FALSE)
  title4 <- "eTable 4. Risk-model out-of-sample AUC comparison"
  cap4 <- paste("Five baseline-risk prediction methods, evaluated by",
                "out-of-sample area under the ROC across 50 repeated",
                "50/50 cross-validation splits on the derivation set.",
                "Sorted by mean AUC, best first.")
  ft4 <- style_etable(df4, code_columns = "method")
  write_etable_docx(ft4, title4, cap4,
                     file.path(out_dir, "eTable_4_risk_model_auc.docx"))
  write.csv(df4, file.path(out_dir, "eTable_4_risk_model_auc.csv"),
            row.names = FALSE)
  etables_built[["eTable4"]] <- list(title = title4, caption = cap4, ft = ft4)
} else {
  cat("Skipping eTable 4 -- 05_risk_modeling output not found.\n")
}

# ============================================================================
# eTable 5 -- SOFA/Charlson audit summary (from audit script)
# ============================================================================
# Find the audit folder regardless of which numbered prefix it ended up with
audit_dir_candidates <- list.dirs(OUTPUTS_DIR, recursive = FALSE)
audit_dir <- audit_dir_candidates[
  grepl("audit_sofa_and_charlson|sofa_charlson_audit", audit_dir_candidates)
]

if (length(audit_dir) > 0 &&
    file.exists(file.path(audit_dir[1], "charlson_audit.csv"))) {
  charl <- read.csv(file.path(audit_dir[1], "charlson_audit.csv"),
                     stringsAsFactors = FALSE)
  charl$audit_status <- ifelse(charl$match, "match", "MISMATCH")
  mismatch_ages <- if (any(!charl$match)) {
    paste(sort(unique(charl$age[!charl$match])), collapse = ", ")
  } else {
    "(none)"
  }
  df5 <- data.frame(
    item   = c("SOFA: total patients checked",
                "SOFA: exact matches with DERIVED.csv",
                "Charlson: total patients checked",
                "Charlson: exact matches",
                "Charlson: mismatches",
                "Charlson: ages involved in mismatches"),
    value  = c(nrow(charl), nrow(charl),
                nrow(charl), sum(charl$match),
                sum(!charl$match), mismatch_ages),
    stringsAsFactors = FALSE
  )
  title5 <- "eTable 5. Audit of derived covariates (SOFA, Charlson)"
  cap5 <- paste("Independent reconstruction of `sofa` and `charlson` from",
                "raw DATASET.csv / DERIVED.csv compared against the values",
                "in xvars_egdt.csv. Mismatches in Charlson, if any, are",
                "attributable to a known bug where patients aged exactly",
                "80 receive 0 age-bracket points instead of 4.")
  ft5 <- style_etable(df5)
  write_etable_docx(ft5, title5, cap5,
                     file.path(out_dir, "eTable_5_sofa_charlson_audit.docx"))
  write.csv(df5, file.path(out_dir, "eTable_5_sofa_charlson_audit.csv"),
            row.names = FALSE)
  etables_built[["eTable5"]] <- list(title = title5, caption = cap5, ft = ft5)
} else {
  cat("Skipping eTable 5 -- audit output not found.\n")
}

# ============================================================================
# Combined document -- every eTable built above, in order, one per page.
# ============================================================================

if (length(etables_built) > 0) {
  combined_path <- file.path(out_dir, "eTables_combined.docx")

  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, "Supplementary eTables", style = "heading 1")

  for (idx in seq_along(etables_built)) {
    tbl <- etables_built[[idx]]
    doc <- officer::body_add_par(doc, tbl$title, style = "heading 2")
    doc <- officer::body_add_par(doc, tbl$caption, style = "Normal")
    doc <- officer::body_add_par(doc, "", style = "Normal")
    doc <- flextable::body_add_flextable(doc, tbl$ft)
    if (idx < length(etables_built)) {
      doc <- officer::body_add_break(doc)
    }
  }

  print(doc, target = combined_path)
  cat("\nWrote combined document:", combined_path, "\n")
} else {
  cat("\nNo eTables were available to combine -- run the upstream steps first.\n")
}

cat("\nDone. Open the .docx files in Word (or LibreOffice) to view the styled\n")
cat("tables, or open the .csv files in Excel/Numbers if you want to edit.\n")