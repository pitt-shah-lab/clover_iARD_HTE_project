# ============================================================================
# 07_cate_modeling.R
#
# PURPOSE (Victor's stage 5): the main event. Mirrors
# ADRENAL_internal_derivation_v2.R exactly, but adapted for CLOVERS:
#
# Runs the CATE comparison TWICE -- once with clinical covariates only,
# once with clinical + biomarkers -- so the project's actual question
# (do biomarkers add HTE-detection power beyond clinical variables alone?)
# gets a clean, directly comparable answer rather than a single run that
# leaves the comparison implicit.
#
# Each run:
#   1. Sources Victor's scoring.methods_v4.R, cross.validation.R, post.process_v6.R
#   2. Loads the imputed dataset and re-splits via ids_internal_list.RDS
#   3. Calls oest.cv.fun with cv_k=100 on the 5 CATE methods Victor
#      actually used for ADRENAL: cf.CATE, llrf.CATE.rlearn, bcf.CATE,
#      pbart.slearn.CATE, ModCovElasticNet.CATE
#   4. Runs each method's CV output through post.cv.fun to get AUTOC,
#      AUQINI, Kendall tau, adj AUQINI
#   5. Writes a comparison matrix CSV and a 5-method TOC plot PNG
#
# Paths come from config/config.R. Inputs from outputs/04_impute/ and
# outputs/03_split_derivation_validation/.
#
# RUNTIME WARNING: this is the slowest script in the project. cv_k=100,
# fitting 5 different CATE methods (including BART and Bayesian causal
# forest), run twice -- realistic expectation is several hours total,
# possibly overnight. Output writes only at the end of each run, so the
# terminal will look idle the whole time. It is not frozen.
#
# Output (written to outputs/07_cate_modeling/):
#   cate_resmat_clinical.csv      -- AUTOC/AUQINI/etc, clinical-only run
#   cate_resmat_with_biomarker.csv -- AUTOC/AUQINI/etc, +biomarker run
#   cate_toc_plot_clinical.png     -- 5-method TOC plot, clinical-only
#   cate_toc_plot_with_biomarker.png -- 5-method TOC plot, +biomarker
#   cate_headline_comparison.csv   -- side-by-side AUTOC comparison,
#                                      the project's actual answer
#   cate_modeling_report.txt       -- full console output captured
# ============================================================================

source(file.path(dirname(getwd()), "config", "config.R"))
out_dir <- make_output_subdir("07_cate_modeling")

report_path <- file.path(out_dir, "cate_modeling_report.txt")
report_con  <- file(report_path, open = "wt")
sink(report_con, split = TRUE)

# ---- Victor's shared toolkit ----
source(SCORING_METHODS_R)
source(CROSS_VALIDATION_R)
source(POST_PROCESS_R)

# ---- load imputed data + the split ----
dsi_all <- readRDS(file.path(OUTPUTS_DIR, "04_impute", "dsi_all.RDS"))
ids     <- readRDS(file.path(OUTPUTS_DIR, "03_split_derivation_validation",
                              "ids_internal_list.RDS"))

train <- dsi_all[dsi_all$id %in% ids$ids_der, ]
rownames(train) <- 1:nrow(train)

cat("Derivation set for CATE:", nrow(train), "patients\n\n")

# ---- column lists matching the flat file builder ----
clinical_norm <- c("age","temp","rr","hr","map","sbp",
                    "sofa","albumin","ln_bili",
                    "cr","bun","ln_g","ln_lac","sqrt_plt",
                    "ln_wbc","hgb","na","bicarb","prefluid",
                    "ln_bmi","gcs","charlson","o2sat","s2f")
clinical_bin  <- c("site_lung","site_abdom","site_urine","mv",
                    "vaso","ards","dial","chf","copd",
                    "liver","kidney")
biomarker_vars <- grep("^ln_(il1|ang1|ang2|tnfr|il6|strem1|kim1|srage)",
                        names(train), value = TRUE)

clinical_only_x  <- c(clinical_norm, clinical_bin)
with_biomarker_x <- c(clinical_norm, clinical_bin, biomarker_vars)

cat("Clinical-only covariate count: ", length(clinical_only_x), "\n")
cat("With-biomarker covariate count:", length(with_biomarker_x),
    "(", length(biomarker_vars), "biomarkers added)\n\n")

# the 5 methods Victor actually used for ADRENAL
mods <- c("cf.CATE",
          "llrf.CATE.rlearn",
          "bcf.CATE",
          "pbart.slearn.CATE",
          "ModCovElasticNet.CATE")

# ============================================================================
# Helper that does ONE full CV + post-processing pass on a given x vector.
# Returns a list with the resmat, the TOC objects (for plotting), and the
# mean.ate (for the dashed reference line in the TOC plot).
# ============================================================================

run_one_pass <- function(x_vars, label) {
  cat("============================================================\n")
  cat("CATE pass:", label, "\n")
  cat("  ", length(x_vars), "covariates\n")
  cat("  cv_k = 100, methods:", paste(mods, collapse = ", "), "\n")
  cat("============================================================\n")

  setA_list <- list(train)

  t0 <- Sys.time()
  cv <- oest.cv.fun(
    setA_list,
    cv.tr.pr = 0.5,
    cv_k = 100,
    y = "inhosp90",
    w = "w",
    x = x_vars,
    pi = 0.5,
    mods = mods,
    cores = 1,
    seed = 15223
  )
  cat("  oest.cv.fun runtime:", format(round(difftime(Sys.time(), t0), 1)), "\n")

  # ---- post.cv.fun for each method ----
  cat("  Running post.cv.fun for each method...\n")
  toc_list <- list()
  for (m in mods) {
    toc_list[[m]] <- post.cv.fun(cv[[m]],
                                  stats = c("stats.bin"),
                                  scale = 20,
                                  cores = 1,
                                  direction = "both",
                                  score = "cate")
  }

  # ---- build resmat (4 metrics x 5 methods) ----
  rates <- sapply(toc_list, function(t) colMeans(t$ratelist[[1]]))
  # rates is now 4 rows x 5 cols
  resmat <- rates
  resmat[c(1, 2, 4), ] <- round(resmat[c(1, 2, 4), ] * 100, 3)
  resmat[-4, ]         <- -resmat[-4, ]   # flip sign so Y(0)-Y(1)
  colnames(resmat) <- c("cf", "llrf-r", "bcf", "pbart", "Tian elastic net")
  rownames(resmat) <- c("autoc", "auqini", "kendall tau", "adj auqini")

  cat("  Resmat:\n")
  print(resmat)

  # ---- mean.ate for the dashed line on the TOC plot ----
  mean.ate <- mean(sapply(cv[[1]], function(x) {
    dat <- x[[1]]
    mean(dat$Y[dat$W == 1]) - mean(dat$Y[dat$W == 0])
  }))

  list(cv = cv, toc_list = toc_list, resmat = resmat, mean.ate = mean.ate)
}

# ============================================================================
# Helper to draw the 5-method TOC plot.
# ============================================================================

draw_toc_plot <- function(pass_result, label, png_path) {
  png(png_path, width = 900, height = 650)
  par(mar = c(5, 5, 4, 2))

  qq <- seq(0, 1, 0.05)[-1]
  qq <- qq[-length(qq)]
  qq <- c(qq, 0.995, 1)

  tl <- pass_result$toc_list
  ma <- pass_result$mean.ate

  get_y <- function(toc) -c(rev(toc$outmat$ad_DropHighFirst), ma) * 100

  plot(x = qq, y = get_y(tl[["cf.CATE"]]),
       type = "l", lwd = 2, col = 1, bty = "l",
       ylim = c(-12, 12),
       xlab = "Proportion treated",
       ylab = "Absolute Risk Difference (%)",
       main = paste("CATE TOC curves --", label),
       cex.axis = 1.2, cex.lab = 1.2)
  points(qq, get_y(tl[["llrf.CATE.rlearn"]]),       type = "l", col = 2, lwd = 2)
  points(qq, get_y(tl[["bcf.CATE"]]),                type = "l", col = 3, lwd = 2)
  points(qq, get_y(tl[["pbart.slearn.CATE"]]),       type = "l", col = 4, lwd = 2)
  points(qq, get_y(tl[["ModCovElasticNet.CATE"]]),   type = "l", col = 5, lwd = 2)
  abline(h = -ma * 100, lty = 2)

  legend("topright",
         col = 1:5, lty = 1, lwd = 3,
         legend = c("Causal Forest",
                    "Local Linear RF R-learner",
                    "Bayesian Causal Forest",
                    "Probit BART S-learner",
                    "Tian Elastic Net"),
         bty = "n")
  dev.off()
  cat("  Wrote TOC plot to:", png_path, "\n")
}

# ============================================================================
# Run BOTH passes and compare.
# ============================================================================

cat("\n>>> Pass 1 of 2: clinical covariates only\n\n")
pass_clinical <- run_one_pass(clinical_only_x, "clinical only")

write.csv(pass_clinical$resmat,
          file.path(out_dir, "cate_resmat_clinical.csv"))
draw_toc_plot(pass_clinical, "clinical only",
              file.path(out_dir, "cate_toc_plot_clinical.png"))

cat("\n>>> Pass 2 of 2: clinical covariates + biomarkers\n\n")
pass_with_bio <- run_one_pass(with_biomarker_x, "clinical + biomarker")

write.csv(pass_with_bio$resmat,
          file.path(out_dir, "cate_resmat_with_biomarker.csv"))
draw_toc_plot(pass_with_bio, "clinical + biomarker",
              file.path(out_dir, "cate_toc_plot_with_biomarker.png"))

# ============================================================================
# Headline comparison: same metrics, two columns side by side.
# This is the actual answer to the project's question.
# ============================================================================

cat("\n============================================================\n")
cat("HEADLINE COMPARISON: clinical vs clinical+biomarker\n")
cat("(AUTOC: higher = more detected heterogeneity; AUQINI similar)\n")
cat("============================================================\n")

headline <- data.frame(
  metric           = rownames(pass_clinical$resmat),
  cf_clinical      = pass_clinical$resmat[, "cf"],
  cf_with_bio      = pass_with_bio$resmat[, "cf"],
  llrf_clinical    = pass_clinical$resmat[, "llrf-r"],
  llrf_with_bio    = pass_with_bio$resmat[, "llrf-r"],
  bcf_clinical     = pass_clinical$resmat[, "bcf"],
  bcf_with_bio     = pass_with_bio$resmat[, "bcf"],
  pbart_clinical   = pass_clinical$resmat[, "pbart"],
  pbart_with_bio   = pass_with_bio$resmat[, "pbart"],
  enet_clinical    = pass_clinical$resmat[, "Tian elastic net"],
  enet_with_bio    = pass_with_bio$resmat[, "Tian elastic net"],
  stringsAsFactors = FALSE
)
print(headline, row.names = FALSE)

write.csv(headline,
          file.path(out_dir, "cate_headline_comparison.csv"),
          row.names = FALSE)

cat("\nAll outputs written to:\n  ", out_dir, "\n")

sink()
close(report_con)
cat("\nFull report saved to:\n", report_path, "\n")