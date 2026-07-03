# ============================================================================
# 05_risk_modeling.R
#
# PURPOSE (Victor's stage 4): run several baseline-risk prediction models
# on the imputed dataset, using HIS own cross-validation framework
# (oest.cv.fun from cate-repo/cross.validation.R), and compare them on
# out-of-sample discrimination (AUC for inhosp90).
#
# Methods compared (six total):
#   Victor's:  logreg.risk, rf.risk, rf.risk.10
#   Added:     enet.risk, xgb.risk, bart.risk
# All have identical calling signatures and plug into oest.cv.fun by name --
# this driver doesn't special-case any of them.
#
# Paths come from config/config.R. Inputs from outputs/04_impute/.
#
# Output (written to outputs/05_risk_modeling/):
#   risk_model_comparison.csv  -- AUC + summary per method
#   risk_model_auc_boxplot.png -- per-fold AUC distribution per method
#   risk_modeling_report.txt   -- full console output captured
# ============================================================================

source(file.path(dirname(getwd()), "config", "config.R"))
out_dir <- make_output_subdir("05_risk_modeling")

report_path <- file.path(out_dir, "risk_modeling_report.txt")
report_con  <- file(report_path, open = "wt")
sink(report_con, split = TRUE)

# ---- source Victor's tools + our additions ---------------------------------
source(SCORING_METHODS_R)                                       # logreg.risk, rf.risk, rf.risk.10
source(CROSS_VALIDATION_R)                                      # oest.cv.fun
source(file.path(PROJECT_ROOT, "preprocessing", "05_additional_risk_models.R"))  # enet.risk, xgb.risk, bart.risk

# ---- load imputed data -----------------------------------------------------
dsi_all <- readRDS(file.path(OUTPUTS_DIR, "04_impute", "dsi_all.RDS"))
cat("Loaded imputed dataset:", nrow(dsi_all), "patients,",
    ncol(dsi_all), "columns\n\n")

# oest.cv.fun expects train to be a LIST of imputed datasets (length M).
# We did single imputation, so M=1 and train is a length-1 list.
train_list <- list(dsi_all)
M <- 1

# Identify column roles
y       <- "inhosp90"
w       <- "w"
non_x   <- c("id", "inhosp90", "sofa_diff", "w")
x       <- setdiff(names(dsi_all), non_x)

cat("Outcome      :", y, "\n")
cat("Treatment    :", w, "\n")
cat("# covariates :", length(x), "\n\n")

# ---- choose the methods to run --------------------------------------------
mods <- c("logreg.risk", "rf.risk", "rf.risk.10",
          "enet.risk", "bart.risk")

# mods <- c("logreg.risk", "rf.risk", "rf.risk.10",
#           "enet.risk",   "xgb.risk", "bart.risk")

# ---- CV settings -----------------------------------------------------------
# cv_k = 50 repeated 50/50 splits. Each method gets fit on half the data and
# evaluated on the other half, 50 times, then we average. Victor used cv_k=100
# for the ADRENAL CATE comparison; risk modeling per iteration is lighter,
# but 50 is plenty stable for an AUC estimate and noticeably faster.
cv_k <- 50

cat("Running cv_k =", cv_k, "repeated 50/50 splits, methods:\n")
cat(" ", paste(mods, collapse = ", "), "\n\n")

# ---- run the comparison ---------------------------------------------------
t0 <- Sys.time()
res <- oest.cv.fun(
  train = train_list,
  cv.tr.pr = 0.5,
  cv_k = cv_k,
  y = y, w = w, x = x,
  pi = 0.5,        # 1:1 randomization in CLOVERS
  mods = mods,
  cores = 1,
  seed = 15223
)
cat("Total runtime:", format(round(difftime(Sys.time(), t0), 1)), "\n\n")

# ---- compute per-fold AUC for each method ---------------------------------
# res is structured as: res[[mod_name]][[k]] = list(tempeval, importance)
# where tempeval has columns Y, W, score. We compute AUC of `score` against
# Y on the eval set for each (mod, k) pair.

simple_auc <- function(pred, truth) {
  # AUC = probability that a random positive has a higher score than a
  # random negative. Equivalent to the Mann-Whitney U normalization.
  pos <- pred[truth == 1]; neg <- pred[truth == 0]
  if (length(pos) == 0 || length(neg) == 0) return(NA)
  n_concordant <- sum(outer(pos, neg, ">")) + 0.5 * sum(outer(pos, neg, "=="))
  n_concordant / (length(pos) * length(neg))
}

auc_mat <- matrix(NA, nrow = cv_k, ncol = length(mods))
colnames(auc_mat) <- mods

for (j in seq_along(mods)) {
  for (k in seq_len(cv_k)) {
    fold      <- res[[mods[j]]][[k]][[1]]
    auc_mat[k, j] <- simple_auc(fold$score, fold$Y)
  }
}

# ---- summarize ------------------------------------------------------------
auc_summary <- data.frame(
  method     = mods,
  mean_auc   = round(colMeans(auc_mat, na.rm = TRUE), 4),
  sd_auc     = round(apply(auc_mat, 2, sd, na.rm = TRUE), 4),
  median_auc = round(apply(auc_mat, 2, median, na.rm = TRUE), 4),
  q25_auc    = round(apply(auc_mat, 2, quantile, 0.25, na.rm = TRUE), 4),
  q75_auc    = round(apply(auc_mat, 2, quantile, 0.75, na.rm = TRUE), 4),
  stringsAsFactors = FALSE
)
auc_summary <- auc_summary[order(-auc_summary$mean_auc), ]

cat("AUC summary across", cv_k, "cross-validation folds, sorted best first:\n\n")
print(auc_summary, row.names = FALSE)

write.csv(auc_summary, file.path(out_dir, "risk_model_comparison.csv"),
          row.names = FALSE)
write.csv(auc_mat, file.path(out_dir, "risk_model_auc_per_fold.csv"),
          row.names = FALSE)

# ---- per-fold boxplot ------------------------------------------------------
plot_path <- file.path(out_dir, "risk_model_auc_boxplot.png")
png(plot_path, width = 900, height = 550)
par(mar = c(8, 4, 4, 2), las = 2)
boxplot(auc_mat,
        main = sprintf("Out-of-sample AUC across %d CV folds, by method", cv_k),
        ylab = "AUC", col = "#cfe2ff",
        names = mods)
dev.off()

cat("\nWrote:\n")
cat("  ", file.path(out_dir, "risk_model_comparison.csv"), "\n")
cat("  ", file.path(out_dir, "risk_model_auc_per_fold.csv"), "\n")
cat("  ", plot_path, "\n")

sink()
close(report_con)
cat("\nFull report saved to:\n  ", report_path, "\n")