# ============================================================================
# 10b_confirmatory_analysis_with_yhat.R
#
# Same as 10_confirmatory_analysis.R but with two changes:
#   1. Risk-model Y.hat is passed to causal_forest() (Victor's request)
#   2. Outputs go to outputs/10b_confirmatory_with_yhat/ (preserves original)
#
# The original results (without Y.hat) remain untouched in
# outputs/10_confirmatory_analysis/.
#
# After this runs, re-run 08b_gates_with_yhat.R on the new predictions.
# ============================================================================

this_script <- sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])
source(file.path(dirname(dirname(this_script)), "config", "config.R"))
out_dir <- make_output_subdir("10b_confirmatory_with_yhat")

report_path <- file.path(out_dir, "confirmatory_report.txt")
report_con  <- file(report_path, open = "wt")
sink(report_con, split = TRUE)

# ---- source tools ----
source(SCORING_METHODS_R)
source(CROSS_VALIDATION_R)
source(POST_PROCESS_R)
source(CONFIRM_ANALYSIS_R)

library(glmnet)
library(grf)

# ---- load data ----
dsi_all <- readRDS(file.path(OUTPUTS_DIR, "04_impute", "dsi_all.RDS"))
ids     <- readRDS(file.path(OUTPUTS_DIR, "03_split_derivation_validation",
                              "ids_internal_list.RDS"))

der <- dsi_all[dsi_all$id %in% ids$ids_der, ]
val <- dsi_all[dsi_all$id %in% ids$ids_val, ]
rownames(der) <- 1:nrow(der)
rownames(val) <- 1:nrow(val)

cat("Derivation set:", nrow(der), "patients\n")
cat("Validation set:", nrow(val), "patients\n\n")

# ---- covariate lists ----
clinical_norm <- c("age","temp","rr","hr","map","sbp",
                    "sofa","albumin","ln_bili",
                    "cr","bun","ln_g","ln_lac","sqrt_plt",
                    "ln_wbc","hgb","na","bicarb","prefluid",
                    "ln_bmi","gcs","charlson","o2sat","s2f")
clinical_bin  <- c("site_lung","site_abdom","site_urine","mv",
                    "vaso","ards","dial","chf","copd",
                    "liver","kidney")
biomarker_vars <- grep("^ln_(il1|ang1|ang2|tnfr|il6|strem1|kim1|srage)",
                        names(der), value = TRUE)

clinical_only_x  <- c(clinical_norm, clinical_bin)
with_biomarker_x <- c(clinical_norm, clinical_bin, biomarker_vars)

# ============================================================================
# STEP 0: Fit risk model and generate Y.hat for BOTH derivation and validation
# ============================================================================

cat("============================================================\n")
cat("STEP 0: Fit risk model for Y.hat\n")
cat("============================================================\n")

non_x  <- c("id", "inhosp90", "sofa_diff", "w")
all_x  <- setdiff(names(der), non_x)

ctrl   <- der[der$w == 0, ]
X_ctrl <- as.matrix(ctrl[, all_x])
y_ctrl <- ctrl$inhosp90

set.seed(15223)
cv_fit <- cv.glmnet(X_ctrl, y_ctrl, family = "binomial", alpha = 0.5,
                    nfolds = 10, type.measure = "auc")

Y_hat_der <- as.numeric(predict(cv_fit, newx = as.matrix(der[, all_x]),
                                 s = "lambda.min", type = "response"))
Y_hat_val <- as.numeric(predict(cv_fit, newx = as.matrix(val[, all_x]),
                                 s = "lambda.min", type = "response"))

cat("  Y.hat derivation: mean =", round(mean(Y_hat_der), 4),
    ", range", round(min(Y_hat_der), 4), "-", round(max(Y_hat_der), 4), "\n")
cat("  Y.hat validation: mean =", round(mean(Y_hat_val), 4),
    ", range", round(min(Y_hat_val), 4), "-", round(max(Y_hat_val), 4), "\n\n")

# ============================================================================
# STEP 1: Train cf.CATE on derivation WITH Y.hat, predict on validation
#
# We bypass confirm_analysis() and call causal_forest directly so we can
# pass Y.hat. This is cleaner than overriding the function.
# ============================================================================

cat("============================================================\n")
cat("STEP 1: CATE predictions on validation set (with Y.hat)\n")
cat("============================================================\n\n")

run_cf_with_yhat <- function(x_vars, label) {
  cat("Running cf.CATE with Y.hat:", label, "\n")
  t0 <- Sys.time()

  X_der <- as.matrix(der[, x_vars])
  X_val <- as.matrix(val[, x_vars])

  cf <- causal_forest(
    X       = X_der,
    Y       = der$inhosp90,
    W       = der$w,
    Y.hat   = Y_hat_der,
    W.hat   = rep(0.5, nrow(der)),
    honesty = TRUE,
    honesty.fraction = 0.5,
    num.trees = 5000,
    sample.fraction = 0.5,
    min.node.size = 5,
    alpha = 0.05,
    imbalance.penalty = 0,
    mtry = ceiling(2 * length(x_vars) / 3)
  )

  # Predict on validation
  cf_pred <- predict(cf, newdata = X_val)$predictions

  # Variable importance
  cf_imp <- variable_importance(cf)

  # Build the (y, w, score, b, s) data frame that gates.fun expects
  preds_df <- data.frame(
    y     = val$inhosp90,
    w     = val$w,
    score = cf_pred,
    b     = Y_hat_val,
    s     = cf_pred
  )

  elapsed <- round(difftime(Sys.time(), t0), 1)
  cat("  Runtime:", format(elapsed), "\n")
  cat("  Score range:", round(min(cf_pred), 4), "to", round(max(cf_pred), 4),
      "  median:", round(median(cf_pred), 4), "\n\n")

  list(preds_df = preds_df, importance = cf_imp, x_vars = x_vars)
}

res_clin <- run_cf_with_yhat(clinical_only_x,  "clinical only (35 covariates)")
res_bio  <- run_cf_with_yhat(with_biomarker_x, "clinical + biomarker (43 covariates)")

# ---- Save predictions ----
write.csv(res_clin$preds_df,
          file.path(out_dir, "confirmatory_preds_clinical.csv"),
          row.names = FALSE)
write.csv(res_bio$preds_df,
          file.path(out_dir, "confirmatory_preds_with_biomarker.csv"),
          row.names = FALSE)

# ---- Save importance ----
for (res in list(list(r = res_clin, name = "clinical"),
                 list(r = res_bio,  name = "with_biomarker"))) {
  imp_df <- data.frame(
    variable   = res$r$x_vars,
    importance = as.numeric(res$r$importance),
    stringsAsFactors = FALSE
  )
  imp_df <- imp_df[order(-imp_df$importance), ]
  write.csv(imp_df,
            file.path(out_dir, paste0("confirmatory_importance_", res$name, ".csv")),
            row.names = FALSE)
  cat("Top 10 variables (", res$name, "):\n")
  print(head(imp_df, 10), row.names = FALSE)
  cat("\n")
}

# ============================================================================
# STEP 2: Subgroup classification (same tertile approach as original)
# ============================================================================

cat("============================================================\n")
cat("STEP 2: Subgroup classification and treatment effect testing\n")
cat("============================================================\n\n")

classify_and_test <- function(preds_df, label) {
  cat("--- ", label, " ---\n")

  score <- preds_df$score
  y_val <- preds_df$y
  w_val <- preds_df$w

  q33 <- quantile(score, 1/3)
  q67 <- quantile(score, 2/3)

  cat("  Score distribution: median =", round(median(score), 4),
      ", Q33 =", round(q33, 4), ", Q67 =", round(q67, 4), "\n")

  sub_ben   <- as.integer(score < q33)
  sub_harm  <- as.integer(score > q67)
  sub_indet <- as.integer(score >= q33 & score <= q67)

  cat("  Benefit group:       n =", sum(sub_ben), "\n")
  cat("  Indeterminate group: n =", sum(sub_indet), "\n")
  cat("  Harm group:          n =", sum(sub_harm), "\n\n")

  test_subgroup <- function(sub_indicator, alt) {
    n_w1 <- sum(sub_indicator * w_val)
    n_w0 <- sum(sub_indicator * (1 - w_val))
    d_w1 <- sum(sub_indicator * w_val * y_val)
    d_w0 <- sum(sub_indicator * (1 - w_val) * y_val)

    if (n_w1 == 0 || n_w0 == 0) {
      return(c(N = sum(sub_indicator), n_w1 = n_w1, n_w0 = n_w0,
               deaths_w1 = d_w1, rate_w1 = NA,
               deaths_w0 = d_w0, rate_w0 = NA,
               rate_diff = NA, pval = NA))
    }

    rate_w1 <- d_w1 / n_w1
    rate_w0 <- d_w0 / n_w0
    rate_diff <- rate_w1 - rate_w0

    pv <- tryCatch(
      prop.test(x = c(d_w1, d_w0), n = c(n_w1, n_w0),
                alternative = alt)$p.value,
      error = function(e) NA
    )

    c(N = sum(sub_indicator), n_w1 = n_w1, n_w0 = n_w0,
      deaths_w1 = d_w1, rate_w1 = round(rate_w1, 4),
      deaths_w0 = d_w0, rate_w0 = round(rate_w0, 4),
      rate_diff = round(rate_diff, 4), pval = round(pv, 4))
  }

  resmat <- rbind(
    benefit       = test_subgroup(sub_ben,   "less"),
    indeterminate = test_subgroup(sub_indet, "two.sided"),
    harm          = test_subgroup(sub_harm,  "greater")
  )

  cat("  Subgroup treatment effects (inhosp90):\n")
  print(as.data.frame(resmat))
  cat("\n")

  as.data.frame(resmat)
}

resmat_clinical <- classify_and_test(res_clin$preds_df, "Clinical only + Y.hat")
resmat_bio      <- classify_and_test(res_bio$preds_df,  "Clinical + biomarker + Y.hat")

resmat_clinical$model    <- "clinical_only_with_yhat"
resmat_bio$model         <- "clinical_plus_biomarker_with_yhat"
resmat_clinical$subgroup <- rownames(resmat_clinical)
resmat_bio$subgroup      <- rownames(resmat_bio)

combined <- rbind(resmat_clinical, resmat_bio)
write.csv(combined, file.path(out_dir, "confirmatory_subgroups.csv"),
          row.names = FALSE)



cat("\nDone. Report saved to:", report_path, "\n")