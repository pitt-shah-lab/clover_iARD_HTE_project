# ============================================================================
# 07b_cate_modeling_with_risk.R
#
# This is the updated CATE modeling script that addresses Victor's two requests:
#
#   1. Feed risk predictions into CATE models (Y.hat)
#      - The selected risk model (enet) is fit on the full derivation
#        control arm and used to predict Y.hat for all derivation patients.
#      - Y.hat is passed to causal_forest() instead of NULL, so the forest
#        separates prognostic from predictive (HTE) effects.
#
#   2. Third covariate-set pass: risk score + biomarkers only
#      - Pass 1: clinical covariates only (35 vars)
#      - Pass 2: clinical + biomarkers (43 vars)
#      - Pass 3: risk score + biomarkers only (9 vars)
#      - Pass 3 answers: "do biomarkers predict who benefits even when all
#        you know about the patient's clinical state is a single summary?"
#
# This script replaces 07_cate_modeling.R. It sources the same shared
# toolkit but overrides cf.CATE to accept Y.hat.
#
# RUNTIME: several hours (3 passes × 100 folds × 5 methods). Run overnight.
#
# Output (written to outputs/07b_cate_modeling_with_risk/):
#   risk_predictions_derivation.csv
#   cate_resmat_clinical.csv
#   cate_resmat_with_biomarker.csv
#   cate_resmat_risk_plus_biomarker.csv
#   cate_toc_plot_*.png  (one per pass)
#   cate_headline_comparison.csv
#   cate_modeling_report.txt
# ============================================================================

this_script <- sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])
source(file.path(dirname(dirname(this_script)), "config", "config.R"))
out_dir <- make_output_subdir("07b_cate_modeling_with_risk")

report_path <- file.path(out_dir, "cate_modeling_report.txt")
report_con  <- file(report_path, open = "wt")
sink(report_con, split = TRUE)

# ---- Victor's shared toolkit ----
source(SCORING_METHODS_R)
source(CROSS_VALIDATION_R)
source(POST_PROCESS_R)
source(file.path(PROJECT_ROOT, "scripts", "05_additional_risk_models.R"))

# ---- load imputed data + the split ----
dsi_all <- readRDS(file.path(OUTPUTS_DIR, "04_impute", "dsi_all.RDS"))
ids     <- readRDS(file.path(OUTPUTS_DIR, "03_split_derivation_validation",
                              "ids_internal_list.RDS"))

train <- dsi_all[dsi_all$id %in% ids$ids_der, ]
rownames(train) <- 1:nrow(train)

cat("Derivation set for CATE:", nrow(train), "patients\n\n")

# ============================================================================
# STEP 1: Generate risk predictions (Y.hat) on the full derivation set.
#
# Fit the best risk model (elastic net) on the control arm of the derivation
# set, then predict for ALL derivation patients. This gives each patient an
# estimated P(Y=1 | X, W=0) — their baseline (untreated) risk.
# ============================================================================

cat("============================================================\n")
cat("STEP 1: Fit risk model and generate Y.hat\n")
cat("============================================================\n")

library(glmnet)

non_x <- c("id", "inhosp90", "sofa_diff", "w")
all_x <- setdiff(names(train), non_x)

# Fit elastic net on control arm only
ctrl <- train[train$w == 0, ]
X_ctrl <- as.matrix(ctrl[, all_x])
y_ctrl <- ctrl$inhosp90

set.seed(15223)
cv_fit <- cv.glmnet(X_ctrl, y_ctrl, family = "binomial", alpha = 0.5,
                    nfolds = 10, type.measure = "auc")

# Predict for ALL derivation patients (both arms)
X_all  <- as.matrix(train[, all_x])
Y_hat  <- as.numeric(predict(cv_fit, newx = X_all, s = "lambda.min",
                              type = "response"))

cat("  Risk model: elastic net (alpha=0.5), fit on", nrow(ctrl),
    "control-arm patients\n")
cat("  Y.hat range:", round(min(Y_hat), 4), "to", round(max(Y_hat), 4),
    "  mean:", round(mean(Y_hat), 4), "\n")

# Add risk score as a column for Pass 3
train$risk_score <- Y_hat

# Save for reproducibility
risk_df <- data.frame(id = train$id, risk_score = Y_hat)
write.csv(risk_df, file.path(out_dir, "risk_predictions_derivation.csv"),
          row.names = FALSE)
cat("  Saved risk predictions to: risk_predictions_derivation.csv\n\n")

# ============================================================================
# STEP 2: Override cf.CATE to accept Y.hat
#
# Victor's cf.CATE hardcodes Y.hat=NULL inside causal_forest(). We redefine
# it here to accept an optional Y.hat vector. When provided, each training
# fold uses the corresponding subset of Y.hat. When NULL, it falls back to
# the original behavior.
# ============================================================================

# Store the Y.hat vector in an environment the function can see
.yhat_env <- new.env(parent = emptyenv())
.yhat_env$Y_hat_full <- Y_hat

cf.CATE <- function(train, cvt, cve, M, y, w, x, pi = 0.5,
                    ntrees = 5000, tune = FALSE, predY = FALSE, ...) {

  cf.pred.all <- cf.imp.all <- cf.predy.all <- c()

  for (m in 1:M) {
    ds  <- train[[m]]
    dst <- ds[cvt, ]
    dse <- ds[cve, ]

    if (tune) {
      hyp <- tune_causal_forest(
        X = as.matrix(dst[, x]), Y = dst[, y], W = dst[, w],
        Y.hat = mean(dst[, y]), W.hat = pi,
        tune.num.trees = 1000,
        tune.parameters = c("sample.fraction", "mtry", "min.node.size",
                            "alpha", "imbalance.penalty"))
    } else {
      hyp <- list("params" = c(5, 0.5, ceiling(2 * length(x) / 3), 0.05, 0))
      names(hyp[[1]]) <- c("min.node.size", "sample.fraction",
                            "mtry", "alpha", "imbalance.penalty")
    }

    # --- THE KEY CHANGE: use Y.hat from the risk model ---
    if (!is.null(.yhat_env$Y_hat_full)) {
      yhat_train <- .yhat_env$Y_hat_full[cvt]
    } else {
      yhat_train <- NULL
    }

    cf <- causal_forest(
      X = as.matrix(dst[, x]),
      Y = dst[, y],
      W = dst[, w],
      Y.hat = yhat_train,
      W.hat = pi,
      honesty = TRUE,
      honesty.fraction = 0.5,
      num.trees = ntrees,
      sample.fraction = as.numeric(hyp$params["sample.fraction"]),
      min.node.size   = as.numeric(hyp$params["min.node.size"]),
      alpha           = as.numeric(hyp$params["alpha"]),
      imbalance.penalty = as.numeric(hyp$params["imbalance.penalty"]),
      mtry            = as.numeric(hyp$params["mtry"]))

    if (identical(dst, dse)) {
      cf.pred <- predict(cf)$predictions
    } else {
      cf.pred <- predict(cf, newdata = as.matrix(dse[, x]))$predictions
    }

    cf.pred.all <- cbind(cf.pred.all, cf.pred)
    cf.imp      <- variable_importance(cf)
    cf.imp.all  <- cbind(cf.imp.all, cf.imp)

    if (predY) {
      a_fromCF <- get_forest_weights(cf, newdata = as.matrix(dse[, x]))
      sum_aY_mu0 <- a_fromCF[dse$w == 0, dst$w == 0] %*%
                    matrix(dst[dst$w == 0, y], , 1)
      sum_a_mu0  <- apply(a_fromCF[dse$w == 0, dst$w == 0], 1, sum)
      preds_mu0  <- as.numeric(sum_aY_mu0 / sum_a_mu0)
      sum_aY_mu1 <- a_fromCF[dse$w == 1, dst$w == 1] %*%
                    matrix(dst[dst$w == 1, y], , 1)
      sum_a_mu1  <- apply(a_fromCF[dse$w == 1, dst$w == 1], 1, sum)
      preds_mu1  <- as.numeric(sum_aY_mu1 / sum_a_mu1)
      pred.y <- c(preds_mu0, preds_mu1)[order(order(dse$w))]
      cf.predy.all <- cbind(cf.predy.all, pred.y)
    } else {
      cf.predy.all <- cbind(cf.predy.all, rep(0, length(cve)))
    }
  }

  list(apply(cf.pred.all, 1, mean),
       apply(cf.imp.all, 1, mean),
       apply(cf.predy.all, 1, mean))
}

cat("cf.CATE overridden to use risk-model Y.hat\n\n")

# ---- column lists ----
clinical_norm <- c("age", "temp", "rr", "hr", "map", "sbp",
                    "sofa", "albumin", "ln_bili",
                    "cr", "bun", "ln_g", "ln_lac", "sqrt_plt",
                    "ln_wbc", "hgb", "na", "bicarb", "prefluid",
                    "ln_bmi", "gcs", "charlson", "o2sat", "s2f")
clinical_bin  <- c("site_lung", "site_abdom", "site_urine", "mv",
                    "vaso", "ards", "dial", "chf", "copd",
                    "liver", "kidney")
biomarker_vars <- grep("^ln_(il1|ang1|ang2|tnfr|il6|strem1|kim1|srage)",
                        names(train), value = TRUE)

clinical_only_x      <- c(clinical_norm, clinical_bin)
with_biomarker_x     <- c(clinical_norm, clinical_bin, biomarker_vars)
risk_plus_biomarker_x <- c("risk_score", biomarker_vars)

cat("Pass 1 — clinical only:        ", length(clinical_only_x), "covariates\n")
cat("Pass 2 — clinical + biomarker: ", length(with_biomarker_x), "covariates\n")
cat("Pass 3 — risk + biomarker:     ", length(risk_plus_biomarker_x), "covariates\n\n")

mods <- c("cf.CATE",
          "llrf.CATE.rlearn",
          "bcf.CATE",
          "pbart.slearn.CATE",
          "ModCovElasticNet.CATE")

# ============================================================================
# run_one_pass — identical to 07_cate_modeling.R
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

  rates <- sapply(toc_list, function(t) colMeans(t$ratelist[[1]]))
  resmat <- rates
  resmat[c(1, 2, 4), ] <- round(resmat[c(1, 2, 4), ] * 100, 3)
  resmat[-4, ]         <- -resmat[-4, ]
  colnames(resmat) <- c("cf", "llrf-r", "bcf", "pbart", "Tian elastic net")
  rownames(resmat) <- c("autoc", "auqini", "kendall tau", "adj auqini")

  cat("  Resmat:\n")
  print(resmat)

  mean.ate <- mean(sapply(cv[[1]], function(x) {
    dat <- x[[1]]
    mean(dat$Y[dat$W == 1]) - mean(dat$Y[dat$W == 0])
  }))

  list(cv = cv, toc_list = toc_list, resmat = resmat, mean.ate = mean.ate)
}

# ============================================================================
# draw_toc_plot — identical to 07_cate_modeling.R
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
# Run ALL THREE passes
# ============================================================================

cat("\n>>> Pass 1 of 3: clinical covariates only (with Y.hat from risk model)\n\n")
pass_clinical <- run_one_pass(clinical_only_x, "clinical only + Y.hat")
write.csv(pass_clinical$resmat,
          file.path(out_dir, "cate_resmat_clinical.csv"))
draw_toc_plot(pass_clinical, "clinical only + Y.hat",
              file.path(out_dir, "cate_toc_plot_clinical.png"))

cat("\n>>> Pass 2 of 3: clinical + biomarkers (with Y.hat from risk model)\n\n")
pass_with_bio <- run_one_pass(with_biomarker_x, "clinical + biomarker + Y.hat")
write.csv(pass_with_bio$resmat,
          file.path(out_dir, "cate_resmat_with_biomarker.csv"))
draw_toc_plot(pass_with_bio, "clinical + biomarker + Y.hat",
              file.path(out_dir, "cate_toc_plot_with_biomarker.png"))

cat("\n>>> Pass 3 of 3: risk score + biomarkers only (Victor's new request)\n\n")
pass_risk_bio <- run_one_pass(risk_plus_biomarker_x, "risk + biomarker only")
write.csv(pass_risk_bio$resmat,
          file.path(out_dir, "cate_resmat_risk_plus_biomarker.csv"))
draw_toc_plot(pass_risk_bio, "risk + biomarker only",
              file.path(out_dir, "cate_toc_plot_risk_plus_biomarker.png"))

# ============================================================================
# Headline comparison: 3-way side by side
# ============================================================================

cat("\n============================================================\n")
cat("HEADLINE COMPARISON: 3-way\n")
cat("(AUTOC: higher = more detected heterogeneity)\n")
cat("============================================================\n")

headline <- data.frame(
  metric              = rownames(pass_clinical$resmat),
  cf_clin             = pass_clinical$resmat[, "cf"],
  cf_bio              = pass_with_bio$resmat[, "cf"],
  cf_risk_bio         = pass_risk_bio$resmat[, "cf"],
  llrf_clin           = pass_clinical$resmat[, "llrf-r"],
  llrf_bio            = pass_with_bio$resmat[, "llrf-r"],
  llrf_risk_bio       = pass_risk_bio$resmat[, "llrf-r"],
  bcf_clin            = pass_clinical$resmat[, "bcf"],
  bcf_bio             = pass_with_bio$resmat[, "bcf"],
  bcf_risk_bio        = pass_risk_bio$resmat[, "bcf"],
  pbart_clin          = pass_clinical$resmat[, "pbart"],
  pbart_bio           = pass_with_bio$resmat[, "pbart"],
  pbart_risk_bio      = pass_risk_bio$resmat[, "pbart"],
  enet_clin           = pass_clinical$resmat[, "Tian elastic net"],
  enet_bio            = pass_with_bio$resmat[, "Tian elastic net"],
  enet_risk_bio       = pass_risk_bio$resmat[, "Tian elastic net"],
  stringsAsFactors    = FALSE
)
print(headline, row.names = FALSE)

write.csv(headline,
          file.path(out_dir, "cate_headline_comparison.csv"),
          row.names = FALSE)

cat("\nAll outputs written to:\n  ", out_dir, "\n")

sink()
close(report_con)
cat("\nFull report saved to:\n", report_path, "\n")