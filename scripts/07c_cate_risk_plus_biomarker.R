# ============================================================================
# 07c_cate_risk_plus_biomarker.R
#
# Victor's request: "include a branch of results from model including only
# predicted risk + biomarkers, but no other clinical markers"
#
# This runs a SINGLE CATE pass using only 9 covariates:
#   1 predicted risk score (from elastic net risk model)
#   8 log-transformed biomarkers
#
# This answers: do the biomarkers predict who benefits even when all you
# know about the patient's clinical state is a single summary risk number?
#
# Also implements Victor's second request: risk predictions (Y.hat) are
# passed to causal_forest() instead of NULL.
#
# Can be run independently — does not re-run Passes 1 or 2.
#
# Output (written to outputs/07c_cate_risk_plus_biomarker/):
#   risk_predictions_derivation.csv
#   cate_resmat_risk_plus_biomarker.csv
#   cate_toc_plot_risk_plus_biomarker.png
#   cate_report.txt
# ============================================================================

this_script <- sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])
source(file.path(dirname(dirname(this_script)), "config", "config.R"))
out_dir <- make_output_subdir("07c_cate_risk_plus_biomarker")

report_path <- file.path(out_dir, "cate_report.txt")
report_con  <- file(report_path, open = "wt")
sink(report_con, split = TRUE)

# ---- Victor's shared toolkit ----
source(SCORING_METHODS_R)
source(CROSS_VALIDATION_R)
source(POST_PROCESS_R)

library(glmnet)

# ---- load imputed data + the split ----
dsi_all <- readRDS(file.path(OUTPUTS_DIR, "04_impute", "dsi_all.RDS"))
ids     <- readRDS(file.path(OUTPUTS_DIR, "03_split_derivation_validation",
                              "ids_internal_list.RDS"))

train <- dsi_all[dsi_all$id %in% ids$ids_der, ]
rownames(train) <- 1:nrow(train)

cat("Derivation set:", nrow(train), "patients\n\n")

# ============================================================================
# STEP 1: Generate risk predictions (Y.hat) on the derivation set
# ============================================================================

cat("============================================================\n")
cat("STEP 1: Fit risk model and generate Y.hat + risk_score\n")
cat("============================================================\n")

non_x <- c("id", "inhosp90", "sofa_diff", "w")
all_x <- setdiff(names(train), non_x)

ctrl   <- train[train$w == 0, ]
X_ctrl <- as.matrix(ctrl[, all_x])
y_ctrl <- ctrl$inhosp90

set.seed(15223)
cv_fit <- cv.glmnet(X_ctrl, y_ctrl, family = "binomial", alpha = 0.5,
                    nfolds = 10, type.measure = "auc")

X_all  <- as.matrix(train[, all_x])
Y_hat  <- as.numeric(predict(cv_fit, newx = X_all, s = "lambda.min",
                              type = "response"))

train$risk_score <- Y_hat

cat("  Risk model: elastic net (alpha=0.5), fit on", nrow(ctrl),
    "control-arm patients\n")
cat("  Y.hat range:", round(min(Y_hat), 4), "to", round(max(Y_hat), 4),
    "  mean:", round(mean(Y_hat), 4), "\n")

write.csv(data.frame(id = train$id, risk_score = Y_hat),
          file.path(out_dir, "risk_predictions_derivation.csv"),
          row.names = FALSE)
cat("  Saved risk predictions\n\n")

# ============================================================================
# STEP 2: Override cf.CATE to use Y.hat from the risk model
# ============================================================================

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
        Y.hat = mean(dst[, y]), W.hat = pi, tune.num.trees = 1000,
        tune.parameters = c("sample.fraction", "mtry", "min.node.size",
                            "alpha", "imbalance.penalty"))
    } else {
      hyp <- list("params" = c(5, 0.5, ceiling(2 * length(x) / 3), 0.05, 0))
      names(hyp[[1]]) <- c("min.node.size", "sample.fraction",
                            "mtry", "alpha", "imbalance.penalty")
    }

    # Use risk-model Y.hat instead of NULL
    yhat_train <- .yhat_env$Y_hat_full[cvt]

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

# ============================================================================
# STEP 3: Run Pass 3 — risk score + biomarkers only
# ============================================================================

biomarker_vars <- grep("^ln_(il1|ang1|ang2|tnfr|il6|strem1|kim1|srage)",
                        names(train), value = TRUE)
risk_plus_biomarker_x <- c("risk_score", biomarker_vars)

cat("============================================================\n")
cat("Pass 3: risk score + biomarkers only\n")
cat("  Covariates (", length(risk_plus_biomarker_x), "):\n")
for (v in risk_plus_biomarker_x) cat("    ", v, "\n")
cat("============================================================\n\n")

mods <- c("cf.CATE",
          "llrf.CATE.rlearn",
          "bcf.CATE",
          "pbart.slearn.CATE",
          "ModCovElasticNet.CATE")

setA_list <- list(train)

t0 <- Sys.time()
cv <- oest.cv.fun(
  setA_list,
  cv.tr.pr = 0.5,
  cv_k = 100,
  y = "inhosp90",
  w = "w",
  x = risk_plus_biomarker_x,
  pi = 0.5,
  mods = mods,
  cores = 1,
  seed = 15223
)
cat("  Runtime:", format(round(difftime(Sys.time(), t0), 1)), "\n\n")

# ---- post-processing ----
cat("Running post.cv.fun for each method...\n")
toc_list <- list()
for (m in mods) {
  toc_list[[m]] <- post.cv.fun(cv[[m]],
                                stats = c("stats.bin"),
                                scale = 20,
                                cores = 1,
                                direction = "both",
                                score = "cate")
}

rates  <- sapply(toc_list, function(t) colMeans(t$ratelist[[1]]))
resmat <- rates
resmat[c(1, 2, 4), ] <- round(resmat[c(1, 2, 4), ] * 100, 3)
resmat[-4, ]         <- -resmat[-4, ]
colnames(resmat) <- c("cf", "llrf-r", "bcf", "pbart", "Tian elastic net")
rownames(resmat) <- c("autoc", "auqini", "kendall tau", "adj auqini")

cat("\nResults — risk score + biomarkers only:\n")
print(resmat)

write.csv(resmat, file.path(out_dir, "cate_resmat_risk_plus_biomarker.csv"))

# ---- TOC plot ----
mean.ate <- mean(sapply(cv[[1]], function(x) {
  dat <- x[[1]]
  mean(dat$Y[dat$W == 1]) - mean(dat$Y[dat$W == 0])
}))

qq <- seq(0, 1, 0.05)[-1]
qq <- qq[-length(qq)]
qq <- c(qq, 0.995, 1)
get_y <- function(toc) -c(rev(toc$outmat$ad_DropHighFirst), mean.ate) * 100

png_path <- file.path(out_dir, "cate_toc_plot_risk_plus_biomarker.png")
png(png_path, width = 900, height = 650)
par(mar = c(5, 5, 4, 2))
plot(qq, get_y(toc_list[["cf.CATE"]]),
     type = "l", lwd = 2, col = 1, bty = "l", ylim = c(-12, 12),
     xlab = "Proportion treated", ylab = "Absolute Risk Difference (%)",
     main = "CATE TOC — risk score + biomarkers only",
     cex.axis = 1.2, cex.lab = 1.2)
points(qq, get_y(toc_list[["llrf.CATE.rlearn"]]),     type = "l", col = 2, lwd = 2)
points(qq, get_y(toc_list[["bcf.CATE"]]),              type = "l", col = 3, lwd = 2)
points(qq, get_y(toc_list[["pbart.slearn.CATE"]]),     type = "l", col = 4, lwd = 2)
points(qq, get_y(toc_list[["ModCovElasticNet.CATE"]]), type = "l", col = 5, lwd = 2)
abline(h = -mean.ate * 100, lty = 2)
legend("topright", col = 1:5, lty = 1, lwd = 3, bty = "n",
       legend = c("Causal Forest", "Local Linear RF R-learner",
                  "Bayesian Causal Forest", "Probit BART S-learner",
                  "Tian Elastic Net"))
dev.off()
cat("Wrote TOC plot:", png_path, "\n")

# ---- 3-way comparison if previous results exist ----
cat("\n============================================================\n")
cat("3-WAY COMPARISON (if previous passes available)\n")
cat("============================================================\n")

prev_clin <- file.path(OUTPUTS_DIR, "07_cate_modeling", "cate_resmat_clinical.csv")
prev_bio  <- file.path(OUTPUTS_DIR, "07_cate_modeling", "cate_resmat_with_biomarker.csv")

if (file.exists(prev_clin) && file.exists(prev_bio)) {
  rc <- as.matrix(read.csv(prev_clin, row.names = 1))
  rb <- as.matrix(read.csv(prev_bio,  row.names = 1))
  rp <- resmat

  cat("\nAUTOC by method (higher = more HTE detected):\n")
  cat(sprintf("%-25s %8s %8s %8s\n", "Method", "Clinical", "+Bio", "Risk+Bio"))
  for (j in 1:ncol(rc)) {
    cat(sprintf("%-25s %8.2f %8.2f %8.2f\n",
        colnames(rc)[j], rc["autoc", j], rb["autoc", j], rp["autoc", j]))
  }

  cat("\nKendall tau by method:\n")
  cat(sprintf("%-25s %8s %8s %8s\n", "Method", "Clinical", "+Bio", "Risk+Bio"))
  for (j in 1:ncol(rc)) {
    cat(sprintf("%-25s %8.2f %8.2f %8.2f\n",
        colnames(rc)[j], rc["kendall tau", j], rb["kendall tau", j], rp["kendall tau", j]))
  }

  # Save combined comparison
  headline <- data.frame(
    metric         = rownames(rc),
    cf_clin        = rc[, "cf"],       cf_bio = rb[, "cf"],       cf_risk_bio = rp[, "cf"],
    llrf_clin      = rc[, "llrf.r"],   llrf_bio = rb[, "llrf.r"], llrf_risk_bio = rp[, "llrf-r"],
    bcf_clin       = rc[, "bcf"],      bcf_bio = rb[, "bcf"],     bcf_risk_bio = rp[, "bcf"],
    pbart_clin     = rc[, "pbart"],    pbart_bio = rb[, "pbart"],  pbart_risk_bio = rp[, "pbart"],
    enet_clin      = rc[, "Tian.elastic.net"], enet_bio = rb[, "Tian.elastic.net"],
    enet_risk_bio  = rp[, "Tian elastic net"],
    stringsAsFactors = FALSE
  )
  write.csv(headline, file.path(out_dir, "cate_headline_comparison_3way.csv"),
            row.names = FALSE)
  cat("\nSaved 3-way comparison to: cate_headline_comparison_3way.csv\n")
} else {
  cat("Previous pass results not found in outputs/07_cate_modeling/.\n")
  cat("Run 07b_cate_modeling_with_risk.R for the full 3-pass comparison.\n")
}

cat("\nAll outputs written to:\n  ", out_dir, "\n")

sink()
close(report_con)
cat("\nDone. Report saved to:", report_path, "\n")