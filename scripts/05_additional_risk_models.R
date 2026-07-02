# ============================================================================
# risk_methods_extra.R
#
# PURPOSE: extra risk-modeling functions that mirror the signature and
# behavior of cate-repo/scoring.methods_v4.R's existing logreg.risk/rf.risk.
#
# These are designed to plug into the same cross-validation framework
# Victor uses (oest.cv.fun in cate-repo/cross.validation.R), so once these
# are sourced alongside scoring.methods_v4.R, you can pass any of them by
# name into the same `mods` mechanism that drives the CATE comparison --
# no other changes to the pipeline needed.
#
# What's added here, all matching logreg.risk's contract exactly
# (signature, control-arm-only fit, MI-aware loop, return shape):
#   enet.risk    : penalized logistic regression via glmnet
#                   (elastic net by default; ridge or LASSO via alpha)
#   xgb.risk     : gradient boosting via xgboost
#   bart.risk    : Bayesian additive regression trees via stochtree
#                   (reuses the dependency already in cate-repo)
#
# Each function returns list(predicted_probabilities_on_eval_set, NA)
# matching logreg.risk's two-element-list output shape, so it slots into
# the same downstream machinery as logreg.risk and rf.risk without any
# special-case handling.
#
# These do NOT live in cate-repo/ because they aren't Victor's; they're
# additions written for CLOVERS specifically. Source them from your driver
# script via:
#   source(file.path(PROJECT_ROOT, "preprocessing", "risk_methods_extra.R"))
# ============================================================================


# ---- shared helper used by all three -----------------------------------
.invlogit <- function(B) exp(B) / (1 + exp(B))


# ============================================================================
# enet.risk: elastic net penalized logistic regression via glmnet
# ----------------------------------------------------------------------------
# `alpha` controls the mixing between ridge (alpha=0) and LASSO (alpha=1).
# Default 0.5 is elastic net. Pass alpha=0 for pure ridge, alpha=1 for
# pure LASSO -- same one function covers all three by changing one arg.
#
# `lambda` is chosen by internal cross-validation (cv.glmnet) on each
# imputation's control-arm-only training fold, which is the standard
# practice for glmnet -- there's no good way to pre-pick lambda by hand.
# ============================================================================
enet.risk <- function(train, cvt, cve, M, y, w, x, pi, alpha = 0.5, ...) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("enet.risk requires the glmnet package -- install.packages('glmnet')")
  }

  iw <- which(colnames(train[[1]]) == w)

  pred.all <- c()
  for (m in 1:M) {
    ds  <- train[[m]]
    dst <- ds[cvt, ]
    i0  <- which(dst[, iw] == 0)   # control arm only
    dst <- dst[i0, ]
    dse <- ds[cve, ]

    X_train <- as.matrix(dst[, x])
    y_train <- dst[, y]
    X_eval  <- as.matrix(dse[, x])

    fit <- glmnet::cv.glmnet(X_train, y_train,
                              family = "binomial", alpha = alpha)
    pred <- predict(fit, newx = X_eval, s = "lambda.min", type = "link")[, 1]

    pred.all <- cbind(pred.all, pred)
  }

  mean.pred <- apply(pred.all, 1, mean)
  list(.invlogit(mean.pred), NA)
}


# ============================================================================
# xgb.risk: gradient boosting via xgboost
# ----------------------------------------------------------------------------
# Uses binary:logistic objective so predictions are already probabilities.
# Hyperparameters are conservative defaults chosen to keep runtime reasonable
# and to avoid catastrophic overfitting on a ~670-patient training fold:
#   nrounds = 200 with early stopping at 20 rounds of no improvement
#   max_depth = 4, eta = 0.05  (shallow trees, slow learning rate)
# Adjust if needed, but these are sensible starting points -- and they
# match what the comparison literature (search results earlier) tends to use
# for clinical risk-prediction benchmarks at this sample size.
# ============================================================================
xgb.risk <- function(train, cvt, cve, M, y, w, x, pi,
                      nrounds = 200, max_depth = 4, eta = 0.05, ...) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("xgb.risk requires the xgboost package -- install.packages('xgboost')")
  }

  iw <- which(colnames(train[[1]]) == w)

  pred.all <- c()
  for (m in 1:M) {
    ds  <- train[[m]]
    dst <- ds[cvt, ]
    i0  <- which(dst[, iw] == 0)
    dst <- dst[i0, ]
    dse <- ds[cve, ]

    dtrain <- xgboost::xgb.DMatrix(as.matrix(dst[, x]), label = dst[, y])
    deval  <- xgboost::xgb.DMatrix(as.matrix(dse[, x]))

    params <- list(
      objective = "binary:logistic",
      max_depth = max_depth,
      eta = eta,
      verbosity = 0
    )

    fit <- xgboost::xgb.train(
      params = params,
      data = dtrain,
      nrounds = nrounds,
      watchlist = list(train = dtrain),
      early_stopping_rounds = 20,
      verbose = 0
    )

    pred <- predict(fit, newdata = deval)   # already on probability scale
    pred.all <- cbind(pred.all, pred)
  }

  mean.pred <- apply(pred.all, 1, mean)
  list(mean.pred, NA)
}


# ============================================================================
# bart.risk: Bayesian additive regression trees via stochtree
# ----------------------------------------------------------------------------
# Uses stochtree::bart() since it's already a dependency in cate-repo for
# pbart.slearn.CATE -- no extra package required if pbart.slearn.CATE
# already runs on your machine.
# ============================================================================
bart.risk <- function(train, cvt, cve, M, y, w, x, pi, ...) {
  if (!requireNamespace("stochtree", quietly = TRUE)) {
    stop("bart.risk requires the stochtree package (already needed for pbart.slearn.CATE)")
  }

  iw <- which(colnames(train[[1]]) == w)

  pred.all <- c()
  for (m in 1:M) {
    ds  <- train[[m]]
    dst <- ds[cvt, ]
    i0  <- which(dst[, iw] == 0)
    dst <- dst[i0, ]
    dse <- ds[cve, ]

    X_train <- as.matrix(dst[, x])
    y_train <- dst[, y]
    X_eval  <- as.matrix(dse[, x])

    fit <- stochtree::bart(
      X_train = X_train,
      y_train = y_train,
      X_test  = X_eval
    )

    # stochtree::bart returns posterior draws of predictions; average them
    # over draws to get the posterior mean probability per eval patient.
    # The element name has historically been y_hat_test in stochtree.
    if (!is.null(fit$y_hat_test)) {
      pred_draws <- fit$y_hat_test
    } else {
      stop("Unexpected stochtree::bart output structure -- expected y_hat_test")
    }
    pred <- rowMeans(pred_draws)

    pred.all <- cbind(pred.all, pred)
  }

  mean.pred <- apply(pred.all, 1, mean)
  # bart's output is already on the probability scale by default in stochtree
  # when y is a 0/1 binary outcome with probit link -- no invlogit needed.
  list(mean.pred, NA)
}