# =============================================================================
# fit_models.R
# Fits the teaching models on the simulated development cohort and returns
# predictions for the internal test set and the two external sites.
#
# Models
#   M1  "Clinical"      : sex, multiple birth, birthweight (spline)       [lm]
#   M3  "Full"          : M1 + log analytes + squared log analytes +
#                         collection age                             [elastic net]
#   RF  "Random forest" : same predictors as M3                         [ranger]
#   SMALL "Small-n OLS" : M3 predictors, unpenalised, fitted on n = 200
#                         (an overfitting demonstration)                  [lm]
# =============================================================================

library(glmnet)
library(splines)

ANALYTES <- analyte_table()$analyte

# ---- Feature construction ---------------------------------------------------
# The birthweight spline basis and the centring constants are learned on the
# TRAINING data and re-used unchanged everywhere else (no peeking).
feature_spec <- function(train) {
  la <- log(as.matrix(train[, ANALYTES]))
  list(
    bw_knots  = quantile(train$bw, c(0.25, 0.50, 0.75)),
    bw_bounds = range(train$bw),
    la_centre = colMeans(la)
  )
}

make_x <- function(d, spec, analytes = TRUE) {
  bw <- pmin(pmax(d$bw, spec$bw_bounds[1]), spec$bw_bounds[2])
  bw_ns <- ns(bw, knots = spec$bw_knots, Boundary.knots = spec$bw_bounds)
  colnames(bw_ns) <- paste0("bw_ns", seq_len(ncol(bw_ns)))
  x <- cbind(male = d$male, plural = d$plural, bw_ns)
  if (analytes) {
    la <- sweep(log(as.matrix(d[, ANALYTES])), 2, spec$la_centre)
    la2 <- la^2
    colnames(la2) <- paste0(colnames(la), "_sq")
    x <- cbind(x, log_age = log(d$age_h), la, la2)
  }
  x
}

# ---- Fit everything ---------------------------------------------------------
fit_all <- function(dev, siteA, siteB, seed = 2024) {
  set.seed(seed)
  # Two-thirds / one-third split, as in the Ontario development work
  idx_train <- sample(nrow(dev), size = round(2 / 3 * nrow(dev)))
  train <- dev[idx_train, ]
  test  <- dev[-idx_train, ]
  spec  <- feature_spec(train)

  x_tr1 <- make_x(train, spec, analytes = FALSE)
  x_tr3 <- make_x(train, spec, analytes = TRUE)

  # M1: clinical only (ordinary least squares)
  m1 <- lm(train$ga_ref ~ x_tr1)

  # M3: elastic net, alpha fixed at 0.5, lambda by 10-fold CV
  m3 <- cv.glmnet(x_tr3, train$ga_ref, alpha = 0.5, nfolds = 10)

  # Random forest on the same predictors (a flexible comparator)
  rf <- NULL
  if (requireNamespace("ranger", quietly = TRUE)) {
    sub <- sample(nrow(train), 15000)          # keeps the demo quick
    rf <- ranger::ranger(x = x_tr3[sub, ], y = train$ga_ref[sub],
                         num.trees = 300, min.node.size = 10,
                         seed = seed, num.threads = 2)
    rf$train_rows <- sub
  }

  # SMALL: unpenalised regression on 200 infants -> overfitting
  small_rows <- sample(nrow(train), 200)
  small <- lm.fit(cbind(1, x_tr3[small_rows, ]), train$ga_ref[small_rows])

  models <- list(spec = spec, m1 = m1, m3 = m3, rf = rf, small = small,
                 small_rows = small_rows, train_rows = idx_train)

  list(models = models,
       train = cbind(train, predict_all(models, train)),
       test  = cbind(test,  predict_all(models, test)),
       siteA = cbind(siteA, predict_all(models, siteA)),
       siteB = cbind(siteB, predict_all(models, siteB)))
}

predict_all <- function(models, d) {
  x1 <- make_x(d, models$spec, analytes = FALSE)
  x3 <- make_x(d, models$spec, analytes = TRUE)
  b_small <- models$small$coefficients
  b_small[is.na(b_small)] <- 0
  out <- data.frame(
    pred_m1 = drop(cbind(1, x1) %*% coef(models$m1)),
    pred_m3 = drop(predict(models$m3, newx = x3, s = "lambda.min")),
    pred_small = drop(cbind(1, x3) %*% b_small)
  )
  if (!is.null(models$rf)) {
    out$pred_rf <- predict(models$rf, data = x3, num.threads = 2)$predictions
  }
  out
}

# ---- From predicted GA to a predicted probability of preterm birth ----------
# If the prediction error is roughly normal with SD sigma, then
#   P(GA < 37 | prediction) = Phi((37 - predicted GA) / sigma)
# sigma is estimated from residuals in the TRAINING data only.
preterm_prob <- function(pred, sigma, cut = 37) {
  pnorm((cut - pred) / sigma)
}

# ---- Class imbalance: upweighting preterm infants in training ---------------
# Graded weights used in the original metabolic GA models (translated from SAS):
#   GA >= 37 -> 1;  22 <= GA < 37 -> 1 + (37 - GA) / 1.5
# so an infant born at 34 wk counts 3x, at 28 wk 7x. Weights are used in
# model TRAINING only; evaluation is always unweighted.
preterm_weight <- function(ga) ifelse(ga >= 37, 1, 1 + (37 - ga) / 1.5)

fit_weighted <- function(res, seed = 2025) {
  set.seed(seed)
  tr <- res$train
  x  <- make_x(tr, res$models$spec)
  w_graded <- preterm_weight(tr$ga_ref)
  w_bal <- ifelse(tr$preterm == 1, (1 - mean(tr$preterm)) / mean(tr$preterm), 1)
  folds <- sample(rep(1:10, length.out = nrow(tr)))   # same folds for all fits
  list(
    w_bal      = unique(w_bal[tr$preterm == 1]),
    cont_graded = cv.glmnet(x, tr$ga_ref, weights = w_graded, alpha = 0.5,
                            foldid = folds),
    bin_unw    = cv.glmnet(x, tr$preterm, alpha = 0.5, foldid = folds,
                           family = "binomial"),
    bin_graded = cv.glmnet(x, tr$preterm, weights = w_graded, alpha = 0.5,
                           foldid = folds, family = "binomial"),
    bin_bal    = cv.glmnet(x, tr$preterm, weights = w_bal, alpha = 0.5,
                           foldid = folds, family = "binomial")
  )
}
