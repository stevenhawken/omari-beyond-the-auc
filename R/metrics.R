# =============================================================================
# metrics.R
# Small, transparent implementations of the performance measures used in the
# seminar. Written in base R on purpose so every formula is visible; for real
# analyses see e.g. {pmcalibration}, {CalibrationCurves}, {pROC}, {dcurves},
# {rms::val.prob}.
# =============================================================================

# ---- Continuous outcome (predicted vs reference GA, in weeks) ---------------
cont_metrics <- function(obs, pred) {
  err <- pred - obs
  cal <- coef(lm(obs ~ pred))
  data.frame(
    n          = length(obs),
    bias_wk    = mean(err),                        # mean(predicted - observed)
    MAE_wk     = mean(abs(err)),
    RMSE_wk    = sqrt(mean(err^2)),
    R2         = 1 - sum(err^2) / sum((obs - mean(obs))^2),
    within_7d  = mean(abs(err) <= 1),
    within_14d = mean(abs(err) <= 2),
    cal_int    = unname(cal[1]),
    cal_slope  = unname(cal[2])
  )
}

# ---- Binary outcome: discrimination ----------------------------------------
# AUC = probability that a randomly chosen event (preterm) has a higher
# predicted risk than a randomly chosen non-event (Mann-Whitney form).
auc <- function(y, p) {
  r  <- rank(p)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# DeLong standard error and 95% CI for a single AUC
auc_ci <- function(y, p, level = 0.95) {
  x1 <- p[y == 1]; x0 <- p[y == 0]
  n1 <- length(x1); n0 <- length(x0)
  a  <- auc(y, p)
  # structural components
  v10 <- vapply(x1, function(s) mean((s > x0) + 0.5 * (s == x0)), numeric(1))
  v01 <- vapply(x0, function(s) mean((x1 > s) + 0.5 * (x1 == s)), numeric(1))
  se <- sqrt(var(v10) / n1 + var(v01) / n0)
  z  <- qnorm(1 - (1 - level) / 2)
  c(AUC = a, lower = max(0, a - z * se), upper = min(1, a + z * se))
}

roc_df <- function(y, p) {
  # vectorised via cumulative sums over sorted predictions
  o <- order(p, decreasing = TRUE)
  ys <- y[o]; ps <- p[o]
  tp <- cumsum(ys == 1); fp <- cumsum(ys == 0)
  keep <- c(diff(ps) != 0, TRUE)
  data.frame(fpr = c(0, fp[keep] / sum(y == 0)),
             tpr = c(0, tp[keep] / sum(y == 1)))
}

# ---- Binary outcome: classification at a threshold -------------------------
class_metrics <- function(y, p, threshold) {
  pos <- p >= threshold
  tp <- sum(pos & y == 1); fp <- sum(pos & y == 0)
  fn <- sum(!pos & y == 1); tn <- sum(!pos & y == 0)
  data.frame(threshold = threshold,
             prevalence = mean(y),
             sensitivity = tp / (tp + fn),
             specificity = tn / (tn + fp),
             PPV = tp / (tp + fp),
             NPV = tn / (tn + fn),
             accuracy = (tp + tn) / length(y))
}

# ---- Binary outcome: overall accuracy and calibration -----------------------
clamp <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)

bin_metrics <- function(y, p) {
  lp <- qlogis(clamp(p))
  brier <- mean((p - y)^2)
  brier_null <- mean((mean(y) - y)^2)    # always predicting the prevalence
  # Calibration-in-the-large: intercept with the linear predictor as offset
  citl  <- coef(glm(y ~ offset(lp), family = binomial))[1]
  # Calibration slope: coefficient on the linear predictor
  slope <- coef(glm(y ~ lp, family = binomial))[2]
  ci <- auc_ci(y, p)
  data.frame(
    n = length(y),
    prevalence = mean(y),
    mean_pred = mean(p),
    O_E = mean(y) / mean(p),
    CITL = unname(citl),
    cal_slope = unname(slope),
    AUC = ci[["AUC"]], AUC_lo = ci[["lower"]], AUC_hi = ci[["upper"]],
    Brier = brier,
    scaled_Brier = 1 - brier / brier_null
  )
}

# Flexible calibration curve: logistic regression of the outcome on a natural
# spline of the logit of predicted risk (smoother and better behaved than
# loess when most predicted risks are close to 0)
cal_curve <- function(y, p, df = 3, grid = NULL) {
  lp <- qlogis(clamp(p))
  fit <- glm(y ~ splines::ns(lp, df = df), family = binomial)
  if (is.null(grid)) grid <- seq(quantile(p, 0.01), quantile(p, 0.99),
                                 length.out = 100)
  grid <- grid[grid >= min(p) & grid <= max(p)]
  data.frame(pred = grid,
             obs = predict(fit, newdata = data.frame(lp = qlogis(clamp(grid))),
                           type = "response"))
}

# ---- Clinical usefulness: net benefit / decision curve ----------------------
# Net benefit at threshold pt = TP/n - FP/n * pt / (1 - pt)
net_benefit <- function(y, p, thresholds = seq(0.01, 0.60, by = 0.01)) {
  n <- length(y); prev <- mean(y)
  do.call(rbind, lapply(thresholds, function(pt) {
    pos <- p >= pt
    tp <- sum(pos & y == 1); fp <- sum(pos & y == 0)
    w <- pt / (1 - pt)
    data.frame(threshold = pt,
               model = tp / n - fp / n * w,
               treat_all = prev - (1 - prev) * w,
               treat_none = 0)
  }))
}

# ---- Model updating ---------------------------------------------------------
# Continuous: re-estimate intercept (and optionally slope) in local data
recal_continuous <- function(obs, pred, slope = TRUE) {
  fit <- if (slope) lm(obs ~ pred) else lm(I(obs - pred) ~ 1)
  function(newpred) {
    if (slope) drop(cbind(1, newpred) %*% coef(fit)) else newpred + coef(fit)
  }
}

# Binary: logistic recalibration of the linear predictor
recal_binary <- function(y, p, slope = TRUE) {
  lp <- qlogis(clamp(p))
  fit <- if (slope) glm(y ~ lp, family = binomial) else
    glm(y ~ offset(lp), family = binomial)
  function(newp) {
    nlp <- qlogis(clamp(newp))
    if (slope) plogis(coef(fit)[1] + coef(fit)[2] * nlp) else
      plogis(coef(fit)[1] + nlp)
  }
}
