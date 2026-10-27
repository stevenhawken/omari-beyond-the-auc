# =============================================================================
# plots.R
# ggplot2 helpers shared by the companion document and the slides.
# Colours: validated colour-blind-safe categorical slots (sites and models are
# always shown with the same colour everywhere); every multi-series plot also
# carries a legend or direct labels so identity is never colour alone.
# =============================================================================

library(ggplot2)

COL_SITE  <- c("Development (test)" = "#2a78d6",
               "External site A"    = "#eb6834",
               "External site B"    = "#1baf7a")
COL_MODEL <- c("Full (clinical + analytes)" = "#4a3aa7",
               "Clinical only"              = "#eda100",
               "Random forest"              = "#008300",
               "Small-n OLS"                = "#e34948")
INK   <- "#0b0b0b"
INK2  <- "#52514e"
GRID  <- "#e6e5e0"

SITE_LABEL <- c(test = "Development (test)", siteA = "External site A",
                siteB = "External site B")

theme_seminar <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(colour = INK),
      axis.text = element_text(colour = INK2),
      axis.title = element_text(colour = INK2),
      panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", colour = INK, hjust = 0),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", size = rel(1.1)),
      plot.subtitle = element_text(colour = INK2),
      plot.caption = element_text(colour = INK2, hjust = 0),
      plot.title.position = "plot"
    )
}

# Stack the per-site data frames into one long frame with a site label
stack_sites <- function(res, sites = c("test", "siteA", "siteB")) {
  do.call(rbind, lapply(sites, function(s) {
    d <- res[[s]]
    d$site_lab <- factor(SITE_LABEL[s], levels = SITE_LABEL)
    d
  }))
}

# ---- Continuous: predicted vs observed with loess calibration ---------------
plot_cont_calibration <- function(long, pred = "pred_m3", n_points = 4000,
                                  seed = 1) {
  set.seed(seed)
  pts <- do.call(rbind, lapply(split(long, long$site_lab), function(d) {
    d[sample(nrow(d), min(nrow(d), n_points)), ]
  }))
  ggplot(long, aes(x = .data[[pred]], y = ga_ref)) +
    geom_abline(slope = 1, intercept = 0, colour = INK2, linetype = "dashed") +
    geom_point(data = pts, aes(colour = site_lab), alpha = 0.18, size = 0.8,
               show.legend = FALSE) +
    geom_smooth(method = "loess", se = FALSE, colour = INK, linewidth = 0.9,
                formula = y ~ x, span = 0.75) +
    facet_wrap(~ site_lab) +
    scale_colour_manual(values = COL_SITE) +
    coord_equal(xlim = c(26, 43), ylim = c(26, 43)) +
    labs(x = "Predicted GA (weeks)", y = "Reference GA (weeks)") +
    theme_seminar()
}

# ---- Continuous: error conditional on truth vs on the prediction ------------
# Left: error plotted against reference GA  -> preterm look "overestimated"
# Right: error plotted against predicted GA -> flat (the model is calibrated)
plot_regression_to_mean <- function(d, pred = "pred_m3") {
  d$err <- d[[pred]] - d$ga_ref
  long <- rbind(
    data.frame(x = d$ga_ref, err = d$err,
               panel = "Error vs REFERENCE GA"),
    data.frame(x = d[[pred]], err = d$err,
               panel = "Error vs PREDICTED GA")
  )
  long$panel <- factor(long$panel, levels = unique(long$panel))
  set.seed(2)
  pts <- long[sample(nrow(long), min(nrow(long), 8000)), ]
  ggplot(long, aes(x, err)) +
    geom_hline(yintercept = 0, colour = INK2, linetype = "dashed") +
    geom_point(data = pts, alpha = 0.12, size = 0.7, colour = COL_SITE[1]) +
    geom_smooth(method = "loess", se = FALSE, colour = INK, linewidth = 0.9,
                formula = y ~ x) +
    facet_wrap(~ panel) +
    coord_cartesian(xlim = c(26, 43), ylim = c(-5, 5)) +
    labs(x = "Gestational age (weeks)",
         y = "Predicted - reference GA (weeks)") +
    theme_seminar()
}

# ---- Binary: ROC curves -----------------------------------------------------
plot_roc <- function(roc_long, colour_var, palette, zoom = FALSE) {
  p <- ggplot(roc_long, aes(fpr, tpr, colour = .data[[colour_var]])) +
    geom_abline(colour = INK2, linetype = "dashed") +
    geom_step(linewidth = 0.9, direction = "vh") +
    scale_colour_manual(values = palette) +
    labs(x = "1 - specificity (false-positive rate)",
         y = "Sensitivity (true-positive rate)") +
    theme_seminar()
  if (zoom) p <- p + coord_equal(xlim = c(0, 0.3), ylim = c(0.7, 1)) else
    p <- p + coord_equal()
  p
}

# ---- Binary: flexible calibration curves -------------------------------------
plot_bin_calibration <- function(cal_long, hist_long = NULL, max_p = 1,
                                 colour_var = "site_lab", palette = COL_SITE) {
  p <- ggplot(cal_long, aes(pred, obs, colour = .data[[colour_var]])) +
    geom_abline(colour = INK2, linetype = "dashed") +
    geom_line(linewidth = 0.9)
  if (!is.null(hist_long)) {
    p <- p + geom_rug(data = hist_long, aes(x = p, colour = .data[[colour_var]]),
                      inherit.aes = FALSE, sides = "b", alpha = 0.15,
                      length = grid::unit(0.03, "npc"), show.legend = FALSE)
  }
  p + scale_colour_manual(values = palette) +
    coord_equal(xlim = c(0, max_p), ylim = c(0, max_p)) +
    labs(x = "Predicted risk of preterm birth",
         y = "Observed proportion preterm") +
    theme_seminar()
}

# ---- Decision curve -----------------------------------------------------------
plot_dca <- function(nb_long, ymin = -0.02) {
  ymax <- max(nb_long$net_benefit, na.rm = TRUE) * 1.08
  ggplot(nb_long, aes(threshold, net_benefit, colour = strategy,
                      linetype = strategy)) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = c(COL_MODEL[c("Full (clinical + analytes)",
                                               "Clinical only")],
                                   "Treat all" = INK2, "Treat none" = INK),
                        breaks = c("Full (clinical + analytes)",
                                   "Clinical only", "Treat all",
                                   "Treat none")) +
    scale_linetype_manual(values = c("Full (clinical + analytes)" = "solid",
                                     "Clinical only" = "solid",
                                     "Treat all" = "dashed",
                                     "Treat none" = "dotted"),
                          breaks = c("Full (clinical + analytes)",
                                     "Clinical only", "Treat all",
                                     "Treat none")) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    labs(x = "Threshold probability (risk at which you would act)",
         y = "Net benefit") +
    theme_seminar()
}

# ---- Class imbalance: unweighted vs weighted models -------------------------
COL_WT <- c("Unweighted"                = "#4a3aa7",
            "Graded preterm weights"    = "#008300",
            "Balanced weights (x11)"    = "#e34948")
LTY_WT <- c("Unweighted" = "solid", "Graded preterm weights" = "solid",
            "Balanced weights (x11)" = "longdash")

# Loess of the error against reference GA (left) and predicted GA (right)
plot_weighting_bias <- function(d, preds) {
  long <- do.call(rbind, lapply(names(preds), function(m) {
    err <- d[[preds[[m]]]] - d$ga_ref
    rbind(data.frame(model = m, x = d$ga_ref, err = err,
                     panel = "Error vs REFERENCE GA"),
          data.frame(model = m, x = d[[preds[[m]]]], err = err,
                     panel = "Error vs PREDICTED GA"))
  }))
  long$panel <- factor(long$panel, levels = unique(long$panel))
  long$model <- factor(long$model, levels = names(preds))
  ggplot(long, aes(x, err, colour = model, linetype = model)) +
    geom_hline(yintercept = 0, colour = INK2, linetype = "dashed") +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.9,
                formula = y ~ x, span = 0.6) +
    facet_wrap(~ panel) +
    scale_colour_manual(values = COL_WT) +
    scale_linetype_manual(values = LTY_WT) +
    coord_cartesian(xlim = c(28, 42), ylim = c(-1.5, 1.5)) +
    labs(x = "Gestational age (weeks)",
         y = "Mean of predicted - reference GA (weeks)") +
    theme_seminar()
}
