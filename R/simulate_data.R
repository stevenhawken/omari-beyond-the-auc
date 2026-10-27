# =============================================================================
# simulate_data.R
# Simulated newborn-screening data for estimating gestational age (GA).
#
# TEACHING DATA ONLY. Nothing here is real. The structure (clinical
# covariates + ~30 screening analytes, a development cohort and two external
# sites) is loosely inspired by the metabolic GA literature, but every
# relationship and parameter value below is invented.
#
# Because the data are simulated we keep the TRUE gestational age as well as
# the REFERENCE GA (what a dating ultrasound would record). In real studies
# only the reference is available; having both lets us show how an imperfect
# reference standard changes apparent model performance.
# =============================================================================

# ---- Analyte definitions ----------------------------------------------------
# Each analyte is generated on a latent standard-normal-ish scale z, then
# transformed to a concentration:  value = median * exp(log_sd * z)
# (haemoglobin fractions use a logit scale so they stay within 0-100%).
#
# z = b_lin * g + b_quad * g^2 + b_growth * growth_z + b_age * log(age_h / 48)
#     + b_male * male + (loadings %*% latent factors) + noise,
# where g = (true GA - 39) / 2.

analyte_table <- function() {
  a <- data.frame(
    analyte = c("HbF", "HbA", "TSH", "OHP17", "BIOT", "GALT",
                "ALA", "ARG", "CIT", "GLY", "LEU", "MET", "ORN", "PHE",
                "TYR", "VAL",
                "C0", "C2", "C3", "C3DC", "C4", "C5", "C5DC", "C6", "C8",
                "C10", "C14", "C16", "C18", "C181"),
    group = c("Haemoglobin", "Haemoglobin", "Hormone", "Hormone",
              "Enzyme", "Enzyme", rep("Amino acid", 10),
              rep("Acylcarnitine", 14)),
    median = c(75, 22, 4, 12, 120, 30,
               250, 15, 15, 400, 150, 20, 90, 55, 90, 140,
               25, 20, 1.5, 0.08, 0.25, 0.12, 0.10, 0.06, 0.08, 0.09,
               0.15, 3, 1.0, 1.4),
    scale = c("logit", "logit", rep("log", 28)),
    log_sd = c(0.25, 0.25, 0.55, 0.50, 0.30, 0.30,
               0.35, 0.40, 0.35, 0.30, 0.30, 0.35, 0.35, 0.25, 0.40, 0.30,
               0.30, 0.35, 0.40, 0.45, 0.40, 0.40, 0.45, 0.40, 0.45, 0.45,
               0.40, 0.35, 0.35, 0.35),
    stringsAsFactors = FALSE
  )
  k <- nrow(a)
  # Weak "background" associations for all analytes, then overwrite the
  # analytes that carry most of the GA signal.
  set.seed(20261028)
  a$b_lin    <- round(runif(k, -0.12, 0.12), 2)
  a$b_quad   <- 0
  a$b_growth <- 0
  a$b_age    <- round(runif(k, -0.08, 0.08), 2)
  a$b_male   <- round(runif(k, -0.05, 0.05), 2)

  set_b <- function(a, col, vals) {
    a[match(names(vals), a$analyte), col] <- vals
    a
  }
  a <- set_b(a, "b_lin", c(HbF = -1.00, HbA = 0.90, OHP17 = -0.75,
                           TYR = -0.45, CIT = 0.35, ARG = 0.30, C0 = 0.30,
                           C3 = 0.30, C16 = 0.25, C181 = -0.30,
                           C5DC = -0.25, TSH = 0.20, ALA = -0.20,
                           PHE = -0.20))
  a <- set_b(a, "b_quad", c(HbF = 0.04, OHP17 = -0.05, TYR = -0.04,
                            C0 = -0.03))
  a <- set_b(a, "b_growth", c(C0 = 0.25, C16 = 0.20, TSH = -0.15,
                              HbF = 0.10, ALA = 0.15, C181 = 0.15))
  a <- set_b(a, "b_age", c(OHP17 = -0.45, TSH = -0.55, TYR = 0.30,
                           C0 = -0.20, C2 = 0.20, PHE = 0.10))
  a <- set_b(a, "b_male", c(OHP17 = 0.10, HbF = -0.05))
  a
}

# ---- Site definitions -------------------------------------------------------
# dev   : large high-income development cohort (Ontario-like)
# siteA : external site, lower birthweight for GA (more SGA), higher preterm
#         prevalence, different laboratory calibration, early-ultrasound
#         reference standard
# siteB : external site, similar population shift, a different lab shift, and
#         a reference standard from ultrasounds done at variable timing
site_settings <- function(site) {
  switch(site,
    dev = list(p_pre = 0.062, p_plural = 0.030, bw_shift = 0.00,
               lab_shift = c(),
               us_timing_probs = c("9-13 wk" = 1, "14-24 wk" = 0, ">24 wk" = 0)),
    siteA = list(p_pre = 0.13, p_plural = 0.020, bw_shift = -0.08,
                 lab_shift = c(HbF = 0.35, OHP17 = -0.30, CIT = 0.40,
                               C3 = -0.30, TSH = 0.50, C0 = -0.25),
                 us_timing_probs = c("9-13 wk" = 1, "14-24 wk" = 0, ">24 wk" = 0)),
    siteB = list(p_pre = 0.12, p_plural = 0.025, bw_shift = -0.05,
                 lab_shift = c(HbF = -0.20, OHP17 = 0.25, ARG = 0.35,
                               C16 = -0.30, TYR = 0.25),
                 us_timing_probs = c("9-13 wk" = 0.15, "14-24 wk" = 0.50,
                                     ">24 wk" = 0.35)),
    stop("Unknown site")
  )
}

# SD (weeks) of ultrasound dating error by timing of the dating scan
us_error_sd <- c("9-13 wk" = 0.45, "14-24 wk" = 1.00, ">24 wk" = 1.80)

# Expected log birthweight (Hadlock-type fetal growth curve, GA in weeks)
expected_log_bw <- function(ga, male) {
  0.578 + 0.332 * ga - 0.00354 * ga^2 + 0.035 * male
}
LOG_BW_SD <- 0.11

# ---- True gestational age ---------------------------------------------------
draw_ga <- function(n, p_pre) {
  preterm <- rbinom(n, 1, p_pre)
  # Term births: truncated normal on [37, 42.5)
  lo <- pnorm(37, 39.4, 1.15); hi <- pnorm(42.5, 39.4, 1.15)
  ga_term <- qnorm(runif(n, lo, hi), 39.4, 1.15)
  # Preterm births: left-skewed tail below 37 weeks (mostly late preterm)
  ga_pre <- 37 - rgamma(n, shape = 1.4, rate = 0.50)
  while (any(bad <- ga_pre < 23.5)) {
    ga_pre[bad] <- 37 - rgamma(sum(bad), shape = 1.4, rate = 0.50)
  }
  ifelse(preterm == 1, ga_pre, ga_term)
}

# ---- Main simulator ---------------------------------------------------------
simulate_site <- function(n, site = c("dev", "siteA", "siteB"), seed = 1) {
  site <- match.arg(site)
  set.seed(seed)
  s  <- site_settings(site)
  at <- analyte_table()

  male   <- rbinom(n, 1, 0.51)
  plural <- rbinom(n, 1, s$p_plural)
  # Multiples are much more often preterm
  p_pre_i <- ifelse(plural == 1, 0.55, s$p_pre)
  ga_true <- draw_ga(n, p_pre_i)

  growth_z <- rnorm(n)
  log_bw <- expected_log_bw(ga_true, male) - 0.12 * plural + s$bw_shift +
    LOG_BW_SD * growth_z
  bw <- round(exp(log_bw))
  # SGA relative to the development-population growth standard (sex-specific)
  sga <- as.integer((log_bw - expected_log_bw(ga_true, male)) <
                      qnorm(0.10) * LOG_BW_SD)

  # Heel-prick collection age (hours); mostly 24-72 h
  age_h <- round(pmin(24 + rgamma(n, shape = 2, rate = 0.09), 168))

  # Analytes
  g <- (ga_true - 39) / 2
  n_fac <- 3
  set.seed(99)                       # loadings are the same at every site
  load <- matrix(rnorm(nrow(at) * n_fac, 0, 0.35), nrow(at), n_fac)
  set.seed(seed + 1000)
  fac <- matrix(rnorm(n * n_fac), n, n_fac)
  shift <- setNames(rep(0, nrow(at)), at$analyte)
  if (length(s$lab_shift)) shift[names(s$lab_shift)] <- s$lab_shift

  analytes <- sapply(seq_len(nrow(at)), function(k) {
    z <- at$b_lin[k] * g + at$b_quad[k] * g^2 +
      at$b_growth[k] * growth_z + at$b_age[k] * log(age_h / 48) +
      at$b_male[k] * male + drop(fac %*% load[k, ]) + rnorm(n) +
      shift[k]
    if (at$scale[k] == "logit") {
      100 * plogis(qlogis(at$median[k] / 100) + at$log_sd[k] * z)
    } else {
      at$median[k] * exp(at$log_sd[k] * z)
    }
  })
  colnames(analytes) <- at$analyte
  analytes <- signif(analytes, 4)

  # Reference standard: dating ultrasound with timing-dependent error
  us_timing <- sample(names(s$us_timing_probs), n, replace = TRUE,
                      prob = s$us_timing_probs)
  ga_ref <- ga_true + rnorm(n, 0, us_error_sd[us_timing])
  ga_ref <- round(ga_ref * 7) / 7   # recorded in whole days

  data.frame(
    id = sprintf("%s_%05d", site, seq_len(n)),
    site = site,
    ga_ref = ga_ref,                 # what a real study would observe
    ga_true = round(ga_true, 3),     # only knowable in a simulation
    us_timing = us_timing,
    preterm = as.integer(ga_ref < 37),
    preterm_true = as.integer(ga_true < 37),
    sga = sga,
    male = male,
    plural = plural,
    bw = bw,
    age_h = age_h,
    analytes,
    stringsAsFactors = FALSE
  )
}

# ---- Build and save the three datasets --------------------------------------
simulate_all <- function(out_dir = "data") {
  dir.create(out_dir, showWarnings = FALSE)
  dev   <- simulate_site(60000, "dev",   seed = 2016)
  siteA <- simulate_site(1500,  "siteA", seed = 2019)
  siteB <- simulate_site(1000,  "siteB", seed = 2022)
  write.csv(dev,   file.path(out_dir, "sim_dev.csv"),   row.names = FALSE)
  write.csv(siteA, file.path(out_dir, "sim_siteA.csv"), row.names = FALSE)
  write.csv(siteB, file.path(out_dir, "sim_siteB.csv"), row.names = FALSE)
  invisible(list(dev = dev, siteA = siteA, siteB = siteB))
}
