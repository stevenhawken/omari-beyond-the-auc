# =============================================================================
# real_data.R
# Helpers for the open Zambia / Bangladesh newborn screening data
# (Dryad, doi:10.5061/dryad.m37pvmd6b) and the saved Ontario models.
# Relies on Steve's original code/ga_functions.R and code/gavars.R.
# =============================================================================

COHORTS <- data.frame(
  file    = c("bang_heel.csv", "bang_cord.csv", "zambia_heel.csv", "zambia_cord.csv"),
  country = c("Bangladesh", "Bangladesh", "Zambia", "Zambia"),
  sample  = c("Heel", "Cord", "Heel", "Cord"),
  stringsAsFactors = FALSE
)
COHORTS$cohort <- paste(COHORTS$country, tolower(COHORTS$sample))

# Read one Dryad file into the structure used by the original pipeline
load_cohort <- function(file, country, sample, dir = "data") {
  d <- read.csv(file.path(dir, file), check.names = FALSE)
  names(d)[names(d) == "_17OHP"] <- "N_17OHP"
  names(d) <- sub("^HGB___", "HGB_", names(d))
  d$study_id <- d$Study_ID
  d$multiple <- ifelse(d$multiple == 1, "Yes", "No")   # gavars() expects Yes/No
  d <- as.data.frame(gavars(d))                         # Steve's derived variables
  id <- as.character(d$study_id)
  d$pair_id <- paste0(substr(id, 1, nchar(id) - 2), substr(id, nchar(id), nchar(id)))
  d$country <- country
  d$sample  <- sample
  d$cohort  <- paste(country, tolower(sample))
  d$preterm01 <- as.integer(d$gestage < 37)
  d
}

load_all_cohorts <- function(dir = "data") {
  out <- lapply(seq_len(nrow(COHORTS)), function(i)
    load_cohort(COHORTS$file[i], COHORTS$country[i], COHORTS$sample[i], dir))
  stats::setNames(out, COHORTS$cohort)
}

# Original (published) scoring: prep.data() standardises the batch on itself,
# missing values set to the batch mean (a simple stand-in for the multiple
# imputation used in the papers), then predict. Knots come from the batch too.
score_site_standardised <- function(d, model, alpha = NULL, seed = 1) {
  set.seed(seed)                      # prep.data() adds a tiny random jitter
  p <- prep.data(d)
  num <- vapply(p, is.numeric, logical(1))
  p[num] <- lapply(p[num], function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x })
  p$Sex <- factor(p$Sex, levels = c("F", "M"))
  p$multiple <- factor(p$multiple, levels = c("multiple", "singleton"))
  pr <- if (is.null(alpha)) stats::predict(model, newdata = p) else
    stats::predict(model, newdata = p, alpha = alpha)
  as.numeric(pr)
}

# Performance with bootstrap percentile CIs (predictions held fixed)
boot_perf <- function(obs, pred, B = 1000, seed = 2026) {
  f <- function(o, p) c(MAE = mean(abs(p - o)), bias = mean(p - o),
                        within_7d = mean(abs(p - o) <= 1),
                        MAE_preterm = mean(abs(p - o)[o < 37]),
                        MAE_term = mean(abs(p - o)[o >= 37]))
  est <- f(obs, pred)
  set.seed(seed)
  bs <- replicate(B, { i <- sample(length(obs), replace = TRUE); f(obs[i], pred[i]) })
  lo <- apply(bs, 1, stats::quantile, 0.025, na.rm = TRUE)
  hi <- apply(bs, 1, stats::quantile, 0.975, na.rm = TRUE)
  data.frame(metric = names(est), estimate = est, lower = lo, upper = hi,
             row.names = NULL)
}

fmt_ci <- function(est, lo, hi, d = 2, pct = FALSE) {
  if (pct) { est <- 100 * est; lo <- 100 * lo; hi <- 100 * hi }
  sprintf(paste0("%.", d, "f (%.", d, "f to %.", d, "f)"), est, lo, hi)
}

# Variable lists used by prep.data() and the saved Ontario models
LOG_VARS <- c("N_17OHP", "TREC", "C4DC", "TYR", "C5DC", "Ala", "Arg", "bio",
              "C0", "C2", "C3", "C4", "C5", "C6", "C8", "C10", "C12", "C14",
              "C16", "C18", "C10_1", "C12_1", "C14OH", "C14_1", "C14_2",
              "C16OH", "C16_1OH", "C18OH", "C18_1", "C18_2", "C18_1OH", "C3DC",
              "C4OH", "C5OH", "C5_1", "C6DC", "C8_1", "Cit", "GALT", "Gly",
              "IRT", "LEU", "MET", "Orn", "PHE", "TSH", "Val")
NOLOG_VARS  <- c("af_ratio", "bweight")
SPLINE_VARS <- c("bweight", "af_ratio", "N_17OHP", "TYR", "C5", "C4DC",
                 "C5DC", "TREC", "Ala")

COL_COUNTRY <- c(Bangladesh = "#eb6834", Zambia = "#1baf7a")
SHAPE_SAMPLE <- c(Heel = 16, Cord = 1)
