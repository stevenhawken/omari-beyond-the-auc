# =============================================================================
# demo_frozen_preprocessing.R
# Runs on the OPEN Dryad data (doi:10.5061/dryad.m37pvmd6b) and the saved
# Ontario models. Run from the seminar folder root:
#   source("code/demo_frozen_preprocessing.R")
#
# 1. Shows the batch dependence of the original pipeline.
# 2. Demonstrates the frozen pipeline, using Bangladesh cord as a STAND-IN
#    reference (the Ontario constants are not available outside ICES yet):
#    exact reproduction, batch invariance, covariate-shift table, clamping.
# =============================================================================
suppressMessages({
  library(dplyr); library(glmnetUtils); library(splines)
})
source("code/ga_functions.R")
source("code/ga_preprocess.R")
full_model <- readRDS("code/full_model_pto_relax.rds")

# ---- Read a Dryad file into the structure prep.data() expects ----------------
read_dryad <- function(file) {
  d <- read.csv(file, check.names = FALSE)
  names(d)[names(d) == "_17OHP"] <- "N_17OHP"
  names(d) <- sub("^HGB___", "HGB_", names(d))
  d |>
    mutate(study_id = Study_ID,
           af_ratio = HGB_A / (HGB_A + HGB_F_F1),
           multiple = ifelse(multiple == 1, "multiple", "singleton"),
           bweight_raw = bweight, gest_wk = floor(gestage),
           sga3 = NA, sga10 = NA, lga90 = NA, premie = NA,
           preterm = ifelse(gestage < 37, "<37 wks", ">=37 wks"),
           lbw = ifelse(bweight < 2500, "<2500 g", ">=2500 g"))
}
as_factors <- function(p) {
  p$Sex <- factor(p$Sex, levels = c("F", "M"))
  p$multiple <- factor(p$multiple, levels = c("multiple", "singleton"))
  p
}
# Original pipeline: prep.data() on the batch, then predict (knots from the batch)
score_original <- function(d) {
  set.seed(1)
  p <- prep.data(d)
  num <- sapply(p, is.numeric)
  p[num] <- lapply(p[num], function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x })
  as.numeric(predict(full_model, newdata = as_factors(p), alpha = 0.75))
}

LOG_VARS <- c("N_17OHP", "TREC", "C4DC", "TYR", "C5DC", "Ala", "Arg", "bio",
              "C0", "C2", "C3", "C4", "C5", "C6", "C8", "C10", "C12", "C14",
              "C16", "C18", "C10_1", "C12_1", "C14OH", "C14_1", "C14_2",
              "C16OH", "C16_1OH", "C18OH", "C18_1", "C18_2", "C18_1OH", "C3DC",
              "C4OH", "C5OH", "C5_1", "C6DC", "C8_1", "Cit", "GALT", "Gly",
              "IRT", "LEU", "MET", "Orn", "PHE", "TSH", "Val")
NOLOG_VARS  <- c("af_ratio", "bweight")
SPLINE_VARS <- c("bweight", "af_ratio", "N_17OHP", "TYR", "C5", "C4DC",
                 "C5DC", "TREC", "Ala")

bang_heel <- read_dryad("data/bang_heel.csv")

# ---- 1. Batch dependence of the ORIGINAL pipeline ----------------------------
full_batch <- score_original(bang_heel)
pre  <- which(bang_heel$gestage < 37)
term <- which(bang_heel$gestage >= 37)
set.seed(12)
idx <- c(sample(pre, 30), sample(term, 30))     # a preterm-enriched batch
small_batch <- score_original(bang_heel[idx, ])
cat("ORIGINAL pipeline, the same 60 infants:\n")
cat("  mean shift when scored as their own batch:",
    round(mean(small_batch - full_batch[idx]), 2), "weeks\n")
cat("  MAE within full cohort:", round(mean(abs(full_batch[idx] - bang_heel$gestage[idx])), 2),
    "| as own batch:", round(mean(abs(small_batch - bang_heel$gestage[idx])), 2), "\n\n")

# ---- 2. Frozen pipeline with a stand-in reference ----------------------------
ref_raw <- read_dryad("data/bang_cord.csv")     # STAND-IN for Ontario training
prep <- ga_prep_fit(ref_raw, LOG_VARS, NOLOG_VARS,
                    reference_label = "Bangladesh cord (stand-in)")
set.seed(1)
ref_pto <- prep.data(ref_raw)                   # what the original code produced
fit_rows <- complete.cases(ref_pto[, c(LOG_VARS, NOLOG_VARS)])
knots <- ga_knots_fit(ref_pto[fit_rows, ], SPLINE_VARS, type = "bs_range")
frozen <- ga_freeze_terms(full_model, knots)
bundle <- ga_bundle(frozen, prep, knots, alpha = 0.75)

# (a) exact reproduction on the reference data
orig <- as.numeric(predict(full_model, newdata = as_factors(ref_pto[fit_rows, ]), alpha = 0.75))
frz  <- predict_ga(bundle, ref_raw[fit_rows, ])
cat("Frozen vs original on the reference data: max |diff| =",
    signif(max(abs(frz - orig)), 2), "weeks (jitter only)\n")

# (b) batch invariance on external data
print(check_batch_invariance(bundle, bang_heel, n = 60))

# (c) covariate shift of an external cohort relative to the reference
zam <- read_dryad("data/zambia_heel.csv")
cat("\nLargest shifts, Zambia heel vs reference (SD units):\n")
print(head(ga_shift(ga_prep_apply(zam, prep)), 8), row.names = FALSE, digits = 2)

# (d) why clamping matters for B-splines
b_none <- bundle; b_none$clamp <- "none"
p_none <- suppressWarnings(predict_ga(b_none, zam, impute = "reference_mean"))
p_clmp <- predict_ga(bundle, zam, impute = "reference_mean")
cat("\nZambia heel MAE without clamping:", round(mean(abs(p_none - zam$gestage)), 2),
    "(worst error", round(max(abs(p_none - zam$gestage)), 1), "wk);",
    "with clamping:", round(mean(abs(p_clmp - zam$gestage)), 2), "\n")
