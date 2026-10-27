# Frozen preprocessing for the metabolic GA models

## The issue

The original pipeline re-estimates three things on whatever data are being scored:

| Step | Where | What is re-estimated |
|---|---|---|
| Floor for values ≤ 0 | `pareto_log()` / `pareto_nolog()` | minimum positive value |
| Winsorisation and Pareto scaling | same | Q3 + 3 IQR cap; mean and sqrt(SD) |
| Spline knots | model formula `bs(x, knots = quantile(x, ...))` | knots **and** boundary knots (glmnetUtils re-evaluates the expression on new data and does not store the training values) |

As a result, an infant's prediction depends on which other infants are in the batch. On the open Bangladesh heel data, scoring the same 60 infants as their own preterm-enriched batch shifted predictions by +0.86 weeks, and MAE rose from 1.10 to 1.66. A single infant cannot be scored on its own.

## The fix: learn once, apply everywhere

`ga_preprocess.R` stores every data-dependent constant from the reference (training) data and applies it unchanged. This is the same idea as `recipes::prep()` / `bake()`.

| Function | Purpose |
|---|---|
| `ga_prep_fit()` | learn floor, cap, centre and scale per variable (mirrors `pareto_log` exactly; verified to 1e-5) |
| `ga_prep_apply()` | apply them to any data; adds a covariate-shift table (`ga_shift()`) |
| `ga_knots_fit()` | learn knots on the **prepared** reference data: `"bs_range"` (reproduces the current models) or `"rcs"` (restricted cubic, linear tails) |
| `ga_freeze_terms()` | retrofit a saved glmnetUtils model: writes the knots into its formula as numbers; coefficients untouched |
| `ga_formula()` | build a new model formula with numeric knots (for refitting) |
| `ga_bundle()`, `predict_ga()` | keep model + constants together; score raw data |
| `check_batch_invariance()` | unit test: same infant, same prediction, any batch, including a batch of one |
| `ga_prep_export()` | write the constants to CSV |

Verified on the open data with a stand-in reference: retrofitted models reproduce the original predictions exactly (difference 0 on the fitting data; 8.6e-5 from the original `jitter()`), and the batch-invariance check passes.

## What to run at ICES (Part A of `run_ga_models_frozen.R`)

1. Identify the raw training frame that produced `ga_train_pto`. The check in step A3 confirms it.
2. Fit the constants, compute the knots from the rows the model was actually fitted on, freeze both models, and verify (A6 must give zero difference).
3. Report the frozen-pipeline MAE on the Ontario test set next to the published value.
4. Save the bundles and the constants CSV.

**Output vetting.** Exact retrofit needs the `bs()` boundary knots, which are the training minimum and maximum of each spline variable. Moving them to percentiles changed predictions by an average of 16 weeks in testing, so they must be exact. The floors are also single observed values. Ask ICES whether these ~60 numbers can be released. If they cannot, refit with `type = "rcs"`, which only uses percentiles, and fixed LOD-based floors (Part C).

## Recommendations for next time

1. **Treat the whole pipeline as the model.** Learn every preprocessing constant in the training data only, store it with the coefficients, and test batch invariance before any external validation.
2. **Write spline knots in as numbers**, or use a framework that stores them (tidymodels recipes, `rms`).
3. **Prefer restricted (natural) cubic splines for models that will travel.** `bs()` extrapolates beyond its boundary knots as an unconstrained cubic. In testing, Zambia heel MAE was 1.18 (worst error 8.3 weeks) without clamping and 0.96 with it. `ns(x, knots = inner, Boundary.knots = outer)` spans the same function space as `rms::rcs(x, 5 knots)` and is linear beyond the outer knots.
4. **Look at covariate shift before outcomes.** The `ga_shift()` table (standardised mean differences vs the reference) shows where lab or population differences lie.
5. **Choose site-level standardisation on purpose.** Re-centring each site on its own data can absorb laboratory calibration differences, which is a legitimate local-adaptation strategy. It also assumes the site's GA and birthweight mix resembles the reference, and it needs a representative batch. If used, describe it as part of the model and test its sensitivity to batch composition.
6. **Justify and tune the class-imbalance weights** (see the seminar notes), and report calibration.

## Other things noticed in the code (minor)

- `mean(min(x[x > 0], na.rm = TRUE), 0)` returns the minimum positive value, because the `0` is taken as `mean()`'s `trim` argument. If half the minimum was intended, it is not doing that.
- `jitter(x, amount = 0.00001)` makes predictions slightly non-deterministic (≤ 1e-4 weeks). It can be dropped with frozen knots.
- The paper describes "restricted cubic splines with 5 knots"; the code uses `bs()` (unrestricted cubic B-splines, 5 interior knots, boundary knots at the data range).
- In `data_prep.R`, `zambia.data <- gavars(Zambia)` is applied to the unfiltered `Zambia`, so the preceding `filter()` (missing GA, birthweight or sex) has no effect.
