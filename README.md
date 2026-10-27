# Beyond the AUC: developing and validating clinical prediction models

**OMARI Seminar Series · Wednesday 28 October 2026 · 5:00–7:00 pm**
Steven Hawken, PhD · Senior Scientist, OHRI · Associate Professor, School of Epidemiology and Public Health, University of Ottawa

A model can have an impressive AUC and still give badly wrong individual estimates, fail in a new population, or offer no clinical benefit. This seminar shows how to develop and evaluate a prediction model properly and honestly, and how to read the performance metrics reported in papers.

The case study is a set of machine learning models that estimate gestational age (GA) from routine newborn screening blood tests. The models were developed on Ontario data and externally validated internationally. Using the real models and open data, we look at why performance changes when a model travels and what to do about it.

This repository has everything you need to follow along and rerun every example.

## Before the seminar (10 minutes)

You don't need to run any code to follow the seminar. If you'd like to:

1. **Read ahead:** open `demo.html` and `demo_real.html` in a browser. They are the full companion documents, with all code, output and figures.
2. **Set up R:** install [R](https://cran.r-project.org/) (4.3 or later), [RStudio](https://posit.co/download/rstudio-desktop/) or Positron, and [Quarto](https://quarto.org/docs/get-started/).
3. **Get the files:** click the green **Code** button above, then **Download ZIP** (or `git clone` the repository), and unzip it.
4. **Install packages:** open the folder in RStudio and run `source("install_packages.R")`.
5. **Check it works:** run `quarto::quarto_render("demo_real.qmd")`. It takes about 30 seconds.

`demo.qmd` takes about 5 minutes the first time, because it simulates the data and fits the models. It caches the fits in `output/` so later runs are quick.

## What's here

| Path | What it is |
|---|---|
| `slides/` | Seminar slides (PDF) |
| `demo.qmd` / `demo.html` | **Part 1, simulated data.** Accuracy for a continuous outcome, discrimination, calibration, net benefit, class-imbalance weighting, internal validation (bootstrap optimism), recalibration, and the effect of an imperfect reference standard |
| `demo_real.qmd` / `demo_real.html` | **Part 2, real open data.** The Zambia and Bangladesh cohorts scored with the fitted Ontario models: exploring the data, external validation with bootstrap CIs, small-for-gestational-age infants, heel vs cord blood, preterm discrimination, batch dependence of the original pipeline, frozen preprocessing, and lessons learned |
| `R/` | Helpers for the companion documents: data simulation, model fitting, performance measures written out in base R, plotting |
| `code/` | The original model code, the two fitted Ontario models (`*.rds`), and the frozen-preprocessing fix (`ga_preprocess.R`, `demo_frozen_preprocessing.R`, `FROZEN_PREPROCESSING_NOTES.md`) |
| `data/` | The four open Zambia and Bangladesh cohort files and their data dictionary. The simulated data are generated automatically the first time you render `demo.qmd` |

## Seminar outline

1. **Start with the question:** what is the model for? That decides how to develop and evaluate it.
2. **Accuracy for continuous predictions:** error metrics, calibration plots, regression to the mean, and subgroups such as small-for-gestational-age infants.
3. **Discrimination:** the AUC and what it hides, choosing thresholds, and why rare outcomes kill "accuracy".
4. **Calibration:** is a predicted risk of 20% really 20%? Common failures and how to fix them.
5. **Clinical usefulness:** net benefit and decision curves.
6. **When models travel:** internal vs external validation; shifts in labs, samples and populations; model updating; and pitfalls hidden in the modelling pipeline.
7. **Discussion and wrap-up.**

The session suits researchers, trainees and clinicians comfortable with basic statistics. No engineering or coding background is needed.

## Data

- **Simulated data (Part 1).** Every number in `demo.qmd` comes from simulated data loosely inspired by the published GA work. The relationships are invented, so the results do not reproduce any published study.
- **Zambia and Bangladesh data (Part 2).** Openly available on Dryad: [doi:10.5061/dryad.m37pvmd6b](https://doi.org/10.5061/dryad.m37pvmd6b). Please cite the Dryad dataset and Hawken et al., PLOS ONE 2023 if you use them.
- **Ontario data.** The Ontario development data cannot be shared. The fitted Ontario models in `code/` contain model coefficients and settings only, not individual-level data.

## Key references

- Hawken S, et al. Development and external validation of machine learning algorithms for postnatal gestational age estimation using clinical data and metabolomic markers. *PLOS ONE* 2023.
- Riley RD, Archer L, et al. Evaluation of clinical prediction models (parts 1–3). *BMJ* 2024.
- Van Calster B, et al. Calibration: the Achilles heel of predictive analytics. *BMC Med* 2019.
- Vickers AJ, Elkin EB. Decision curve analysis. *Med Decis Making* 2006.
- van den Goorbergh R, et al. The harm of class imbalance corrections for risk prediction models. *JAMIA* 2022.
- Collins GS, et al. TRIPOD+AI statement. *BMJ* 2024.

Fuller reference lists are at the end of each companion document.

## Licence

Code is released under the MIT licence (see `LICENSE`). Slides and companion documents: © Steven Hawken, shared for teaching. The Dryad data keep their original licence.

Questions: Steven Hawken · shawken@ohri.ca
