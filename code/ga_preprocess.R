# =============================================================================
# ga_preprocess.R
# "Learn once, apply everywhere" preprocessing for the metabolic GA models.
#
# Problem this solves
#   In the original pipeline, three things were re-estimated on whatever data
#   were being scored:
#     1. the floor used to replace values <= 0 (minimum positive value),
#     2. the winsorisation cap (Q3 + 3 IQR) and the Pareto mean / sqrt(SD),
#     3. the spline knots, written as knots = quantile(x, ...) inside the model
#        formula (glmnetUtils re-evaluates that expression on new data).
#   As a result, an infant's predicted GA depended on which other infants were
#   in the same batch, and a single infant could not be scored on its own.
#
# The fix
#   Estimate every data-dependent constant ONCE in the reference (training)
#   data, store the constants with the model, and apply them unchanged to any
#   new data. Same idea as recipes::prep() / bake() in tidymodels.
#
# Main functions
#   ga_prep_fit()          learn floor, cap, centre, scale for each variable
#   ga_prep_apply()        apply stored constants to new data (+ shift diagnostics)
#   ga_knots_fit()         learn spline knots / boundary knots on PREPARED data
#   ga_formula()           build a model formula with the knots written in as numbers
#   ga_freeze_terms()      retrofit an existing glmnetUtils model to fixed knots
#   ga_bundle(), predict_ga()   store model + constants together; score new data
#   check_batch_invariance()    unit test: same infant, same prediction, any batch
#   ga_prep_export()       write the constants to CSV (for documentation / vetting)
#
# Requires: splines, glmnetUtils (for the model), base R otherwise.
# =============================================================================

# ---- 1. Learn preprocessing constants in the reference data ------------------
#' @param data       reference (training) data, RAW analyte scale
#' @param log_vars   variables that are log-transformed (as in pareto_log)
#' @param nolog_vars variables that are not log-transformed (as in pareto_nolog)
#' @param winsor_k   upper winsorisation at Q3 + k * IQR (on the log scale for log_vars)
#' @param scaling    "pareto": (x - mean) / sqrt(sd)   [original pipeline]
#'                   "unit":   (x - mean) / sd
#'                   "none":   no centring or scaling
#' @param floor_values optional named vector of fixed floors (e.g. half the
#'                   laboratory limit of detection). Default: minimum positive
#'                   reference value, which reproduces the original pipeline.
ga_prep_fit <- function(data, log_vars, nolog_vars = character(0),
                        winsor_k = 3, scaling = c("pareto", "unit", "none"),
                        floor_values = NULL, reference_label = NULL) {
  scaling <- match.arg(scaling)
  vars <- c(log_vars, nolog_vars)
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars)) stop("Not in data: ", paste(missing_vars, collapse = ", "))

  params <- do.call(rbind, lapply(vars, function(v) {
    x <- as.numeric(data[[v]])
    is_log <- v %in% log_vars
    flo <- if (!is.null(floor_values) && v %in% names(floor_values)) {
      unname(floor_values[[v]])
    } else {
      min(x[x > 0], na.rm = TRUE)
    }
    x[!is.na(x) & x <= 0] <- flo
    if (is_log) x <- log(x)
    q <- stats::quantile(x, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
    cap <- q[2] + winsor_k * (q[2] - q[1])
    x <- pmin(x, cap)
    centre <- if (scaling == "none") 0 else mean(x, na.rm = TRUE)
    s <- stats::sd(x, na.rm = TRUE)
    scale <- switch(scaling, pareto = sqrt(s), unit = s, none = 1)
    z <- (x - centre) / scale
    data.frame(variable = v, log = is_log, floor = flo, cap = cap,
               centre = centre, scale = scale,
               ref_mean = mean(x, na.rm = TRUE), ref_sd = s,
               ref_min = min(z, na.rm = TRUE), ref_max = max(z, na.rm = TRUE),
               ref_p005 = stats::quantile(z, 0.005, na.rm = TRUE, names = FALSE),
               ref_p995 = stats::quantile(z, 0.995, na.rm = TRUE, names = FALSE),
               n = sum(!is.na(x)), pct_missing = mean(is.na(x)),
               stringsAsFactors = FALSE)
  }))
  structure(list(params = params, scaling = scaling, winsor_k = winsor_k,
                 n_ref = nrow(data), reference = reference_label,
                 created = Sys.time()),
            class = "ga_prep")
}

print.ga_prep <- function(x, ...) {
  cat("<ga_prep> ", nrow(x$params), " variables; scaling = ", x$scaling,
      "; winsorisation at Q3 + ", x$winsor_k, " IQR; reference n = ", x$n_ref,
      if (!is.null(x$reference)) paste0(" (", x$reference, ")"), "\n", sep = "")
  print(utils::head(x$params[, c("variable", "log", "floor", "cap", "centre",
                                 "scale")], 10), row.names = FALSE, digits = 4)
  if (nrow(x$params) > 10) cat("  ... ", nrow(x$params) - 10, " more\n", sep = "")
  invisible(x)
}

# ---- 2. Apply stored constants to any data -----------------------------------
#' @param clamp   "none": leave values outside the reference range as they are
#'                "range": clamp to the reference min/max (keeps B-splines from
#'                  extrapolating as cubic polynomials)
#' @param impute  "none" or "reference_mean" (sets NA to the reference mean,
#'                i.e. 0 on the centred scale). Use proper multiple imputation
#'                for real analyses; this is a simple, batch-independent default.
#' Returns the data with the variables replaced by their prepared values and a
#' "shift" attribute: for each variable, the mean of the prepared values in the
#' new data (in reference SD units) and the % outside the reference range.
ga_prep_apply <- function(data, prep, clamp = c("none", "range"),
                          impute = c("none", "reference_mean")) {
  clamp <- match.arg(clamp); impute <- match.arg(impute)
  stopifnot(inherits(prep, "ga_prep"))
  p <- prep$params
  miss <- setdiff(p$variable, names(data))
  if (length(miss)) stop("Not in new data: ", paste(miss, collapse = ", "))
  shift <- vector("list", nrow(p))
  for (i in seq_len(nrow(p))) {
    v <- p$variable[i]
    x <- as.numeric(data[[v]])
    x[!is.na(x) & x <= 0] <- p$floor[i]
    if (p$log[i]) x <- log(x)
    x <- pmin(x, p$cap[i])
    z <- (x - p$centre[i]) / p$scale[i]
    shift[[i]] <- data.frame(
      variable = v,
      # standardised mean difference vs reference (in reference SD units)
      smd_vs_reference = (mean(x, na.rm = TRUE) - p$ref_mean[i]) / p$ref_sd[i],
      pct_below_ref = mean(z < p$ref_min[i], na.rm = TRUE),
      pct_above_ref = mean(z > p$ref_max[i], na.rm = TRUE),
      pct_missing = mean(is.na(z)))
    if (clamp == "range") z <- pmin(pmax(z, p$ref_min[i]), p$ref_max[i])
    if (impute == "reference_mean") z[is.na(z)] <- (p$ref_mean[i] - p$centre[i]) / p$scale[i]
    data[[v]] <- z
  }
  attr(data, "shift") <- do.call(rbind, shift)
  data
}

# Convenience: the covariate-shift table produced by ga_prep_apply()
ga_shift <- function(prepped) {
  s <- attr(prepped, "shift")
  s[order(-abs(s$smd_vs_reference)), ]
}

# ---- 3. Learn spline knots on the PREPARED reference data --------------------
#' @param type "bs_range": interior knots at `probs`, boundary knots at the
#'               reference min/max. Reproduces bs(x, knots = quantile(x, probs))
#'               as used in the original models (needed to retrofit them).
#'             "rcs": restricted (natural) cubic spline with knots at `probs`:
#'               outer knots become the boundary knots, inner knots are interior.
#'               ns(x, knots = inner, Boundary.knots = outer) spans the same
#'               function space as rms::rcs(x, knots) and is LINEAR beyond the
#'               outer knots, which is safer for external data.
ga_knots_fit <- function(prepped_ref, spline_vars,
                         probs = c(0.05, 0.275, 0.5, 0.725, 0.95),
                         type = c("bs_range", "rcs")) {
  type <- match.arg(type)
  out <- lapply(spline_vars, function(v) {
    x <- prepped_ref[[v]]
    q <- stats::quantile(x, probs, na.rm = TRUE, names = FALSE)
    if (type == "bs_range") {
      list(fun = "bs", knots = q, Boundary.knots = range(x, na.rm = TRUE))
    } else {
      list(fun = "ns", knots = q[-c(1, length(q))],
           Boundary.knots = q[c(1, length(q))])
    }
  })
  structure(stats::setNames(out, spline_vars), class = "ga_knots",
            probs = probs, type = type)
}

print.ga_knots <- function(x, ...) {
  cat("<ga_knots> type = ", attr(x, "type"), "; probs = ",
      paste(attr(x, "probs"), collapse = ", "), "\n", sep = "")
  for (v in names(x)) {
    cat(sprintf("  %-9s %s(knots = [%s], boundary = [%s])\n", v, x[[v]]$fun,
                paste(sprintf("%.4g", x[[v]]$knots), collapse = ", "),
                paste(sprintf("%.4g", x[[v]]$Boundary.knots), collapse = ", ")))
  }
  invisible(x)
}

# Build a spline call with numeric knots, e.g. bs(bweight, knots = c(...), Boundary.knots = c(...))
.spline_call <- function(v, k) {
  as.call(list(as.name(k$fun), as.name(v), knots = k$knots,
               Boundary.knots = k$Boundary.knots))
}

# ---- 4. Model formula with the knots written in as numbers --------------------
#' Mirrors the structure of fmla5: spline terms + (linear terms)^2.
#' Because the knots are constants, predict() on new data uses exactly the
#' same basis functions the model was trained with.
ga_formula <- function(knots, linear_vars, outcome = "gestage",
                       interactions = TRUE) {
  spline_terms <- lapply(names(knots), function(v) .spline_call(v, knots[[v]]))
  lin <- Reduce(function(a, b) call("+", a, b), lapply(linear_vars, as.name))
  lin <- if (interactions) call("^", call("(", lin), 2) else lin
  rhs <- Reduce(function(a, b) call("+", a, b), c(spline_terms, list(lin)))
  f <- call("~", as.name(outcome), rhs)
  eval(f, envir = asNamespace("splines"))   # bs/ns resolve from {splines}
}

# ---- 5. Retrofit an existing glmnetUtils model to fixed knots ----------------
#' Replaces every bs(v, knots = quantile(v, ...)) in the stored terms with
#' bs(v, knots = <numbers>, Boundary.knots = <numbers>). Coefficients are
#' untouched. `knots` must come from ga_knots_fit(<the data the model was
#' TRAINED on>, type = "bs_range"), so the basis is identical to training.
ga_freeze_terms <- function(model, knots) {
  n_replaced <- 0
  walk <- function(e) {
    if (is.call(e)) {
      fn <- as.character(e[[1]])[1]
      if (fn %in% c("bs", "ns") && is.name(e[[2]]) &&
          as.character(e[[2]]) %in% names(knots)) {
        v <- as.character(e[[2]])
        k <- knots[[v]]
        if (k$fun != fn) stop("Knot type for ", v, " (", k$fun,
                              ") does not match model term (", fn, ")")
        n_replaced <<- n_replaced + 1
        return(.spline_call(v, k))
      }
      for (j in seq_along(e)[-1]) if (!is.null(e[[j]])) e[[j]] <- walk(e[[j]])
    }
    e
  }
  model$terms <- walk(model$terms)
  remaining <- grep("quantile\\(", deparse(model$terms, width.cutoff = 500L),
                    value = TRUE)
  if (length(remaining)) stop("Some spline terms still use data-dependent knots.")
  attr(model, "frozen_knots") <- knots
  message("Froze ", n_replaced, " spline term(s).")
  model
}

# ---- 6. Bundle model + constants; score new data ------------------------------
#' clamp: default "range" for B-spline (bs) models. bs() extrapolates beyond its
#' boundary knots as unconstrained cubic polynomials, so any external value
#' outside the reference range can produce absurd predictions. Natural splines
#' (type = "rcs") are linear beyond the outer knots and are much safer.
ga_bundle <- function(model, prep, knots, alpha = NULL,
                      factor_levels = list(Sex = c("F", "M"),
                                           multiple = c("multiple", "singleton")),
                      clamp = NULL, meta = list()) {
  if (is.null(clamp)) {
    clamp <- if (any(vapply(knots, function(k) k$fun == "bs", logical(1)))) "range" else "none"
  }
  structure(list(model = model, prep = prep, knots = knots, alpha = alpha,
                 factor_levels = factor_levels, clamp = clamp,
                 meta = c(meta, list(created = Sys.time(),
                                     R = R.version.string,
                                     glmnet = as.character(utils::packageVersion("glmnet")),
                                     glmnetUtils = as.character(utils::packageVersion("glmnetUtils"))))),
            class = "ga_bundle")
}

#' @param prepared set TRUE if `newdata` has already been through
#'   ga_prep_apply() (e.g. prepared once, then multiply imputed with
#'   aregImpute); the stored constants are then not applied a second time.
predict_ga <- function(bundle, newdata, impute = "none", prepared = FALSE, ...) {
  p <- if (prepared) newdata else
    ga_prep_apply(newdata, bundle$prep, clamp = bundle$clamp, impute = impute)
  # B-spline terms: keep inputs inside the boundary knots used at training
  if (identical(bundle$clamp, "range")) {
    for (v in names(bundle$knots)) {
      k <- bundle$knots[[v]]
      if (k$fun == "bs" && v %in% names(p)) {
        p[[v]] <- pmin(pmax(p[[v]], k$Boundary.knots[1]), k$Boundary.knots[2])
      }
    }
  }
  for (f in names(bundle$factor_levels)) {
    if (f %in% names(p)) p[[f]] <- factor(p[[f]], levels = bundle$factor_levels[[f]])
  }
  pr <- if (is.null(bundle$alpha)) {
    stats::predict(bundle$model, newdata = p, ...)
  } else {
    stats::predict(bundle$model, newdata = p, alpha = bundle$alpha, ...)
  }
  out <- as.numeric(pr)
  attr(out, "shift") <- attr(p, "shift")
  out
}

# ---- 7. Unit test: the same infant gets the same prediction in any batch -------
check_batch_invariance <- function(bundle, newdata, n = 50, reps = 5,
                                   seed = 1, impute = "reference_mean", ...) {
  set.seed(seed)
  full <- predict_ga(bundle, newdata, impute = impute, ...)
  diffs <- sapply(seq_len(reps), function(r) {
    idx <- sample(nrow(newdata), min(n, nrow(newdata)))
    sub <- predict_ga(bundle, newdata[idx, , drop = FALSE], impute = impute, ...)
    max(abs(sub - full[idx]), na.rm = TRUE)
  })
  single <- predict_ga(bundle, newdata[1, , drop = FALSE], impute = impute, ...)
  res <- c(max_abs_diff_subsets = max(diffs),
           single_infant_diff = abs(single - full[1]))
  if (all(res < 1e-8)) message("PASS: predictions do not depend on the batch.") else
    warning("FAIL: predictions change with batch composition.")
  res
}

# ---- 8. Export constants (documentation / output vetting) ---------------------
#' Note: `floor`, `ref_min`/`ref_max` and "bs_range" boundary knots are single
#' observed values (minimum / maximum). Organisations such as ICES may treat
#' these as individual-level values; check vetting rules. Using fixed floors
#' (e.g. half the LOD) and percentile-based boundary knots (type = "rcs")
#' avoids reporting extreme individual values.
ga_prep_export <- function(prep, knots = NULL, file) {
  p <- prep$params
  if (!is.null(knots)) {
    k <- do.call(rbind, lapply(names(knots), function(v) {
      data.frame(variable = v, spline = knots[[v]]$fun,
                 knots = paste(signif(knots[[v]]$knots, 10), collapse = ";"),
                 boundary_knots = paste(signif(knots[[v]]$Boundary.knots, 10),
                                        collapse = ";"))
    }))
    p <- merge(p, k, by = "variable", all.x = TRUE, sort = FALSE)
  }
  utils::write.csv(p, file, row.names = FALSE)
  invisible(p)
}
