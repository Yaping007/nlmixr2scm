## scorecard_stats.R
## -----------------------------------------------------------------------------
## Inferential layer for the scorecard: attaches a p-value to every
## method-vs-reference comparison in the heatmap, using the RAW per-replicate
## data (dataset_id) rather than the pre-aggregated medians.
##
## DESIGN (see docs/METHODS_manuscript.md, "Statistical comparison of methods")
##   * Paired by construction: the same simulated datasets are analysed by all
##     three workflows, so we pair on the key (scenario, dataset_id) and use the
##     dataset_id INTERSECTION shared by the two methods being compared.
##   * Stratified by structure: the model is fit separately within Analytic and
##     within ODE (linCmt/ADVAN4 vs ODE are different experiments; PsN collapses
##     on ODE), giving one CI per heatmap cell.
##   * Reference method = PsN-SCM; every non-reference method gets a
##     "method vs PsN" effect, so the inference matches the heatmap's colour.
##
## PRIMARY analysis  (when lme4 is installed): a mixed-effects model
##   binary   metrics : glmer(y ~ method + (1|scenario) + (1|scenario:dataset_id),
##                             family = binomial)          -> log-odds / OR
##   Runtime          : lmer(log(minutes) ~ method + (1|scenario)
##                             + (1|scenario:dataset_id))  -> log-fold / fold-change
##   Accuracy (MARE)  : lmer(MARE ~ method + (1|scenario) + (1|scenario:dataset_id))
##   The `method` fixed-effect coefficient is the shrunk, dataset-paired MEAN
##   effect of the method vs PsN on the link scale (log-odds for binary,
##   log-fold for runtime, mean difference for MARE), adjusting for how hard
##   each scenario is and respecting the pairing of the three fits of one data.
##
## FALLBACK analysis (used when lme4 is unavailable, e.g. offline): the paired
##   exact tests that the mixed model generalises -
##   binary     : exact McNemar (binomial on the discordant pairs)
##   continuous : paired Wilcoxon signed-rank
##   These are also reported as the sensitivity column even when lme4 is present.
##
## The DISPLAYED cell number stays the descriptive MEDIAN across scenarios
## (see scorecard.R); this script only supplies the SIGNIFICANCE MARKER. Number
## (median, natural units) and marker (model effect, link scale) are different
## estimators of the same "vs PsN" question - they agree in direction but not in
## value; the caption states this explicitly.
##
## OUTPUT: build_scorecard_stats() returns a tidy tibble keyed by
##   (method, structure, metric) with columns:
##     estimate, ci_low, ci_high (effect-size, on the reported scale),
##     effect_scale, p_raw, p_adj (Holm), test, n_pairs, tested, marker, better
## and, with save = TRUE, writes output/figures/scorecard/scorecard_stats.csv.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
})

if (!exists(".SC_SOURCES")) source("script/viz/scorecard.R")

.ST_HAVE_LME4 <- requireNamespace("lme4", quietly = TRUE)

## metrics that are Bernoulli per dataset vs continuous per dataset
.ST_BINARY <- c("Convergence", "CovStep", "Power", "FalsePos")
.ST_CONT   <- c("Accuracy", "Runtime")

## cells we deliberately DO NOT test (structurally non-comparable), reported as
## descriptive only. Runtime IS tested head-to-head for every method: the
## benchmark endpoint is "wall-clock time until covariate selection is complete"
## -- the full stepwise pipeline for SCM, the single joint fit for VAE -- so a
## one-pass VAE is a fair competitor, not an exclusion. No cell is excluded at
## present; the hook is retained for future structurally-incomparable cells.
.st_not_tested <- function(method, structure, metric) {
  FALSE
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

## ---- per-dataset outcome frame for ONE method x structure ------------------
## one row per simulated dataset, all six metrics as columns.

.st_col <- function(d, nm) if (nm %in% names(d)) d[[nm]] else NA

.st_dataset_metrics <- function(method, structure, sample_N = 80) {
  src <- .SC_SOURCES[[method]]
  struct_keep <- if (identical(structure, "Analytic"))
                   c("linCmt", "advan4") else "ode"

  dl <- .sc_cell(.sc_read(src, src$diag_long), src, sample_N, struct_keep)
  wall_col <- if ("wall_total_sec" %in% names(dl)) "wall_total_sec" else
              if ("fit_runtime_sec" %in% names(dl)) "fit_runtime_sec" else NA
  fp <- .st_col(dl, "n_false_pos")

  base <- tibble::tibble(
    scenario    = dl$scenario,
    dataset_id  = dl$dataset_id,
    Convergence = as.integer(.st_col(dl, "converged") %in% TRUE),
    CovStep     = as.integer(.st_col(dl, "cov_ok")    %in% TRUE),
    Power       = as.integer(.st_col(dl, "exact_match") %in% TRUE),
    FalsePos    = as.integer((fp %||% 0) > 0),
    Runtime     = if (!is.na(wall_col)) dl[[wall_col]] / 60 else NA_real_)

  ## MARE per dataset = mean |rel err| over the four covariate betas (x100)
  rl <- .sc_cell(.sc_read(src, src$rse_long), src, sample_N, struct_keep)
  mare <- rl |>
    dplyr::filter(.data$parameter %in% .SC_BETAS) |>
    dplyr::group_by(.data$scenario, .data$dataset_id) |>
    dplyr::summarise(Accuracy = mean(.data$abs_rel, na.rm = TRUE) * 100,
                     .groups = "drop")

  base |> dplyr::left_join(mare, by = c("scenario", "dataset_id"))
}

## ---- one comparison: method vs reference, one metric, one structure --------

.st_one <- function(method, structure, metric, sample_N = 80,
                    reference = "PsN-SCM") {
  na_row <- tibble::tibble(
    method = method, structure = structure, metric = metric,
    estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    effect_scale = NA_character_, p_raw = NA_real_, test = NA_character_,
    n_pairs = 0L, tested = FALSE, better = NA)

  if (identical(method, reference)) return(na_row)          # reference row
  if (.st_not_tested(method, structure, metric)) return(na_row)

  m <- .st_dataset_metrics(method,    structure, sample_N)
  r <- .st_dataset_metrics(reference, structure, sample_N)
  key <- c("scenario", "dataset_id")

  paired <- dplyr::inner_join(
    dplyr::transmute(m, scenario, dataset_id, y_m = .data[[metric]]),
    dplyr::transmute(r, scenario, dataset_id, y_r = .data[[metric]]),
    by = key) |>
    dplyr::filter(is.finite(.data$y_m), is.finite(.data$y_r))

  if (nrow(paired) < 3L) return(na_row)

  binary <- metric %in% .ST_BINARY
  better <- if (binary) {
    (mean(paired$y_m) - mean(paired$y_r)) *
      (.SC_META$dir[.SC_META$metric == metric]) > 0
  } else {
    (stats::median(paired$y_m) - stats::median(paired$y_r)) *
      (.SC_META$dir[.SC_META$metric == metric]) > 0
  }

  out <- na_row
  out$n_pairs <- nrow(paired)
  out$better  <- better

  ## ---- PRIMARY: mixed-effects model (if lme4 available) ------------------
  fit_mixed <- function() {
    dd <- dplyr::bind_rows(
      dplyr::transmute(paired, scenario, dataset_id,
                       method = "cmp", y = y_m),
      dplyr::transmute(paired, scenario, dataset_id,
                       method = "ref", y = y_r)) |>
      dplyr::mutate(method = factor(method, levels = c("ref", "cmp")),
                    grp = interaction(scenario, dataset_id, drop = TRUE))
    form_re <- y ~ method + (1 | scenario) + (1 | grp)
    if (binary) {
      fm <- lme4::glmer(form_re, data = dd, family = stats::binomial,
                        control = lme4::glmerControl(optimizer = "bobyqa"))
      sc <- "log-odds (OR)"
    } else if (identical(metric, "Runtime")) {
      dd$y <- log(pmax(dd$y, .Machine$double.eps))
      fm <- lme4::lmer(y ~ method + (1 | scenario) + (1 | grp), data = dd)
      sc <- "log-fold (x)"
    } else {                                        # Accuracy (MARE)
      fm <- lme4::lmer(form_re, data = dd)
      sc <- "mean diff (%)"
    }
    co <- summary(fm)$coefficients
    est <- co["methodcmp", "Estimate"]
    se  <- co["methodcmp", "Std. Error"]
    z   <- est / se
    p   <- 2 * stats::pnorm(-abs(z))
    list(estimate = est, ci_low = est - 1.96 * se, ci_high = est + 1.96 * se,
         effect_scale = sc, p = p, test = "mixed-effects")
  }

  ## ---- FALLBACK / sensitivity: exact paired tests ------------------------
  fit_paired <- function() {
    if (binary) {
      b <- sum(paired$y_m == 1 & paired$y_r == 0)   # method better
      c <- sum(paired$y_m == 0 & paired$y_r == 1)   # method worse
      disc <- b + c
      if (disc == 0L)
        return(list(estimate = 0, ci_low = NA, ci_high = NA,
                    effect_scale = "risk diff", p = NA_real_,
                    test = "McNemar (no discordant)"))
      p <- stats::binom.test(b, disc, 0.5)$p.value    # exact McNemar
      list(estimate = (b - c) / nrow(paired),
           ci_low = NA_real_, ci_high = NA_real_,
           effect_scale = "risk diff", p = p, test = "McNemar (exact)")
    } else {
      wt <- suppressWarnings(
        stats::wilcox.test(paired$y_m, paired$y_r, paired = TRUE,
                           conf.int = TRUE, exact = FALSE))
      list(estimate = unname(wt$estimate),
           ci_low = wt$conf.int[1], ci_high = wt$conf.int[2],
           effect_scale = "Hodges-Lehmann", p = wt$p.value,
           test = "Wilcoxon (paired)")
    }
  }

  res <- tryCatch(
    if (.ST_HAVE_LME4) fit_mixed() else fit_paired(),
    error = function(e) fit_paired())

  ## Separation guard for ceiling binary metrics (e.g. convergence at 100%):
  ## complete/quasi-separation inflates the glmer Wald SE so the p-value drifts
  ## to ~1 despite an obvious difference. Detect the degenerate fit (implausible
  ## log-odds or non-significant p alongside a clear raw gap) and fall back to
  ## the exact McNemar test, which stays valid on the discordant pairs.
  if (binary && identical(res$test, "mixed-effects")) {
    degenerate <- !is.finite(res$estimate) || abs(res$estimate) > 10 ||
      !is.finite(res$p) ||
      (res$p > 0.05 && abs(mean(paired$y_m) - mean(paired$y_r)) > 0.05)
    if (isTRUE(degenerate)) {
      mc <- fit_paired()
      mc$test <- paste0(mc$test, " [separation fallback]")
      res <- mc
    }
  }

  out$estimate     <- res$estimate
  out$ci_low       <- res$ci_low
  out$ci_high      <- res$ci_high
  out$effect_scale <- res$effect_scale
  out$p_raw        <- res$p
  out$test         <- res$test
  out$tested       <- is.finite(res$p)
  out
}

## ---- driver: full stats table over method x structure x metric -------------

build_scorecard_stats <- function(methods = .SC_LEVELS,
                                   structures = c("Analytic", "ODE"),
                                   metrics = .SC_METRIC_ORDER,
                                   reference = "PsN-SCM", sample_N = 80,
                                   adjust = "holm", save = FALSE,
                                   out_dir = "output/figures/scorecard") {
  grid <- expand.grid(method = methods, structure = structures,
                      metric = metrics, stringsAsFactors = FALSE)
  st <- dplyr::bind_rows(lapply(seq_len(nrow(grid)), function(i)
    .st_one(grid$method[i], grid$structure[i], grid$metric[i],
            sample_N, reference)))

  ## Holm across the actually-tested comparisons only
  st <- st |>
    dplyr::mutate(p_adj = ifelse(tested, p_raw, NA_real_))
  tested_idx <- which(st$tested)
  st$p_adj[tested_idx] <- stats::p.adjust(st$p_raw[tested_idx], method = adjust)

  st <- st |>
    dplyr::mutate(
      marker = .st_marker(p_adj, tested, method, reference,
                          structure, metric),
      method = factor(method, levels = methods),
      structure = factor(structure, levels = structures),
      metric = factor(metric, levels = metrics))

  if (isTRUE(save)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    csv <- file.path(out_dir, "scorecard_stats.csv")
    readr::write_csv(st, csv)
    message("saved: ", csv)
  }
  st
}

## significance glyph from an adjusted p-value
.st_marker <- function(p_adj, tested, method, reference, structure, metric) {
  is_ref <- as.character(method) == reference
  is_nt  <- mapply(.st_not_tested, method, structure, metric)
  out <- rep("", length(p_adj))
  out[is_ref]          <- ""            # reference is the baseline, no test
  out[!is_ref & is_nt] <- "\u2020"      # dagger = descriptive, not tested
  ok <- !is_ref & !is_nt & tested
  out[ok] <- ifelse(p_adj[ok] < 0.001, "***",
              ifelse(p_adj[ok] < 0.01,  "**",
              ifelse(p_adj[ok] < 0.05,  "*", "ns")))
  ## tested-but-degenerate (e.g. no discordant pairs) -> ns
  out[!is_ref & !is_nt & !tested] <- "ns"
  out
}
