## scorecard.R
## Result SYNTHESIS across the six evidence streams (convergence, covariance
## step, power, false-positive, estimation accuracy, runtime) into an objective,
## glanceable scorecard.
##
## Frame (stated explicitly for honesty):
##   * reference cell   : N = 80, per structure (Analytic + ODE)
##   * scenario summary : MEDIAN over the 16 scenarios (robust, one number)
##   * reference method : PsN-SCM (regulatory-standard baseline)
##
## Reduction rule (identical for every method x structure):
##   Convergence  median % minimised            (higher better)
##   CovStep      median % covariance-step OK    (higher better)
##   Power        median exact-match Power %     (higher better)
##   FalsePos     median distractor FP rate %    (lower  better)
##   Accuracy     median MARE % of the 4 covariate betas, True-selection regime
##                                                (lower  better)
##   Runtime      median wall-min per fit        (lower  better)
##
## OUTPUTS
##   scorecard_table.{csv,html}       Layer A: six raw scalars, natural units
##   fig_scorecard_heat.{png,pdf}     Layer B: Delta-vs-PsN diverging heatmap
##   fig_scorecard_radar.{png,pdf}    petal plot, six desirability spokes
##
## Caveats encoded (not hidden):
##   * VAE runtime = single fit (no stepwise search)  -> dagger footnote
##   * PsN-ODE covariance/CN partially n/a            -> grey n/a cell
##   * Analytic = linCmt (nlmixr2) vs ADVAN4 (PsN)    -> caption

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

if (!exists(".fix_covsel_denom")) source("script/viz/fig_covsel_heatmap.R")

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

## ---- per-method source files + fixed estimator cell ------------------------

.SC_SOURCES <- list(
  "nlmixr2-SCM" = list(
    dir = "output/scm_focei_bobyqa_703est_0805_aggregated",
    diag_rates = "scm_diag_rates.csv", power = "scm_power.csv",
    by_covar = "scm_covsel_by_covar.csv", estim_cond = "scm_estim_cond.csv",
    diag_long = "scm_diag_long.csv", rse_long = "scm_rse_long.csv",
    wall = "wall_total_sec",
    estimator = "focei", outer_opt = "bobyqa"),
  "PsN-SCM" = list(
    dir = "output/psn_scm_combined_aggregated",
    diag_rates = "scm_diag_rates.csv", power = "scm_power.csv",
    by_covar = "scm_covsel_by_covar.csv", estim_cond = "scm_estim_cond.csv",
    diag_long = "scm_diag_long.csv", rse_long = "scm_rse_long.csv",
    wall = "wall_total_sec",
    estimator = "nonmem_scm", outer_opt = "focei"),
  "nlmixr2-VAE" = list(
    dir = "output/vae_covsel_full0807_est703_aggregated",
    diag_rates = "vae_diag_rates.csv", power = "vae_power.csv",
    by_covar = "vae_covsel_by_covar.csv", estim_cond = "vae_estim_cond.csv",
    diag_long = "vae_diag_long.csv", rse_long = "vae_rse_long.csv",
    wall = "fit_runtime_sec",
    estimator = NULL, outer_opt = NULL)
)

.SC_LEVELS  <- c("nlmixr2-SCM", "PsN-SCM", "nlmixr2-VAE")
.SC_BETAS   <- c("CLBW", "CLcrCL", "VcBW", "VcSEX")  # covariate effects

## metric metadata: label, direction (+1 higher better / -1 lower better), unit
.SC_META <- tibble::tribble(
  ~metric,       ~label,          ~dir, ~unit,
  "Convergence", "Convergence",    1,   "%",
  "CovStep",     "Cov. step",      1,   "%",
  "Power",       "Power",          1,   "%",
  "FalsePos",    "FP",            -1,   "%",
  "Accuracy",    "MARE",          -1,   "%",
  "Runtime",     "Runtime",       -1,   "min")
.SC_METRIC_ORDER <- .SC_META$metric

## ---- helpers: filter one method's file to its estimator cell + N + structure

.sc_read <- function(src, file) {
  p <- file.path(src$dir, file)
  if (!file.exists(p)) stop("scorecard: file not found: ", p)
  suppressMessages(readr::read_csv(p, show_col_types = FALSE))
}
.sc_cell <- function(d, src, sample_N, struct_keep) {
  if (!is.null(src$estimator) && "estimator" %in% names(d))
    d <- dplyr::filter(d, .data$estimator == src$estimator)
  if (!is.null(src$outer_opt) && "outer_opt" %in% names(d))
    d <- dplyr::filter(d, .data$outer_opt == src$outer_opt)
  d <- dplyr::filter(d, .data$sample_N == !!sample_N,
                        .data$structure %in% struct_keep)
  d
}

## ---- reduce ONE method x structure to the six scalars ----------------------

.sc_one <- function(method, structure, sample_N = 80) {
  src <- .SC_SOURCES[[method]]
  struct_keep <- if (identical(structure, "Analytic"))
                   c("linCmt", "advan4") else "ode"

  ## convergence + covariance step (already 0-100), median over scenarios
  dr <- .sc_cell(.sc_read(src, src$diag_rates), src, sample_N, struct_keep)
  conv    <- stats::median(dr$Converged_pct, na.rm = TRUE)
  covstep <- if ("CovStep_pct" %in% names(dr))
               stats::median(dr$CovStep_pct, na.rm = TRUE) else NA_real_

  ## power (fraction 0-1) -> %, median over scenarios
  pw <- .sc_cell(.sc_read(src, src$power), src, sample_N, struct_keep)
  power <- stats::median(pw$Power, na.rm = TRUE) * 100

  ## false-positive rate: distractor cells (is_true == FALSE), n_FP/n_datasets
  bc <- .sc_cell(.sc_read(src, src$by_covar), src, sample_N, struct_keep) |>
    .fix_covsel_denom()
  fp_rows <- dplyr::filter(bc, !is_true, n_datasets > 0)
  falsepos <- if (nrow(fp_rows))
    stats::median(fp_rows$n_FP / fp_rows$n_datasets, na.rm = TRUE) * 100 else NA_real_

  ## accuracy: median MARE of the 4 covariate betas, True-selection (cond) regime
  ec <- .sc_cell(.sc_read(src, src$estim_cond), src, sample_N, struct_keep)
  betas <- dplyr::filter(ec, parameter %in% .SC_BETAS)
  accuracy <- stats::median(betas$MARE_pct, na.rm = TRUE)

  ## runtime: median wall-min per fit
  dl <- .sc_cell(.sc_read(src, src$diag_long), src, sample_N, struct_keep)
  wall_col <- src$wall %||%
    (if ("wall_total_sec" %in% names(dl)) "wall_total_sec" else "fit_runtime_sec")
  runtime <- stats::median(dl[[wall_col]] / 60, na.rm = TRUE)

  tibble::tibble(
    method = method, structure = structure,
    Convergence = conv, CovStep = covstep, Power = power,
    FalsePos = falsepos, Accuracy = accuracy, Runtime = runtime)
}

## ---- build the full scorecard (Layer A raw table) --------------------------

build_scorecard <- function(sources = .SC_SOURCES, methods = .SC_LEVELS,
                            structures = c("Analytic", "ODE"), sample_N = 80,
                            save = FALSE,
                            out_dir = "output/figures/scorecard") {
  grid <- expand.grid(method = methods, structure = structures,
                      stringsAsFactors = FALSE)
  tab <- purrr_safe_map(grid, sample_N)
  tab <- tab |>
    dplyr::mutate(method = factor(method, levels = methods),
                  structure = factor(structure, levels = structures)) |>
    dplyr::arrange(structure, method)

  if (isTRUE(save)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    csv <- file.path(out_dir, "scorecard_table.csv")
    readr::write_csv(tab, csv)
    message("saved: ", csv)
    .sc_write_html(tab, file.path(out_dir, "scorecard_table.html"), sample_N)
  }
  tab
}

## small internal: map over the method x structure grid (no external deps)
purrr_safe_map <- function(grid, sample_N) {
  dplyr::bind_rows(lapply(seq_len(nrow(grid)), function(i)
    .sc_one(grid$method[i], grid$structure[i], sample_N)))
}

## Layer A -> minimal styled HTML table
.sc_write_html <- function(tab, path, sample_N) {
  fmt <- function(x, u) ifelse(is.na(x), "n/a",
    paste0(formatC(x, format = "f", digits = 1), u))
  units <- setNames(.SC_META$unit, .SC_META$metric)
  body <- apply(tab, 1, function(r) {
    cells <- vapply(.SC_METRIC_ORDER, function(m)
      sprintf("<td style='text-align:right'>%s</td>",
              fmt(as.numeric(r[[m]]), units[[m]])), character(1))
    sprintf("<tr><td>%s</td><td>%s</td>%s</tr>",
            r[["method"]], r[["structure"]], paste(cells, collapse = ""))
  })
  hdr <- paste0("<th>", c("Method", "Structure", .SC_META$label),
                "</th>", collapse = "")
  html <- sprintf(paste0(
    "<html><head><meta charset='utf-8'><style>",
    "table{border-collapse:collapse;font-family:sans-serif;font-size:14px}",
    "th,td{border:1px solid #ccc;padding:6px 10px}",
    "th{background:#f0f0f0}</style></head><body>",
    "<h3>Scorecard (N = %d, median over 16 scenarios)</h3>",
    "<table><thead><tr>%s</tr></thead><tbody>%s</tbody></table>",
    "<p style='font-size:12px;color:#555'>Analytic = linCmt (nlmixr2) / ",
    "ADVAN4 (NONMEM).  VAE runtime = single fit (no stepwise search).  ",
    "n/a = metric unavailable for that cell.</p></body></html>"),
    sample_N, hdr, paste(body, collapse = ""))
  writeLines(html, path)
  message("saved: ", path)
}

## ---- Layer B : Delta-vs-PsN diverging heatmap ------------------------------

fig_scorecard_heat <- function(tab = NULL, stats = NULL, reference = "PsN-SCM",
                               sample_N = 80, show_sig = TRUE,
                               practical_delta = 0.10, absolute_floor = 2,
                               save = FALSE,
                               out_dir = "output/figures/scorecard",
                               base_size = 14) {
  if (is.null(tab)) tab <- build_scorecard(sample_N = sample_N)

  ## significance markers from the paired mixed-effects / McNemar table
  if (isTRUE(show_sig)) {
    if (is.null(stats) && exists("build_scorecard_stats")) {
      stats <- tryCatch(
        build_scorecard_stats(reference = reference, sample_N = sample_N),
        error = function(e) { message("scorecard stats skipped: ", e$message)
                              NULL })
    }
  } else {
    stats <- NULL
  }

  long <- tab |>
    tidyr::pivot_longer(dplyr::all_of(.SC_METRIC_ORDER),
                        names_to = "metric", values_to = "value") |>
    dplyr::left_join(.SC_META, by = "metric")

  ## direction-corrected Delta vs the reference method, WITHIN each structure
  ref <- long |>
    dplyr::filter(method == reference) |>
    dplyr::select(structure, metric, ref_value = value)

  long <- long |>
    dplyr::left_join(ref, by = c("structure", "metric")) |>
    dplyr::mutate(
      ## improvement = dir * (method - ref): positive = BETTER than PsN
      improve = dir * (value - ref_value),
      metric  = factor(metric, levels = .SC_METRIC_ORDER),
      metric_lab = factor(label, levels = .SC_META$label),
      method  = factor(method, levels = .SC_LEVELS),
      row_lab = factor(paste0(method, "  |  ", structure)),
      cell_txt = ifelse(is.na(value), "n/a",
                        paste0(formatC(value, format = "f", digits = 1), unit)),
      ## colour = RELATIVE change vs the reference (magnitude-faithful, unit-free)
      ## so a 0.5pp MARE gap stays pale while a 20pp power gap saturates. A tiny
      ## but statistically-significant difference no longer looks catastrophic.
      ## Zero-reference cells (e.g. PsN cov-step = 0%) can't form a ratio, so a
      ## non-zero method value is treated as a full-scale improvement.
      improve_rel = dplyr::case_when(
        is.na(improve)                 ~ NA_real_,
        abs(ref_value) < .Machine$double.eps & improve == 0 ~ 0,
        abs(ref_value) < .Machine$double.eps ~ sign(improve) * 1,
        TRUE                           ~ improve / abs(ref_value)),
      ## clip to +/-100% so extreme ratios don't dominate the palette
      improve_s = pmax(pmin(improve_rel, 1), -1))

  lim <- 1

  ## attach significance markers (blank for the reference method's own rows)
  if (!is.null(stats)) {
    mk <- stats |>
      dplyr::mutate(metric = as.character(metric),
                    method = as.character(method),
                    structure = as.character(structure)) |>
      dplyr::select(method, structure, metric, marker)
    long <- long |>
      dplyr::mutate(metric = as.character(metric),
                    method = as.character(method),
                    structure = as.character(structure)) |>
      dplyr::left_join(mk, by = c("method", "structure", "metric")) |>
      dplyr::mutate(
        metric  = factor(metric, levels = .SC_METRIC_ORDER),
        marker  = ifelse(is.na(marker), "", marker),
        ## PRACTICAL-EQUIVALENCE band: a difference that is statistically
        ## significant (starred) but too small to matter is reliable-but-
        ## negligible, not "better/worse", so it is re-labelled "~".
        ## A star never implies importance. Negligible is judged by EITHER:
        ##   * relative:  within +/-practical_delta of the reference
        ##                (right yardstick for large baselines, e.g. MARE ~12%)
        ##   * absolute:  within +/-absolute_floor natural units for %-metrics
        ##                (guards near-zero baselines where a 1pp FP change is a
        ##                 huge relative % but practically trivial)
        negligible = (is.finite(improve_rel) & abs(improve_rel) < practical_delta) |
          (unit == "%" & is.finite(improve) & abs(improve) < absolute_floor),
        marker = ifelse(grepl("\\*", marker) & negligible, "\u2248", marker))
  } else {
    long$marker <- ""
  }

  ## faceting by structure inserts the blank margin between Analytic and ODE;
  ## within a facet the row label is just the method
  long <- long |>
    dplyr::mutate(
      method    = factor(method, levels = rev(.SC_LEVELS)),
      structure = factor(structure, levels = c("Analytic", "ODE")))

  p <- ggplot(long, aes(metric_lab, method, fill = improve_s)) +
    geom_tile(colour = "grey85", linewidth = 0.4) +
    geom_text(aes(label = cell_txt), size = 4, fontface = "bold",
              nudge_y = 0.14,
              colour = ifelse(!is.na(long$improve_s) &
                              abs(long$improve_s) > 0.6, "white", "grey15")) +
    geom_text(aes(label = marker), size = 4.6, fontface = "bold",
              nudge_y = -0.20,
              colour = ifelse(!is.na(long$improve_s) &
                              abs(long$improve_s) > 0.6, "white", "grey25")) +
    facet_grid(structure ~ ., scales = "free_y", space = "free_y",
               switch = "y") +
    scale_fill_gradient2(low = "#B30000", mid = "white", high = "#08519C",
      midpoint = 0, limits = c(-1, 1), na.value = "grey88",
      breaks = c(-1, 0, 1),
      labels = c("\u2264 -100%", paste0("= ", reference), "\u2265 +100%"),
      name = paste0("Change vs ", reference, "\n(blue better, red worse)")) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL,
         caption = paste0(
           "N = ", sample_N, ", median over 16 scenarios.  ",
           "NUMBER = raw value (natural units); COLOUR = RELATIVE change vs ",
           reference, " (\u00b1100% clipped), so small gaps stay pale and large ",
           "gaps saturate.  ",
           "MARKER = paired significance vs ", reference,
           " (mixed-effects / exact McNemar, Holm-adjusted): ",
           "*** p<0.001, ** p<0.01, * p<0.05, ns not sig., \u2020 not tested.  ",
           "\u2248 = statistically significant but within \u00b1",
           round(100 * practical_delta), "% of ", reference,
           " or \u00b1", absolute_floor, "pp (reliable but practically ",
           "negligible).  ",
           "MARE cell = median over scenarios of the per-scenario median ",
           "absolute relative error of the 4 covariate \u03b2s (robust display); ",
           "the significance marker instead tests the per-replicate MEAN of ",
           "the 4 \u03b2s in the mixed model - same accuracy construct, ",
           "median (display) vs mean (test) collapse over the \u03b2s.  ",
           "Runtime = wall-clock to completed covariate selection ",
           "(full stepwise pipeline for SCM; single joint fit for VAE).  ",
           "Analytic = linCmt / ADVAN4.")) +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.spacing.y = grid::unit(10, "pt"),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold",
                                       size = base_size),
      legend.position = "right",
      axis.text.x.top = element_text(size = base_size, face = "bold"),
      axis.text.y     = element_text(size = base_size),
      plot.caption    = element_text(size = base_size - 4, hjust = 0,
                                     colour = "grey35"))

  if (isTRUE(save)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    stub <- file.path(out_dir, "fig_scorecard_heat")
    ggsave(paste0(stub, ".png"), p, width = 11, height = 5, dpi = 200)
    ggsave(paste0(stub, ".pdf"), p, width = 11, height = 5)
    message("saved: ", stub, ".{png,pdf}")
  }
  p
}

## ---- radar / petal : six desirability spokes per method --------------------
## desirability in [0,1] via direction-corrected min-max across the 6 rows for
## each metric (best cell = 1, worst = 0).  NO single weighted composite.

fig_scorecard_radar <- function(tab = NULL, reference = "PsN-SCM",
                                sample_N = 80, save = FALSE,
                                out_dir = "output/figures/scorecard",
                                base_size = 13) {
  if (is.null(tab)) tab <- build_scorecard(sample_N = sample_N)

  long <- tab |>
    tidyr::pivot_longer(dplyr::all_of(.SC_METRIC_ORDER),
                        names_to = "metric", values_to = "value") |>
    dplyr::left_join(.SC_META, by = "metric")

  ## reference value per structure x metric (same frame as the heatmap)
  ref <- long |>
    dplyr::filter(method == reference) |>
    dplyr::select(structure, metric, ref_value = value)

  ## RELATIVE scale (not absolute): performance vs the reference method, scaled
  ## within each metric column, then mapped to [0, 1] with reference = 0.5.
  ## outer ring = clearly better than reference, centre = clearly worse.
  long <- long |>
    dplyr::left_join(ref, by = c("structure", "metric")) |>
    dplyr::group_by(metric) |>
    dplyr::mutate(
      improve = dir * (value - ref_value),
      col_lim = max(abs(improve), na.rm = TRUE),
      desir   = ifelse(is.na(value), NA_real_,
                       0.5 + 0.5 * ifelse(col_lim > 0, improve / col_lim, 0))) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      metric_lab = factor(label, levels = .SC_META$label),
      method  = factor(method, levels = .SC_LEVELS),
      structure = factor(structure, levels = c("Analytic", "ODE")))

  pal_struct <- c("Analytic" = "#1F77B4", "ODE" = "#D62728")

  p <- ggplot(long, aes(metric_lab, desir,
                        colour = structure, group = structure)) +
    geom_hline(yintercept = 0.5, colour = "grey70", linetype = 2,
               linewidth = 0.3) +
    geom_polygon(aes(fill = structure), alpha = 0.12, linewidth = 0.8) +
    geom_point(size = 2) +
    coord_radar() +
    facet_wrap(~ method, nrow = 1) +
    scale_colour_manual(values = pal_struct, name = "Model type") +
    scale_fill_manual(values = pal_struct, guide = "none") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1),
                       labels = c("worse", reference, "better")) +
    labs(x = NULL, y = NULL,
         caption = paste0(
           "RELATIVE scale: each spoke is performance vs ", reference,
           " for that metric (dashed ring = ", reference,
           "; outside = better, inside = worse), scaled within each metric.  ",
           "N = ", sample_N, ", median over scenarios.  No weighted composite.")) +
    theme_minimal(base_size = base_size) +
    theme(legend.position = "top",
          strip.text = element_text(face = "bold", size = base_size + 1),
          axis.text.x = element_text(size = base_size - 3),
          plot.caption = element_text(size = base_size - 4, hjust = 0,
                                      colour = "grey35"))

  if (isTRUE(save)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    stub <- file.path(out_dir, "fig_scorecard_radar")
    ggsave(paste0(stub, ".png"), p, width = 12, height = 5, dpi = 200)
    ggsave(paste0(stub, ".pdf"), p, width = 12, height = 5)
    message("saved: ", stub, ".{png,pdf}")
  }
  p
}

## minimal radar coordinate system (closed polar) for ggplot
coord_radar <- function(theta = "x", start = 0, direction = 1) {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggproto("CoordRadar", CoordPolar, theta = theta, r = r,
          start = start, direction = sign(direction),
          is_linear = function(coord) TRUE)
}
