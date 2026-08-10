## fig_runtime_compare.R
## Cross-method WALL-TIME comparison (nlmixr2-SCM / PsN-SCM / nlmixr2-VAE).
##
##   fig_runtime_bymethod()  per-method impact of parameterization x sample size
##                           x scenario: x = scenario, median line + IQR ribbon,
##                           facet = structure (cols) x sample_N (rows), ONE
##                           method (colour fixed).  Run per method to see how
##                           each factor drives its runtime.
##
##   fig_runtime_threeway()  three methods OVERLAID at a fixed sample size
##                           (default N = 80): x = scenario, median line + IQR
##                           ribbon, colour = METHOD, facet = structure
##                           (Analytic | ODE).  ODE is where PsN diverges.
##
## COMPARABILITY CAVEAT (annotated on every figure):
##   nlmixr2-SCM / PsN-SCM wall = wall_total_sec  (FULL pipeline: base fit +
##     stepwise covariate search + tight-tol refit).
##   nlmixr2-VAE       wall = fit_runtime_sec (SINGLE fit; the VAE delivers a
##     covariate model in one pass, no stepwise search).
##   => a fair "time to deliver a covariate model" comparison, NOT identical
##      workloads.  Set wall = "wall_base_sec" for a base-fit-only view.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

## ---- method sources --------------------------------------------------------

.RT_SOURCES <- list(
  "nlmixr2-SCM" = list(
    path      = "output/scm_bench_rescue_winner_aggregated/scm_diag_long.csv",
    wall      = "wall_total_sec",
    estimator = "focei",
    outer_opt = "bobyqa"
  ),
  "PsN-SCM" = list(
    path      = "output/psn_scm_combined_aggregated/scm_diag_long.csv",
    wall      = "wall_total_sec",
    estimator = "nonmem_scm",
    outer_opt = "focei"
  ),
  "nlmixr2-VAE" = list(
    path      = "output/vae_covsel_full0729_est702_aggregated/vae_diag_long.csv",
    wall      = "fit_runtime_sec",
    estimator = NULL,
    outer_opt = NULL
  )
)

.RT_LEVELS  <- c("nlmixr2-SCM", "PsN-SCM", "nlmixr2-VAE")
.PAL_METHOD_R <- c(
  "nlmixr2-SCM"  = "#0F8B8D",   # teal
  "PsN-SCM"      = "#2B2B2B",   # near-black
  "nlmixr2-VAE"  = "#A64AC9"    # purple
)

.RT_CAVEAT <- paste(
  "SCM wall = full pipeline (base + stepwise search + refit);",
  "VAE wall = single fit (no search).  Time to a covariate model, not identical workloads.")

## ---- loader ----------------------------------------------------------------
## returns per-fit rows harmonised to (method, structure, sample_N, scenario,
## wall_min); structure folded to Analytic / ODE so linCmt & advan4 align.

.rt_load <- function(sources, methods, sample_N = NULL, structure = NULL) {
  rows <- lapply(methods, function(m) {
    s <- sources[[m]]
    d <- suppressMessages(readr::read_csv(s$path, show_col_types = FALSE))
    if (!is.null(s$estimator) && "estimator" %in% names(d))
      d <- dplyr::filter(d, .data$estimator == s$estimator)
    if (!is.null(s$outer_opt) && "outer_opt" %in% names(d))
      d <- dplyr::filter(d, .data$outer_opt == s$outer_opt)
    wall_col <- s$wall %||%
      (if ("wall_total_sec" %in% names(d)) "wall_total_sec" else "fit_runtime_sec")
    d |>
      dplyr::transmute(
        method    = m,
        structure = structure,               # raw CSV structure (linCmt/advan4/ode)
        sample_N  = sample_N,
        scenario  = as.integer(scenario),
        wall_min  = .data[[wall_col]] / 60
      )
  })

  dat <- dplyr::bind_rows(rows) |>
    dplyr::filter(is.finite(wall_min))

  ## filter on RAW structure/sample_N first, THEN fold to Analytic/ODE
  if (!is.null(sample_N)) dat <- dplyr::filter(dat, sample_N %in% !!sample_N)
  if (!is.null(structure)) {
    keep <- if (identical(structure, "linCmt")) c("linCmt", "advan4") else structure
    dat <- dplyr::filter(dat, structure %in% keep)
  }

  dat |>
    dplyr::mutate(
      structure = dplyr::case_when(
        structure %in% c("linCmt", "advan4") ~ "Analytic",
        structure == "ode"                    ~ "ODE",
        TRUE                                  ~ structure),
      structure = factor(structure, levels = c("Analytic", "ODE")),
      method    = factor(method, levels = methods))
}

## per (method, structure, sample_N, scenario) median + IQR band
.rt_summ <- function(dat) {
  dat |>
    dplyr::group_by(method, structure, sample_N, scenario) |>
    dplyr::summarise(
      med = stats::median(wall_min),
      q25 = stats::quantile(wall_min, 0.25, names = FALSE),
      q75 = stats::quantile(wall_min, 0.75, names = FALSE),
      .groups = "drop")
}

.rt_theme <- function(base_size = 15) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "top",
      strip.text       = element_text(size = base_size),
      strip.background = element_rect(fill = "grey92", colour = NA),
      plot.caption     = element_text(size = base_size - 6, hjust = 0,
                                      colour = "grey35"))
}

## ---- Figure A : per-method impact (structure x sample_N grid) --------------

fig_runtime_bymethod <- function(method    = "PsN-SCM",
                                 sources   = .RT_SOURCES,
                                 log_y     = TRUE,
                                 save      = FALSE,
                                 out_dir   = "output/figures/runtime_compare",
                                 base_size = 15) {
  ## load ALL N and BOTH structures for this one method
  raw <- lapply(c(40L, 80L, 300L), function(n) {
    lapply(c("linCmt", "ode"), function(st)
      .rt_load(sources, method, sample_N = n, structure = st)) |>
      dplyr::bind_rows()
  }) |> dplyr::bind_rows()
  summ <- .rt_summ(raw) |>
    dplyr::mutate(N_lab = factor(paste0("N = ", sample_N),
                                 levels = paste0("N = ", c(40, 80, 300))))

  col <- .PAL_METHOD_R[[method]]
  p <- ggplot(summ, aes(scenario, med)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, fill = col) +
    geom_line(colour = col, linewidth = 0.7) +
    geom_point(colour = col, size = 1.4) +
    facet_grid(N_lab ~ structure, scales = "free_y") +
    scale_x_continuous(breaks = seq(2, 16, 2)) +
    labs(title = sprintf("%s wall time by scenario", method),
         x = "Simulation scenario",
         y = "Median wall time per fit (min), band = IQR",
         caption = .RT_CAVEAT) +
    .rt_theme(base_size)
  if (log_y) p <- p + scale_y_continuous(trans = "log2")

  if (isTRUE(save)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    tag  <- gsub("[^A-Za-z0-9]+", "", method)
    stub <- file.path(out_dir, sprintf("fig_runtime_bymethod_%s", tag))
    ggsave(paste0(stub, ".png"), p, width = 10, height = 7, dpi = 200)
    ggsave(paste0(stub, ".pdf"), p, width = 10, height = 7)
    message("saved: ", stub, ".{png,pdf}")
  }
  p
}

## ---- Figure B : three-way overlay at fixed N (default 80) ------------------

fig_runtime_threeway <- function(sources   = .RT_SOURCES,
                                 methods   = .RT_LEVELS,
                                 sample_N  = 80,
                                 log_y     = TRUE,
                                 save      = FALSE,
                                 out_dir   = "output/figures/runtime_compare",
                                 base_size = 15) {
  raw <- lapply(c("linCmt", "ode"), function(st)
    .rt_load(sources, methods, sample_N = sample_N, structure = st)) |>
    dplyr::bind_rows()
  summ <- .rt_summ(raw)

  p <- ggplot(summ, aes(scenario, med, colour = method, fill = method,
                        group = method)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    facet_wrap(~ structure, nrow = 1) +
    scale_colour_manual(values = .PAL_METHOD_R[methods], name = "Method") +
    scale_fill_manual(values = .PAL_METHOD_R[methods], guide = "none") +
    scale_x_continuous(breaks = 1:16) +
    labs(title = sprintf("Wall time by scenario (N = %d)", sample_N),
         x = "Simulation scenario",
         y = "Median wall time per fit (min), band = IQR",
         caption = .RT_CAVEAT) +
    .rt_theme(base_size)
  if (log_y) p <- p + scale_y_continuous(trans = "log2")

  if (isTRUE(save)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    stub <- file.path(out_dir,
      sprintf("fig_runtime_threeway_N%d", sample_N))
    ggsave(paste0(stub, ".png"), p, width = 12, height = 5, dpi = 200)
    ggsave(paste0(stub, ".pdf"), p, width = 12, height = 5)
    message("saved: ", stub, ".{png,pdf}")
  }
  p
}

## ---- Figure C : three-way, sample_N (rows) x structure (cols) grid ---------
## all sample sizes at once; three methods overlaid (colour) within each panel.

fig_runtime_grid <- function(sources   = .RT_SOURCES,
                             methods   = .RT_LEVELS,
                             sample_N  = c(40L, 80L, 300L),
                             log_y     = TRUE,
                             save      = FALSE,
                             out_dir   = "output/figures/runtime_compare",
                             base_size = 15) {
  raw <- lapply(sample_N, function(n)
    lapply(c("linCmt", "ode"), function(st)
      .rt_load(sources, methods, sample_N = n, structure = st)) |>
      dplyr::bind_rows()) |>
    dplyr::bind_rows()
  summ <- .rt_summ(raw) |>
    dplyr::mutate(N_lab = factor(paste0("N = ", sample_N),
                                 levels = paste0("N = ", sort(unique(sample_N)))))

  p <- ggplot(summ, aes(scenario, med, colour = method, fill = method,
                        group = method)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.6) +
    facet_grid(N_lab ~ structure) +
    scale_colour_manual(values = .PAL_METHOD_R[methods], name = "Method") +
    scale_fill_manual(values = .PAL_METHOD_R[methods], guide = "none") +
    scale_x_continuous(breaks = seq(2, 16, 2)) +
    labs(title = "Wall time by scenario, sample size and model type",
         x = "Simulation scenario",
         y = "Median wall time per fit (min), band = IQR",
         caption = .RT_CAVEAT) +
    .rt_theme(base_size)
  if (log_y) p <- p + scale_y_continuous(trans = "log2")

  if (isTRUE(save)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    stub <- file.path(out_dir,
      paste0("fig_runtime_grid_allN_", if (log_y) "log2" else "linear"))
    ggsave(paste0(stub, ".png"), p, width = 12, height = 8, dpi = 200)
    ggsave(paste0(stub, ".pdf"), p, width = 12, height = 8)
    message("saved: ", stub, ".{png,pdf}")
  }
  p
}
