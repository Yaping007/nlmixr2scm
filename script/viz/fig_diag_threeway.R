## fig_diag_threeway.R
## Three-way method comparison of SCM RUN-QUALITY (diagnostic) rates.
##
##   Figure A  fig_diag_threeway_N80()   three methods OVERLAID as coloured
##                                       lines for ONE metric at a fixed N;
##                                       facets = Analytic | ODE.
##   Figure B  fig_diag_threeway_grid()  method columns x structure rows;
##                                       colour = sample size,
##                                       linetype = metric (Converged / CN).
##
## Diagnostic CSVs (columns already 0-100 percentages):
##   nlmixr2-SCM : output/scm_bench_rescue_winner_aggregated/scm_diag_rates.csv
##                 (filter estimator == "focei", outer_opt == "bobyqa")
##   PsN-SCM     : output/psn_scm_combined_aggregated/scm_diag_rates.csv
##                 (filter estimator == "nonmem_scm", outer_opt == "focei")
##   nlmixr2-VAE : output/vae_covsel_full0729_est702_aggregated/vae_diag_rates.csv
##                 (no estimator / outer_opt columns)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

## ---- configuration ---------------------------------------------------------

.TWD_SOURCES <- list(
  "nlmixr2-SCM" = list(
    path      = "output/scm_bench_rescue_winner_aggregated/scm_diag_rates.csv",
    estimator = "focei",
    outer_opt = "bobyqa"
  ),
  "PsN-SCM" = list(
    path      = "output/psn_scm_combined_aggregated/scm_diag_rates.csv",
    estimator = "nonmem_scm",
    outer_opt = "focei"
  ),
  "nlmixr2-VAE" = list(
    path      = "output/vae_covsel_full0729_est702_aggregated/vae_diag_rates.csv",
    estimator = NULL,
    outer_opt = NULL
  )
)

.TWD_LEVELS  <- c("nlmixr2-SCM", "PsN-SCM", "nlmixr2-VAE")
.PAL_METHOD_D <- c(
  "nlmixr2-SCM"  = "#0F8B8D",   # teal
  "PsN-SCM"      = "#2B2B2B",   # near-black
  "nlmixr2-VAE"  = "#A64AC9"    # purple
)

## sample-size palette (identical to fig_power.R: grey / orange / blue)
.TWD_PAL_N <- c("40" = "#7F7F7F", "80" = "#E8820C", "300" = "#1F77B4")

## metric -> CSV column, label, linetype
.TWD_METRIC_COL <- c(Converged     = "Converged_pct",
                     CNBelowCutoff = "CNBelowCutoff_pct")
.TWD_METRIC_LAB <- c(Converged     = "Converged",
                     CNBelowCutoff = "Cond# < 1000")
.TWD_METRIC_LT  <- c(Converged = "solid", CNBelowCutoff = "22")

## ---- loaders ---------------------------------------------------------------

.twd_prep <- function(sources, metrics) {
  cols <- unname(.TWD_METRIC_COL[metrics])
  rows <- lapply(names(sources), function(m) {
    s <- sources[[m]]
    d <- readr::read_csv(s$path, show_col_types = FALSE)

    if (!is.null(s$estimator) && "estimator" %in% names(d)) {
      d <- dplyr::filter(d, .data$estimator == s$estimator)
    }
    if (!is.null(s$outer_opt) && "outer_opt" %in% names(d)) {
      d <- dplyr::filter(d, .data$outer_opt == s$outer_opt)
    }

    have <- intersect(cols, names(d))
    d |>
      dplyr::select(dplyr::all_of(c("sample_N", "scenario",
                                    "structure", have))) |>
      tidyr::pivot_longer(dplyr::all_of(have),
                          names_to = "metric_col", values_to = "rate_pct") |>
      dplyr::mutate(method = m)
  })

  dplyr::bind_rows(rows) |>
    dplyr::filter(!is.na(.data$rate_pct)) |>
    dplyr::mutate(
      structure = dplyr::case_when(
        structure %in% c("linCmt", "advan4") ~ "Analytic",
        structure == "ode"                    ~ "ODE",
        TRUE                                  ~ structure
      ),
      structure = factor(structure, levels = c("Analytic", "ODE")),
      method    = factor(method, levels = .TWD_LEVELS),
      metric    = names(.TWD_METRIC_COL)[match(.data$metric_col,
                                               .TWD_METRIC_COL)]
    )
}

## ---- shared theme (mirrors fig_power_threeway.R) ---------------------------

.twd_theme <- function(base_size = 16) {
  theme_bw(base_size = base_size) +
    theme(
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text       = element_text(size = base_size),
      plot.title       = element_text(size = base_size - 2),
      axis.title       = element_text(size = base_size),
      axis.text.x      = element_text(size = base_size - 4),
      axis.text.y      = element_text(size = base_size - 2),
      legend.title     = element_text(size = base_size - 4),
      legend.text      = element_text(size = base_size - 4),
      legend.position  = "top",
      legend.box       = "horizontal",
      panel.grid.minor = element_blank()
    )
}

## ---- Figure A : overlaid methods, one metric, Analytic | ODE ---------------

fig_diag_threeway_N80 <- function(sources   = .TWD_SOURCES,
                                  sample_N  = 80,
                                  metric    = "Converged",
                                  save      = FALSE,
                                  out_dir   = "output/figures/threeway_N80",
                                  base_size = 16) {
  metric <- match.arg(metric, names(.TWD_METRIC_COL))
  dat <- .twd_prep(sources, metric) |>
    dplyr::filter(.data$sample_N == !!sample_N)

  p <- ggplot(dat, aes(scenario, rate_pct,
                       colour = method, group = method)) +
    geom_hline(yintercept = 80, linetype = "dashed", colour = "grey60") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    facet_wrap(~ structure, nrow = 1) +
    scale_colour_manual(values = .PAL_METHOD_D, name = "Method") +
    scale_x_continuous(breaks = 1:16) +
    scale_y_continuous(limits = c(0, 100),
                       labels = function(x) paste0(x, "%")) +
    labs(x = "Scenario", y = unname(.TWD_METRIC_LAB[metric])) +
    .twd_theme(base_size) +
    guides(colour = guide_legend(nrow = 1))

  if (save) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    stub <- file.path(out_dir,
      sprintf("fig_diag_threeway_N%d_%s", sample_N, metric))
    ggsave(paste0(stub, ".png"), p, width = 12, height = 4.5, dpi = 300)
    ggsave(paste0(stub, ".pdf"), p, width = 12, height = 4.5)
    message("saved: ", stub, ".{png,pdf}")
  }
  p
}

## ---- Figure B : 3x2 grid, ALL sample sizes + BOTH metrics ------------------
## method columns x structure rows; colour = sample size, linetype = metric.

fig_diag_threeway_grid <- function(sources   = .TWD_SOURCES,
                                   metrics   = c("Converged", "CNBelowCutoff"),
                                   save      = FALSE,
                                   out_dir   = "output/figures/threeway_N80",
                                   base_size = 16) {
  metrics <- intersect(metrics, names(.TWD_METRIC_COL))
  dat <- .twd_prep(sources, metrics) |>
    dplyr::mutate(
      N_lab  = factor(as.character(.data$sample_N),
                      levels = as.character(sort(unique(.data$sample_N)))),
      metric = factor(.data$metric, levels = names(.TWD_METRIC_COL))
    )

  n_lvls  <- levels(dat$N_lab)
  pal_N   <- .TWD_PAL_N[n_lvls]
  lt_lvls <- levels(droplevels(dat$metric))
  lt_vals <- .TWD_METRIC_LT[lt_lvls]

  p <- ggplot(dat, aes(scenario, rate_pct,
                       colour   = N_lab,
                       linetype = metric,
                       group    = interaction(N_lab, metric))) +
    geom_hline(yintercept = 80, linetype = "dashed", colour = "grey70") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.1) +
    facet_grid(structure ~ method) +
    scale_colour_manual(values = pal_N, name = "Sample size",
                        labels = function(x) paste0(x, " subj")) +
    scale_linetype_manual(values = lt_vals, name = "Metric",
                          labels = .TWD_METRIC_LAB[lt_lvls]) +
    scale_x_continuous(breaks = seq(2, 16, by = 2)) +
    scale_y_continuous(limits = c(0, 100),
                       labels = function(x) paste0(x, "%")) +
    labs(x = "Simulation scenario", y = "Rate") +
    .twd_theme(base_size) +
    guides(colour   = guide_legend(nrow = 1, order = 1),
           linetype = guide_legend(nrow = 1, order = 2))

  if (save) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    stub <- file.path(out_dir, "fig_diag_threeway_grid_all")
    ggsave(paste0(stub, ".png"), p, width = 12, height = 6, dpi = 300)
    ggsave(paste0(stub, ".pdf"), p, width = 12, height = 6)
    message("saved: ", stub, ".{png,pdf}")
  }
  p
}
