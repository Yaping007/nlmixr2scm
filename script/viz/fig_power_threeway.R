## fig_power_threeway.R
## Three-way method comparison of SCM covariate-selection POWER at a fixed
## sample size (default N = 80).
##
##   Figure A  fig_power_threeway_N80()   three methods OVERLAID as coloured
##                                        lines; facets = Analytic | ODE.
##   Figure B  fig_power_threeway_grid()  3x2 small-multiples; methods as
##                                        COLUMNS, structure as ROWS.
##
## Methods
##   nlmixr2-FOCEi : output/scm_bench_rescue_winner_aggregated/scm_power.csv
##                   (filter estimator == "focei", outer_opt == "bobyqa")
##   NONMEM-PsN    : output/psn_scm_combined_aggregated/scm_power.csv
##                   (filter estimator == "nonmem_scm", outer_opt == "focei")
##   nlmixr2-VAE   : output/vae_covsel_full0729_est702_aggregated/vae_power.csv
##                   (no estimator / outer_opt columns)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

## ---- configuration ---------------------------------------------------------

.TW_SOURCES <- list(
  "nlmixr2-SCM" = list(
    path      = "output/scm_bench_rescue_winner_aggregated/scm_power.csv",
    estimator = "focei",
    outer_opt = "bobyqa"
  ),
  "PsN-SCM" = list(
    path      = "output/psn_scm_combined_aggregated/scm_power.csv",
    estimator = "nonmem_scm",
    outer_opt = "focei"
  ),
  "nlmixr2-VAE" = list(
    path      = "output/vae_covsel_full0729_est702_aggregated/vae_power.csv",
    estimator = NULL,
    outer_opt = NULL
  )
)

## fixed method order + palette
.TW_LEVELS  <- c("nlmixr2-SCM", "PsN-SCM", "nlmixr2-VAE")
.PAL_METHOD <- c(
  "nlmixr2-SCM"  = "#0F8B8D",   # teal
  "PsN-SCM"      = "#2B2B2B",   # near-black
  "nlmixr2-VAE"  = "#A64AC9"    # purple
)

## ---- loader ----------------------------------------------------------------

.tw_load <- function(sources = .TW_SOURCES,
                     sample_N = 80,
                     metric   = "Power") {
  rows <- lapply(names(sources), function(m) {
    s <- sources[[m]]
    d <- readr::read_csv(s$path, show_col_types = FALSE)

    if (!is.null(s$estimator) && "estimator" %in% names(d)) {
      d <- dplyr::filter(d, .data$estimator == s$estimator)
    }
    if (!is.null(s$outer_opt) && "outer_opt" %in% names(d)) {
      d <- dplyr::filter(d, .data$outer_opt == s$outer_opt)
    }
    d <- dplyr::filter(d, .data$sample_N == !!sample_N)

    tibble::tibble(
      method    = m,
      scenario  = d$scenario,
      structure = d$structure,
      power_pct = d[[metric]] * 100
    )
  })

  dplyr::bind_rows(rows) |>
    dplyr::mutate(
      structure = dplyr::case_when(
        structure %in% c("linCmt", "advan4") ~ "Analytic",
        structure == "ode"                    ~ "ODE",
        TRUE                                  ~ structure
      ),
      structure = factor(structure, levels = c("Analytic", "ODE")),
      method    = factor(method, levels = .TW_LEVELS)
    )
}

## power-definition labels / linetypes (match fig_power.R convention)
.TW_METRIC_LAB <- c(Power       = "Power",
                    PowerCN     = "PowerCN",
                    PowerMinSuc = "PowerMinSuc")
.TW_METRIC_LT  <- c(Power = "solid", PowerCN = "22", PowerMinSuc = "42")

## sample-size palette (identical to fig_power.R: grey / orange / blue)
.TW_PAL_N <- c("40" = "#7F7F7F", "80" = "#E8820C", "300" = "#1F77B4")

## full loader: ALL sample sizes, ALL power definitions (long format)
.tw_load_full <- function(sources = .TW_SOURCES,
                          metrics = c("Power", "PowerCN", "PowerMinSuc")) {
  rows <- lapply(names(sources), function(m) {
    s <- sources[[m]]
    d <- readr::read_csv(s$path, show_col_types = FALSE)

    if (!is.null(s$estimator) && "estimator" %in% names(d)) {
      d <- dplyr::filter(d, .data$estimator == s$estimator)
    }
    if (!is.null(s$outer_opt) && "outer_opt" %in% names(d)) {
      d <- dplyr::filter(d, .data$outer_opt == s$outer_opt)
    }

    have <- intersect(metrics, names(d))
    d |>
      dplyr::select(dplyr::all_of(c("sample_N", "scenario",
                                    "structure", have))) |>
      tidyr::pivot_longer(dplyr::all_of(have),
                          names_to = "metric", values_to = "power") |>
      dplyr::mutate(method = m)
  })

  dplyr::bind_rows(rows) |>
    dplyr::filter(!is.na(.data$power)) |>
    dplyr::mutate(
      structure = dplyr::case_when(
        structure %in% c("linCmt", "advan4") ~ "Analytic",
        structure == "ode"                    ~ "ODE",
        TRUE                                  ~ structure
      ),
      structure = factor(structure, levels = c("Analytic", "ODE")),
      method    = factor(method, levels = .TW_LEVELS),
      power_pct = .data$power * 100,
      metric    = factor(unname(.TW_METRIC_LAB[.data$metric]),
                         levels = unname(.TW_METRIC_LAB)),
      N_lab     = factor(as.character(.data$sample_N),
                         levels = as.character(sort(unique(.data$sample_N))))
    )
}

## ---- shared theme ----------------------------------------------------------

.tw_theme <- function(base_size = 16) {
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

## ---- Figure A : overlaid lines, Analytic | ODE facets ----------------------

fig_power_threeway_N80 <- function(sources   = .TW_SOURCES,
                                   sample_N  = 80,
                                   metric    = "Power",
                                   save      = FALSE,
                                   out_dir   = "output/figures/threeway_N80",
                                   base_size = 16) {
  dat <- .tw_load(sources, sample_N, metric)

  p <- ggplot(dat, aes(scenario, power_pct,
                       colour = method, group = method)) +
    geom_hline(yintercept = 80, linetype = "dashed",
               colour = "grey60") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    facet_wrap(~ structure, nrow = 1) +
    scale_colour_manual(values = .PAL_METHOD, name = "Method") +
    scale_x_continuous(breaks = 1:16) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(
      x     = "Scenario",
      y     = "Power (%)"
    ) +
    .tw_theme(base_size) +
    guides(colour = guide_legend(nrow = 1))

  if (save) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    stub <- file.path(out_dir,
      sprintf("fig_power_threeway_N%d_%s", sample_N, metric))
    ggsave(paste0(stub, ".png"), p, width = 12, height = 4.5, dpi = 300)
    ggsave(paste0(stub, ".pdf"), p, width = 12, height = 4.5)
    message("saved: ", stub, ".{png,pdf}")
  }
  p
}

## ---- Figure B : 3x2 grid, ALL sample sizes + ALL power definitions ----------
## method columns x structure rows; colour = sample size, linetype = power def.

fig_power_threeway_grid <- function(sources   = .TW_SOURCES,
                                    metrics   = c("Power", "PowerCN",
                                                  "PowerMinSuc"),
                                    save      = FALSE,
                                    out_dir   = "output/figures/threeway_N80",
                                    base_size = 16) {
  dat <- .tw_load_full(sources, metrics)

  n_lvls  <- levels(dat$N_lab)
  pal_N   <- .TW_PAL_N[n_lvls]
  lt_lvls <- levels(dat$metric)
  lt_vals <- setNames(.TW_METRIC_LT[names(.TW_METRIC_LAB)], .TW_METRIC_LAB)[lt_lvls]

  p <- ggplot(dat, aes(scenario, power_pct,
                       colour   = N_lab,
                       linetype = metric,
                       group    = interaction(N_lab, metric))) +
    geom_hline(yintercept = 80, linetype = "dashed",
               colour = "grey70") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.1) +
    facet_grid(structure ~ method) +
    scale_colour_manual(values = pal_N, name = "Sample size",
                        labels = function(x) paste0(x, " subj")) +
    scale_linetype_manual(values = lt_vals, name = "Power definition") +
    scale_x_continuous(breaks = seq(2, 16, by = 2)) +
    scale_y_continuous(limits = c(0, 100),
                       labels = function(x) paste0(x, "%")) +
    labs(
      title = "Covariate-selection power",
      x     = "Simulation scenario",
      y     = "Power"
    ) +
    .tw_theme(base_size) +
    guides(colour   = guide_legend(nrow = 1, order = 1),
           linetype = guide_legend(nrow = 1, order = 2))

  if (save) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    stub <- file.path(out_dir, "fig_power_threeway_grid_all")
    ggsave(paste0(stub, ".png"), p, width = 12, height = 6, dpi = 300)
    ggsave(paste0(stub, ".pdf"), p, width = 12, height = 6)
    message("saved: ", stub, ".{png,pdf}")
  }
  p
}
