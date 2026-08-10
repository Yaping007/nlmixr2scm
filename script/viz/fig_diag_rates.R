# ============================================================================
# fig_diag_rates.R  --  VAE covariate-selection RUN-QUALITY curves
# ----------------------------------------------------------------------------
# Per-scenario diagnostic rates across the 16 scenarios, faceted by model
# parameterization (structure: linCmt / ode), coloured by sample size
# (sample_N: 40 / 80 / 300).  Companion to fig_power.R (same house style).
#
# Reads the aggregated CSV written by aggregate_vae_covsel.R:
#   <agg_dir>/vae_diag_rates.csv
#   cols: sample_N, scenario, structure, n_total,
#         Converged_pct        = % of fits that minimised (min. successful)
#         CNBelowCutoff_pct    = % with cond#(cor) < 1000  (well-conditioned)
#         ConvergedStrict_pct  = % converged AND CN < 1000 (both gates)
#         CovStep_pct, MinSuc_pct, ... (also selectable)
#   NOTE: these columns are already PERCENTAGES (0-100), unlike fig_power.R
#         whose Power* columns are fractions (0-1).
#
# METRICS (any one or several):
#   "Converged"       convergence / minimisation-success rate
#   "CNBelowCutoff"   condition-number stability rate (CN(cor) < 1000)
#   "ConvergedStrict" both gates jointly
#   "CovStep"         covariance-step success rate
#   "MinSuc"          minimisation-success flag rate
#     * ONE metric  -> colour = sample_N (single line per N).
#     * MANY metrics-> colour = sample_N, LINETYPE = metric.
#
# KNOBS (all optional):
#   agg_dir    : aggregation folder (default the current full run)
#   metric     : which rate(s) to plot (default c("Converged","CNBelowCutoff"))
#   sample_N   : NULL = ALL (40, 80, 300) | subset e.g. c(80, 300)
#   structure  : NULL = ALL (linCmt, ode) | subset e.g. "linCmt"
#   save       : FALSE (return only) | TRUE (also ggsave PNG + PDF)
#
# Returns the ggplot object.
#
# USAGE:
#   source("script/viz/fig_diag_rates.R")
#   fig_diag_rates()                                   # Converged + CN<1000
#   fig_diag_rates(metric = "Converged")
#   fig_diag_rates(metric = "CNBelowCutoff", structure = "linCmt")
#   fig_diag_rates(save = TRUE)
# ============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

# ---- shared look (mirrors fig_power.R) -------------------------------------
.PAL_N <- c("40" = "#7F7F7F", "80" = "#E8820C", "300" = "#1F77B4")  # grey/orange/blue

theme_scm <- function(base_size = 14) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "top",
      legend.title     = ggplot2::element_text(size = base_size),
      legend.text      = ggplot2::element_text(size = base_size - 1),
      strip.text       = ggplot2::element_text(size = base_size),
      axis.title       = ggplot2::element_text(size = base_size),
      axis.text.x      = ggplot2::element_text(size = base_size - 4),
      axis.text.y      = ggplot2::element_text(size = base_size - 2),
      plot.title       = ggplot2::element_text(size = base_size - 2)
    )
}

# advan4 FIRST so it is the top facet panel (analytic-on-top convention, matching
# linCmt-on-top for the nlmixr2 figures). PsN data carries advan4 + ode; nlmixr2
# data carries linCmt + ode -- unused levels are dropped by the facet.
.STRUCT_LAB <- c(advan4 = "ADVAN4 (analytic)", linCmt = "linCmt (analytic)",
                 ode = "ADVAN13 (ODE)")

# NOTE: diag-UNIQUE names (.DIAG_*) so sourcing fig_power.R -- which defines its
# own global `.METRIC_LAB` / `.PAL_METRIC_LT` with DIFFERENT (Power*) keys -- can
# never clobber these. A shared name previously caused fig_diag_rates() to look
# up Converged/CNBelowCutoff/ConvergedStrict in the Power palette (-> all NA
# linetypes -> every connecting line dropped) when both files were sourced.
.DIAG_METRIC_LAB <- c(
  Converged       = "Converged (%)",
  CNBelowCutoff   = "Cond# < 1000 (%)",
  ConvergedStrict = "Converged & Cond# < 1000 (%)",
  CovStep         = "Covariance step OK (%)",
  MinSuc          = "Minimisation success (%)"
)

# column in the CSV backing each metric (all already 0-100 percentages)
.METRIC_COL <- c(
  Converged       = "Converged_pct",
  CNBelowCutoff   = "CNBelowCutoff_pct",
  ConvergedStrict = "ConvergedStrict_pct",
  CovStep         = "CovStep_pct",
  MinSuc          = "MinSuc_pct"
)

# metric linetypes (multi-metric overlay mode only); diag-UNIQUE name.
.DIAG_METRIC_LT <- c(Converged = "solid", CNBelowCutoff = "22",
                     ConvergedStrict = "42", CovStep = "44", MinSuc = "13")

# ---- main ------------------------------------------------------------------
fig_diag_rates <- function(agg_dir   = "output/vae_covsel_aggregated",
                           metric    = c("Converged", "CNBelowCutoff"),
                           sample_N  = NULL,   # NULL = all
                           structure = NULL,   # NULL = all
                           estimator = NULL,   # SCM only (e.g. "ifocei"); needs col
                           outer_opt = NULL,   # SCM only (e.g. "bobyqa"); needs col
                           ylim      = c(0, 100),  # zoom e.g. c(90, 100)
                           save      = FALSE,
                           out_dir   = "output/figures/vae_covsel",
                           csv_name     = "vae_diag_rates.csv",
                           title_prefix = "VAE covariate selection") {

  valid_metrics <- names(.METRIC_COL)
  if (missing(metric)) metric <- c("Converged", "CNBelowCutoff")
  metric <- unique(metric)
  bad <- setdiff(metric, valid_metrics)
  if (length(bad)) stop("unknown metric(s): ", paste(bad, collapse = ", "),
                        " -- choose from ", paste(valid_metrics, collapse = ", "))
  multi <- length(metric) > 1L

  csv <- file.path(trimws(agg_dir), csv_name)
  if (!file.exists(csv)) stop("diag-rates CSV not found: ", csv)
  dat <- readr::read_csv(csv, show_col_types = FALSE)

  if (!is.null(sample_N))  dat <- dplyr::filter(dat, sample_N  %in% !!sample_N)
  if (!is.null(structure)) dat <- dplyr::filter(dat, structure %in% !!structure)
  # estimator / outer_opt filters (SCM diag CSV only; VAE CSV lacks these cols)
  if (!is.null(estimator)) {
    if (!"estimator" %in% names(dat))
      stop("estimator filter requested but 'estimator' column absent in ", csv_name)
    dat <- dplyr::filter(dat, estimator %in% !!estimator)
  }
  if (!is.null(outer_opt)) {
    if (!"outer_opt" %in% names(dat))
      stop("outer_opt filter requested but 'outer_opt' column absent in ", csv_name)
    dat <- dplyr::filter(dat, outer_opt %in% !!outer_opt)
  }
  if (nrow(dat) == 0L) stop("no rows after filtering (check sample_N / structure / estimator / outer_opt)")
  # guard: multiple estimator x outer_opt cells would blend distinct estimators
  # onto one line (grouping is only sample_N x scenario x structure).
  if (all(c("estimator", "outer_opt") %in% names(dat))) {
    ncell <- nrow(unique(dat[, c("estimator", "outer_opt")]))
    if (ncell > 1L)
      stop("diag CSV holds ", ncell, " estimator x outer_opt cells; ",
           "pass estimator= / outer_opt= to select one (got: ",
           paste(unique(paste(dat$estimator, dat$outer_opt, sep = "_")),
                 collapse = ", "), ")")
  }

  # platform-aware structure labels: PsN/NONMEM (estimator == "nonmem_scm")
  # renders ode as "ADVAN13 (ODE)"; nlmixr2 (focei/vae) renders plain "ODE".
  is_nonmem  <- "estimator" %in% names(dat) && any(dat$estimator == "nonmem_scm")
  struct_lab <- .STRUCT_LAB
  if (!is_nonmem) struct_lab["ode"] <- "ODE"

  metric_cols <- unname(.METRIC_COL[metric])
  plot_df <- dat |>
    dplyr::select(sample_N, scenario, structure, dplyr::all_of(metric_cols)) |>
    tidyr::pivot_longer(dplyr::all_of(metric_cols),
                        names_to = "metric_col", values_to = "rate_pct") |>
    dplyr::mutate(
      metric    = factor(names(.METRIC_COL)[match(metric_col, .METRIC_COL)],
                         levels = valid_metrics),
      scenario  = factor(scenario, levels = sort(unique(scenario))),
      sample_N  = factor(sample_N, levels = c(40, 80, 300)),
      structure = factor(structure, levels = names(struct_lab))
    )

  y_lab <- if (multi) "Rate (%)" else .DIAG_METRIC_LAB[[metric]]

  ylim   <- range(ylim)
  .span  <- diff(ylim)
  .step  <- if (.span <= 15) 2 else if (.span <= 40) 5 else 20
  .brks  <- seq(floor(ylim[1] / .step) * .step,
                ceiling(ylim[2] / .step) * .step, by = .step)

  aes_base <- if (multi) {
    ggplot2::aes(scenario, rate_pct, colour = sample_N,
                 linetype = metric, group = interaction(sample_N, metric))
  } else {
    ggplot2::aes(scenario, rate_pct, colour = sample_N, group = sample_N)
  }

  p <- ggplot2::ggplot(plot_df, aes_base) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ structure, nrow = 1,
                        labeller = ggplot2::labeller(structure = struct_lab)) +
    ggplot2::scale_colour_manual(values = .PAL_N, name = "Sample size",
                                 labels = function(x) paste0(x, " subj")) +
    ggplot2::scale_y_continuous(breaks = .brks,
                                labels = function(x) paste0(x, "%")) +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::labs(x = "Simulation scenario", y = y_lab) +
    theme_scm()

  if (ylim[1] <= 80 && ylim[2] >= 80) {
    p <- p + ggplot2::geom_hline(yintercept = 80, linetype = "dashed",
                                 colour = "grey40", linewidth = 0.4)
  }

  if (multi) {
    p <- p + ggplot2::scale_linetype_manual(
      values = .DIAG_METRIC_LT[metric], name = "Metric",
      labels = .DIAG_METRIC_LAB[metric])
  }

  # Put the two top legends (Sample size + Metric) SIDE BY SIDE on a single row.
  # (Kept identical to fig_power.R so convergence and power legends match.)
  p <- p + ggplot2::theme(
      legend.box          = "horizontal",
      legend.box.just     = "left",
      legend.box.spacing  = ggplot2::unit(2, "pt"),
      legend.spacing.x    = ggplot2::unit(4, "pt"),
      legend.margin       = ggplot2::margin(0, 0, 0, 0),
      legend.key.size     = ggplot2::unit(14, "pt"),
      legend.title        = ggplot2::element_text(size = 12),
      legend.text         = ggplot2::element_text(size = 12)
    ) +
    ggplot2::guides(colour   = ggplot2::guide_legend(nrow = 1),
                    linetype = ggplot2::guide_legend(nrow = 1))

  if (isTRUE(save)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    tag   <- if (multi) paste(metric, collapse = "-") else metric
    stru  <- if (is.null(structure)) "all" else paste(structure, collapse = "-")
    est_tag <- if (all(c("estimator", "outer_opt") %in% names(dat)))
                 paste0("_", dat$estimator[1], "_", dat$outer_opt[1]) else ""
    stub  <- file.path(out_dir, sprintf("fig_diag_rates_%s_%s%s", stru, tag, est_tag))
    ggplot2::ggsave(paste0(stub, ".png"), p, width = 9, height = 4.5, dpi = 150)
    ggplot2::ggsave(paste0(stub, ".pdf"), p, width = 9, height = 4.5)
    message("saved: ", stub, ".{png,pdf}")
  }

  p
}
