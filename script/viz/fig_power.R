# ============================================================================
# fig_power.R  --  VAE covariate-selection POWER curves
# ----------------------------------------------------------------------------
# Power curve across the 16 scenarios, faceted by model parameterization
# (structure: linCmt / ode), coloured by sample size (sample_N: 40 / 80 / 300).
#
# Reads the aggregated CSV written by aggregate_vae_covsel.R:
#   output/vae_covsel_aggregated/vae_power.csv
#   cols: sample_N, scenario, structure, N, n_exact, ...,
#         Power, PowerCN, PowerMinSuc   (Power* are FRACTIONS 0-1)
#
# KNOBS (all optional):
#   metric     : which power definition(s) to plot. Accepts ONE or SEVERAL of
#                "Power"       = exact-match / N              (default)
#                "PowerCN"     = exact-match among CN < cutoff fits
#                "PowerMinSuc" = exact-match among converged fits
#                * ONE metric  -> colour = sample_N (single line per N).
#                * MANY metrics-> colour = sample_N, LINETYPE = metric. When
#                  the CN / convergence gates are non-binding the three
#                  linetypes overlap exactly (visual proof they are equal).
#                  Pass metric = c("Power","PowerCN","PowerMinSuc") for all.
#   sample_N   : NULL = ALL (40, 80, 300)  |  subset e.g. c(80, 300)
#   structure  : NULL = ALL (linCmt, ode)  |  subset e.g. "linCmt"
#   save       : FALSE (return only) | TRUE (also ggsave PNG+PDF)
#
# DEFAULT: fig_power() shows ALL sample_N x ALL structure (6 series total:
#          3 coloured lines x 2 facet rows), metric = Power.
#
# Returns the ggplot object (testing phase).
#
# USAGE (interactive):
#   source("script/viz/fig_power.R")
#   fig_power()                                   # all N, both structures, Power
#   fig_power(metric = "PowerMinSuc")
#   fig_power(metric = c("Power","PowerCN","PowerMinSuc"))  # overlay by linetype
#   fig_power(sample_N = c(80, 300))              # focus on larger N
#   fig_power(structure = "linCmt")              # focus one parameterization
#   fig_power(save = TRUE)                        # write PNG + PDF
# ============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

# ---- shared look (local for now; promotable to viz_common.R later) ---------
.PAL_N <- c("40" = "#7F7F7F", "80" = "#E8820C", "300" = "#1F77B4")  # grey/orange/blue

theme_scm <- function(base_size = 14) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "top",
      legend.title     = ggplot2::element_text(size = base_size),
      legend.text      = ggplot2::element_text(size = base_size - 1),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      strip.text       = ggplot2::element_text(size = base_size),
      axis.title       = ggplot2::element_text(size = base_size),
      axis.text.x      = ggplot2::element_text(size = base_size - 4),
      axis.text.y      = ggplot2::element_text(size = base_size - 2),
      plot.title       = ggplot2::element_text(size = base_size - 2)
    )
}

# advan4 FIRST -> top facet panel (analytic-on-top convention).
.STRUCT_LAB <- c(advan4 = "ADVAN4 (analytic)", linCmt = "linCmt (analytic)",
                 ode = "ADVAN13 (ODE)")

.METRIC_LAB <- c(Power       = "Power",
                 PowerCN     = "PowerCN",
                 PowerMinSuc = "PowerMinSuc")

# metric linetypes (used only in multi-metric overlay mode)
.PAL_METRIC_LT <- c(Power = "solid", PowerCN = "22", PowerMinSuc = "42")

# ---- main ------------------------------------------------------------------
fig_power <- function(agg_dir   = "output/vae_covsel_aggregated",
                      metric    = c("Power", "PowerCN", "PowerMinSuc"),
                      sample_N  = NULL,   # NULL = all
                      structure = NULL,   # NULL = all
                      estimator = NULL,   # NULL = all (e.g. "focei"); needs col
                      outer_opt = NULL,   # NULL = all (e.g. "bobyqa"); needs col
                      save      = FALSE,
                      out_dir   = "output/figures/vae_covsel",
                      csv_name     = "vae_power.csv",  # "scm_power.csv" for SCM
                      title_prefix = "VAE covariate selection") {

  # metric may be ONE (default "Power") or SEVERAL. Do NOT use match.arg here so
  # that a length>1 vector is allowed; validate manually instead.
  valid_metrics <- c("Power", "PowerCN", "PowerMinSuc")
  if (missing(metric)) metric <- "Power"
  metric <- unique(metric)
  bad <- setdiff(metric, valid_metrics)
  if (length(bad)) stop("unknown metric(s): ", paste(bad, collapse = ", "),
                        " -- choose from ", paste(valid_metrics, collapse = ", "))
  multi <- length(metric) > 1L

  csv    <- file.path(trimws(agg_dir), csv_name)
  if (!file.exists(csv)) stop("power CSV not found: ", csv)

  dat <- readr::read_csv(csv, show_col_types = FALSE)

  # optional filtering knobs (default NULL => keep everything)
  if (!is.null(sample_N))  dat <- dplyr::filter(dat, sample_N  %in% !!sample_N)
  if (!is.null(structure)) dat <- dplyr::filter(dat, structure %in% !!structure)
  # estimator / outer_opt filters (SCM power CSV only; VAE CSV lacks these cols)
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
  if (nrow(dat) == 0L)
    stop("no rows after filtering (check sample_N / structure / estimator / outer_opt)")
  # guard: if multiple estimator/outer_opt cells remain, the (N,scenario,
  # structure) grouping would double-plot -- force the caller to disambiguate.
  if (all(c("estimator", "outer_opt") %in% names(dat))) {
    ncell <- nrow(unique(dat[, c("estimator", "outer_opt")]))
    if (ncell > 1L)
      stop("power CSV holds ", ncell, " estimator x outer_opt cells; ",
           "pass estimator= / outer_opt= to select one (got: ",
           paste(unique(paste(dat$estimator, dat$outer_opt, sep = "_")),
                 collapse = ", "), ")")
  }

  # platform-aware structure labels: PsN/NONMEM (estimator == "nonmem_scm")
  # renders ode as "ADVAN13 (ODE)"; nlmixr2 (focei/vae) renders plain "ODE".
  is_nonmem  <- "estimator" %in% names(dat) && any(dat$estimator == "nonmem_scm")
  struct_lab <- .STRUCT_LAB
  if (!is_nonmem) struct_lab["ode"] <- "ODE"

  # long over the requested metrics
  plot_df <- dat |>
    dplyr::select(sample_N, scenario, structure, dplyr::all_of(metric)) |>
    tidyr::pivot_longer(dplyr::all_of(metric),
                        names_to = "metric", values_to = "power") |>
    dplyr::transmute(
      scenario  = factor(scenario, levels = sort(unique(scenario))),
      sample_N  = factor(sample_N, levels = c(40, 80, 300)),
      structure = factor(structure, levels = names(struct_lab)),
      metric    = factor(metric, levels = valid_metrics),
      power_pct = power * 100
    )

  y_lab   <- if (multi) "Power" else .METRIC_LAB[[metric]]
  title_m <- if (multi) "Power (all definitions)" else .METRIC_LAB[[metric]]

  aes_base <- if (multi) {
    ggplot2::aes(scenario, power_pct, colour = sample_N,
                 linetype = metric, group = interaction(sample_N, metric))
  } else {
    ggplot2::aes(scenario, power_pct, colour = sample_N, group = sample_N)
  }

  p <- ggplot2::ggplot(plot_df, aes_base) +
    ggplot2::geom_hline(yintercept = 80, linetype = "dashed",
                        colour = "grey40", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ structure, nrow = 1,
                        labeller = ggplot2::labeller(structure = struct_lab)) +
    ggplot2::scale_colour_manual(values = .PAL_N, name = "Sample size",
                                 labels = function(x) paste0(x, " subj")) +
    ggplot2::scale_y_continuous(limits = c(0, 100),
                                breaks = seq(0, 100, 20),
                                labels = function(x) paste0(x, "%")) +
    ggplot2::labs(
      x = "Simulation scenario", y = y_lab
    ) +
    theme_scm()

  if (multi) {
    p <- p + ggplot2::scale_linetype_manual(
      values = .PAL_METRIC_LT, name = "Power definition",
      labels = .METRIC_LAB[metric])
  }

  # Put the two top legends (Sample size + Power definition) SIDE BY SIDE on a
  # single row (matching the convergence figure).
  p <- p + ggplot2::theme(
      legend.box          = "horizontal",
      legend.box.just     = "left",
      legend.box.spacing  = ggplot2::unit(2, "pt"),
      legend.spacing.x    = ggplot2::unit(4, "pt"),
      legend.margin       = ggplot2::margin(0, 0, 0, 0),
      legend.key.size     = ggplot2::unit(14, "pt"),
      legend.title        = ggplot2::element_text(size = 11.5),
      legend.text         = ggplot2::element_text(size = 11.5)
    ) +
    ggplot2::guides(colour   = ggplot2::guide_legend(nrow = 1),
                    linetype = ggplot2::guide_legend(nrow = 1))

  if (isTRUE(save)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    tag  <- if (multi) "multi" else metric
    cell <- paste(c(estimator, outer_opt), collapse = "_")
    if (nzchar(cell)) tag <- paste(cell, tag, sep = "_")
    stub <- file.path(out_dir, sprintf("fig_power_%s", tag))
    ggplot2::ggsave(paste0(stub, ".png"), p, width = 12, height = 4.5, dpi = 150)
    ggplot2::ggsave(paste0(stub, ".pdf"), p, width = 12, height = 4.5)
    message("saved: ", stub, ".{png,pdf}")
  }

  p
}
