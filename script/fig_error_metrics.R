# ============================================================================
# fig_error_metrics.R  --  per-parameter estimation accuracy (RMRSE / MARE),
#                          conditioned vs unconditioned on TRUE covariate
#                          selection.
# ----------------------------------------------------------------------------
# Answers: "for a given scenario, how accurately is EACH population parameter /
#           covariate effect estimated, and does RESTRICTING to datasets where
#           the covariate model was selected CORRECTLY change that accuracy?"
#
# Reads the schema-2.1 estimation roll-ups written by aggregate_vae_covsel.R:
#   <agg_dir>/vae_estim_all.csv    -- UNCONDITIONED (every usable fit)
#   <agg_dir>/vae_estim_cond.csv   -- CONDITIONED on exact_match == TRUE
#                                     (selected covariate set == truth: no FP,
#                                      no FN).  Same schema, fewer datasets.
#   cols: sample_N, scenario, structure, parameter, param_class, true_value,
#         n_used, MedRE_pct, MeanRE_pct, MARE_pct, RMRSE_pct, P90AbsRE_pct
#
# METRICS (both in %, higher = worse), per (cell, parameter):
#   RMRSE_pct = 100 * sqrt(mean(rel_err^2))   -- outlier-sensitive accuracy
#   MARE_pct  = 100 * median(|rel_err|)       -- robust central accuracy
#   where rel_err = (theta_hat - theta_true) / theta_true.
#   NOTE: the repo has NO metric literally named "MASE"; MARE_pct is the robust
#   median-relative-error metric and is what "MASE" refers to here.
#
# PARAMETER GROUPS (param_class):
#   covariate_beta       : CLBW, CLcrCL, VcBW, VcSEX  -- the 4 true covariate
#                          effects (headline; centering-invariant).
#   structural_intercept : TVCL, TVVc  -- CAVEAT: these absorb the covariate
#                          centering shift, so their rel_err is NOT pure bias.
#   structural_other     : TVQ, TVVp, TVKA, var_CL, var_Vc, cov_VcCL, ResErr
#                          -- covariate-free baseline (TVKA fixed -> 0, dropped).
#
# LAYOUT
#   y      = parameter (grouped by param_class down the rows)
#   x      = metric value (%)
#   colour = conditioning (Unconditioned vs True selection)
#   facet  = param_class (rows, free y) x metric (cols: RMRSE | MARE)
#           If several sample sizes / structures are present they are added as
#           extra facet columns automatically.
#
# House style matches fig_covsel_heatmap.R (theme_scm, return-first, save=).
#
# KNOBS (all optional):
#   agg_dir    : aggregation folder (default the 0722 VAE run)
#   scenario   : ONE scenario to show (default 16 = all-effects-active)
#   sample_N   : NULL = all present | subset e.g. 300
#   structure  : NULL = all present | "ode" | "linCmt"
#   metrics    : which metrics to draw (default c("RMRSE","MARE"))
#   drop_fixed : drop parameters that are structurally fixed (TVKA -> 0 error)
#   labels     : TRUE = print the value at each point (default FALSE)
#   save       : FALSE (return only) | TRUE (also ggsave PNG + PDF)
#
# Returns the ggplot object.
#
# USAGE:
#   source("script/viz/fig_error_metrics.R")
#   fig_error_metrics()                       # scn16, 0722 run
#   fig_error_metrics(metrics = "MARE")       # robust metric only
#   fig_error_metrics(labels = TRUE, save = TRUE)
# ============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

.STRUCT_LAB <- c(linCmt = "linCmt (analytic)", ode = "ODE")

# pretty parameter labels (kept close to the raw names for traceability)
.PARAM_LAB <- c(
  CLBW = "CL~BW", CLcrCL = "CL~CrCL", VcBW = "Vc~BW", VcSEX = "Vc~SEX",
  TVCL = "TVCL", TVVc = "TVVc",
  TVQ = "TVQ", TVVp = "TVVp", TVKA = "TVKA",
  var_CL = "var(CL)", var_Vc = "var(Vc)", cov_VcCL = "cov(Vc,CL)",
  ResErr = "ResErr"
)

# canonical top-down order within each param_class block
.PARAM_ORD <- c(
  # covariate_beta (headline)
  "CLBW", "CLcrCL", "VcBW", "VcSEX",
  # structural_intercept
  "TVCL", "TVVc",
  # structural_other
  "TVQ", "TVVp", "var_CL", "var_Vc", "cov_VcCL", "ResErr", "TVKA"
)

.CLASS_LAB <- c(
  covariate_beta       = "Covariate effects",
  structural_intercept = "Structural intercepts",
  structural_other     = "Structural / variance"
)

# two-colour qualitative pair (deliberately NOT the blue/red FP-FN semantics)
.COND_COL <- c("Unconditioned" = "#8C8C8C", "True selection" = "#2C7FB8")

theme_scm <- function(base_size = 15) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      legend.position    = "top",
      strip.text         = ggplot2::element_text(face = "bold"),
      plot.title         = ggplot2::element_text(face = "bold"),
      plot.caption       = ggplot2::element_text(hjust = 0, colour = "grey35"),
      axis.text.x        = ggplot2::element_text(size = base_size - 3),
      axis.text.y        = ggplot2::element_text(size = base_size - 2)
    )
}

# ---- main ------------------------------------------------------------------
fig_error_metrics <- function(
    agg_dir    = "output/vae_covsel_aggregated0722",
    scenario   = 16,
    sample_N   = NULL,     # NULL = all present
    structure  = NULL,     # NULL = all present ("ode" | "linCmt")
    metrics    = c("RMRSE", "MARE"),
    drop_fixed = TRUE,     # drop TVKA (fixed -> 0 error, uninformative)
    labels     = FALSE,
    save       = FALSE,
    out_dir    = "output/figures/vae_covsel") {

  metrics <- match.arg(metrics, c("RMRSE", "MARE"), several.ok = TRUE)

  # ---- load both conditioning regimes and stack -----------------------------
  read_estim <- function(file, cond_label) {
    path <- file.path(trimws(agg_dir), file)
    if (!file.exists(path)) stop("estim CSV not found: ", path)
    readr::read_csv(path, show_col_types = FALSE) |>
      dplyr::mutate(conditioning = cond_label)
  }
  dat <- dplyr::bind_rows(
    read_estim("vae_estim_all.csv",  "Unconditioned"),
    read_estim("vae_estim_cond.csv", "True selection")
  )

  scenario <- scenario[1]
  dat <- dplyr::filter(dat, scenario == !!scenario)
  if (!is.null(sample_N))  dat <- dplyr::filter(dat, sample_N %in% !!sample_N)
  if (!is.null(structure)) dat <- dplyr::filter(dat, structure %in% !!structure)
  if (drop_fixed)          dat <- dplyr::filter(dat, parameter != "TVKA")
  if (nrow(dat) == 0L) stop("no rows after filtering (check scenario / sample_N / structure)")

  # ---- long over the requested metrics -------------------------------------
  metric_cols <- c(RMRSE = "RMRSE_pct", MARE = "MARE_pct")[metrics]
  plot_df <- dat |>
    dplyr::select(sample_N, scenario, structure, parameter, param_class,
                  conditioning, dplyr::all_of(unname(metric_cols))) |>
    tidyr::pivot_longer(dplyr::all_of(unname(metric_cols)),
                        names_to = "metric", values_to = "value") |>
    dplyr::mutate(
      metric = factor(names(metric_cols)[match(metric, metric_cols)],
                      levels = metrics),
      param_lab   = dplyr::coalesce(.PARAM_LAB[parameter], parameter),
      class_lab   = dplyr::coalesce(.CLASS_LAB[param_class], param_class),
      class_lab   = factor(class_lab, levels = unname(.CLASS_LAB)),
      conditioning = factor(conditioning,
                            levels = c("Unconditioned", "True selection")),
      struct_lab  = dplyr::coalesce(.STRUCT_LAB[structure], structure),
      sampN_lab   = factor(paste0("N = ", sample_N),
                           levels = paste0("N = ", c(40, 80, 300)))
    )

  # top-down parameter order within each class block (ggplot draws first level
  # at BOTTOM -> reverse so the canonical order reads top-to-bottom)
  p_present <- .PARAM_ORD[.PARAM_ORD %in% unique(plot_df$parameter)]
  lab_present <- dplyr::coalesce(.PARAM_LAB[p_present], p_present)
  plot_df <- dplyr::mutate(plot_df,
    param_lab = factor(param_lab, levels = rev(lab_present)))

  # ---- build faceting: param_class (rows) x metric [x N x structure] (cols) --
  n_N      <- dplyr::n_distinct(plot_df$sample_N)
  n_struct <- dplyr::n_distinct(plot_df$structure)
  col_terms <- c("metric",
                 if (n_N > 1)      "sampN_lab",
                 if (n_struct > 1) "struct_lab")
  facet_spec <- as.formula(paste("class_lab ~", paste(col_terms, collapse = " + ")))

  # ---- plot: lollipop (linerange 0->value + point), dodged by conditioning --
  dodge <- ggplot2::position_dodge(width = 0.6)
  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = value, y = param_lab,
                                    colour = conditioning, group = conditioning)) +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = 0, xmax = value),
      position = dodge, linewidth = 0.7, alpha = 0.55) +
    ggplot2::geom_point(position = dodge, size = 2.6) +
    ggplot2::scale_colour_manual(values = .COND_COL, name = NULL) +
    ggplot2::facet_grid(facet_spec, scales = "free_y", space = "free_y") +
    ggplot2::scale_x_continuous(labels = function(x) paste0(x, "%"),
                                expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::labs(
      title = sprintf("Estimation accuracy by parameter -- scenario %s", scenario),
      subtitle = paste0(
        "Lower is better.  Colour = whether datasets are restricted to those ",
        "with the CORRECT covariate model.\n",
        "RMRSE = 100\u00b7\u221amean(rel.err\u00b2) (outlier-sensitive);  ",
        "MARE = 100\u00b7median|rel.err| (robust)."),
      x = "Relative error metric (%)", y = NULL,
      caption = paste0(
        "Structural intercepts (TVCL, TVVc) are back-transformed to the ",
        "data-generating reference (BW = 70, CrCL = 95) using each fit's own ",
        "power betas, so their error is comparable to truth.")
    ) +
    theme_scm()

  if (isTRUE(labels)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = round(value)),
      position = dodge, hjust = -0.3, size = 3, show.legend = FALSE)
  }

  if (isTRUE(save)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    n_tag <- if (is.null(sample_N)) "allN" else paste0("N", paste(sample_N, collapse = "-"))
    s_tag <- if (is.null(structure)) "allStruct" else paste(structure, collapse = "-")
    m_tag <- paste(metrics, collapse = "-")
    stub  <- file.path(out_dir,
                       sprintf("fig_error_metrics_scn%s_%s_%s_%s",
                               scenario, n_tag, s_tag, m_tag))
    # width scales with number of facet columns; height with parameter count
    n_cols <- length(metrics) * max(1, n_N) * max(1, n_struct)
    w <- 4 + 3.2 * n_cols
    ggplot2::ggsave(paste0(stub, ".png"), p, width = w, height = 9, dpi = 200)
    ggplot2::ggsave(paste0(stub, ".pdf"), p, width = w, height = 9)
    message("saved: ", stub, ".{png,pdf}")
  }

  p
}
