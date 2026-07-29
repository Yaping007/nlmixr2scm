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
    estimator  = NULL,     # SCM only (e.g. "ifocei"); needs col
    outer_opt  = NULL,     # SCM only (e.g. "bobyqa"); needs col
    metrics    = c("RMRSE", "MARE"),
    drop_fixed = TRUE,     # drop TVKA (fixed -> 0 error, uninformative)
    labels     = FALSE,
    save       = FALSE,
    out_dir    = "output/figures/vae_covsel",
    csv_all    = "vae_estim_all.csv",   # SCM: "scm_estim_all.csv"
    csv_cond   = "vae_estim_cond.csv") {  # SCM: "scm_estim_cond.csv"

  metrics <- match.arg(metrics, c("RMRSE", "MARE"), several.ok = TRUE)

  # ---- load both conditioning regimes and stack -----------------------------
  read_estim <- function(file, cond_label) {
    path <- file.path(trimws(agg_dir), file)
    if (!file.exists(path)) stop("estim CSV not found: ", path)
    readr::read_csv(path, show_col_types = FALSE) |>
      dplyr::mutate(conditioning = cond_label)
  }
  dat <- dplyr::bind_rows(
    read_estim(csv_all,  "Unconditioned"),
    read_estim(csv_cond, "True selection")
  )

  scenario <- scenario[1]
  dat <- dplyr::filter(dat, scenario == !!scenario)
  if (!is.null(sample_N))  dat <- dplyr::filter(dat, sample_N %in% !!sample_N)
  if (!is.null(structure)) dat <- dplyr::filter(dat, structure %in% !!structure)
  # estimator / outer_opt filters (SCM estim CSVs only; VAE CSVs lack these cols)
  if (!is.null(estimator)) {
    if (!"estimator" %in% names(dat))
      stop("estimator filter requested but 'estimator' column absent in ", csv_all)
    dat <- dplyr::filter(dat, estimator %in% !!estimator)
  }
  if (!is.null(outer_opt)) {
    if (!"outer_opt" %in% names(dat))
      stop("outer_opt filter requested but 'outer_opt' column absent in ", csv_all)
    dat <- dplyr::filter(dat, outer_opt %in% !!outer_opt)
  }
  if (drop_fixed)          dat <- dplyr::filter(dat, parameter != "TVKA")
  if (nrow(dat) == 0L) stop("no rows after filtering (check scenario / sample_N / structure / estimator / outer_opt)")
  # guard: multiple estimator x outer_opt cells would blend distinct estimators
  # at the same (parameter, conditioning) point.
  if (all(c("estimator", "outer_opt") %in% names(dat))) {
    ncell <- nrow(unique(dat[, c("estimator", "outer_opt")]))
    if (ncell > 1L)
      stop("estim CSV holds ", ncell, " estimator x outer_opt cells; ",
           "pass estimator= / outer_opt= to select one (got: ",
           paste(unique(paste(dat$estimator, dat$outer_opt, sep = "_")),
                 collapse = ", "), ")")
  }

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

  # ---- build faceting: metric (rows) x [N x structure] (cols) ---------------
  # metric is the ONLY row facet so RMRSE stacks cleanly ABOVE MARE; the
  # param_class grouping is preserved by the y-axis ordering (no right strip).
  n_N      <- dplyr::n_distinct(plot_df$sample_N)
  n_struct <- dplyr::n_distinct(plot_df$structure)
  col_terms <- c(if (n_N > 1)      "sampN_lab",
                 if (n_struct > 1) "struct_lab")
  if (length(col_terms) == 0L) col_terms <- "."
  facet_spec <- as.formula(paste("metric ~", paste(col_terms, collapse = " + ")))

  # ---- per-metric subplots stacked vertically -------------------------------
  # Each metric (RMRSE / MARE) is its own panel with the metric name as a
  # CENTRED TOP title, so the two blocks separate cleanly; N (x structure) are
  # the columns.  Panels are stitched with patchwork (shared legend on top).
  col_spec <- as.formula(paste("~", paste(col_terms, collapse = " + ")))
  .METRIC_SUB <- c(
    RMRSE = "100\u00b7\u221amean(rel.err\u00b2)  \u2014 outlier-sensitive",
    MARE  = "100\u00b7median|rel.err|  \u2014 robust")
  dodge <- ggplot2::position_dodge(width = 0.6)

  build_one <- function(m, first = FALSE, last = FALSE) {
    d <- dplyr::filter(plot_df, metric == m)
    g <- ggplot2::ggplot(d,
                         ggplot2::aes(x = value, y = param_lab,
                                      colour = conditioning, group = conditioning)) +
      ggplot2::geom_linerange(
        ggplot2::aes(xmin = 0, xmax = value),
        position = dodge, linewidth = 0.7, alpha = 0.55) +
      ggplot2::geom_point(position = dodge, size = 2.6) +
      ggplot2::scale_colour_manual(values = .COND_COL, name = NULL) +
      ggplot2::facet_grid(col_spec, scales = "free_y", space = "free_y") +
      ggplot2::scale_x_continuous(labels = function(x) paste0(x, "%"),
                                  limits = c(0, 80),
                                  breaks = seq(0, 80, 20),
                                  expand = ggplot2::expansion(mult = c(0, 0.02))) +
      ggplot2::labs(title = as.character(m),
                    subtitle = unname(.METRIC_SUB[as.character(m)]),
                    x = if (last) "Relative error metric (%)" else NULL,
                    y = NULL) +
      theme_scm() +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5,
                                              size = ggplot2::rel(1.15)),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey35"),
        # add a gap between the N columns so adjacent 100%/0% ticks don't collide
        panel.spacing.x = ggplot2::unit(1.4, "lines"),
        # show the N (column) headers only on the TOP block to avoid repetition
        strip.text.x  = if (first) ggplot2::element_text(face = "bold")
                        else ggplot2::element_blank())
    if (isTRUE(labels)) {
      g <- g + ggplot2::geom_text(
        ggplot2::aes(label = round(value)),
        position = dodge, hjust = -0.3, size = 3, show.legend = FALSE)
    }
    g
  }

  n_m    <- length(metrics)
  panels <- lapply(seq_along(metrics), function(i)
    build_one(metrics[i], first = (i == 1L), last = (i == n_m)))

  .cap <- paste0(
    "Structural intercepts (TVCL, TVVc) are back-transformed to the ",
    "data-generating reference (BW = 70, CrCL = 95) using each fit's own ",
    "power betas, so their error is comparable to truth.")

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p <- patchwork::wrap_plots(panels, ncol = 1) +
      patchwork::plot_layout(guides = "collect") +
      patchwork::plot_annotation(
        theme = ggplot2::theme(legend.position = "top"))
  } else {
    # fallback (no patchwork): single grid with metric row strips
    p <- ggplot2::ggplot(plot_df,
                         ggplot2::aes(x = value, y = param_lab,
                                      colour = conditioning, group = conditioning)) +
      ggplot2::geom_linerange(ggplot2::aes(xmin = 0, xmax = value),
                              position = dodge, linewidth = 0.7, alpha = 0.55) +
      ggplot2::geom_point(position = dodge, size = 2.6) +
      ggplot2::scale_colour_manual(values = .COND_COL, name = NULL) +
      ggplot2::facet_grid(facet_spec, scales = "free_y", space = "free_y") +
      ggplot2::scale_x_continuous(labels = function(x) paste0(x, "%"),
                                  expand = ggplot2::expansion(mult = c(0, 0.08))) +
      ggplot2::labs(
        title = sprintf("Estimation accuracy by parameter -- scenario %s", scenario),
        x = "Relative error metric (%)", y = NULL, caption = .cap) +
      theme_scm() +
      ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0, face = "bold"))
    if (isTRUE(labels)) {
      p <- p + ggplot2::geom_text(ggplot2::aes(label = round(value)),
        position = dodge, hjust = -0.3, size = 3, show.legend = FALSE)
    }
  }

  if (isTRUE(save)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    n_tag <- if (is.null(sample_N)) "allN" else paste0("N", paste(sample_N, collapse = "-"))
    s_tag <- if (is.null(structure)) "allStruct" else paste(structure, collapse = "-")
    m_tag <- paste(metrics, collapse = "-")
    est_tag <- if (all(c("estimator", "outer_opt") %in% names(dat)))
                 paste0("_", dat$estimator[1], "_", dat$outer_opt[1]) else ""
    stub  <- file.path(out_dir,
                       sprintf("fig_error_metrics_scn%s_%s_%s_%s%s",
                               scenario, n_tag, s_tag, m_tag, est_tag))
    # width scales with N x structure cols; height with metric (row) count
    n_cols <- max(1, n_N) * max(1, n_struct)
    n_rows <- length(metrics)
    w <- 4 + 3.2 * n_cols
    h <- 5 + 3.0 * n_rows
    ggplot2::ggsave(paste0(stub, ".png"), p, width = w, height = h, dpi = 200)
    ggplot2::ggsave(paste0(stub, ".pdf"), p, width = w, height = h)
    message("saved: ", stub, ".{png,pdf}")
  }

  p
}
