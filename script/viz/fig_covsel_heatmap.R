# ============================================================================
# fig_covsel_heatmap.R  --  VAE covariate-selection error-pattern heatmap
# ----------------------------------------------------------------------------
# Answers: "under each scenario, which covariate EFFECTS tend to become
#           false positive (FP) or false negative (FN)?"
#
# Reads the per-covariate roll-up written by aggregate_vae_covsel.R:
#   output/vae_covsel_aggregated/vae_covsel_by_covar.csv
#   cols: sample_N, scenario, structure, var, covar, shape, n_datasets,
#         is_true, n_detected, detection_rate, n_TP, n_FN, n_FP
#
# SHAPE is now a scored dimension (parity with SCM/PsN): each (var, covar) has a
# TRUE shape (continuous -> "power"; categorical -> "cat") plus DISTRACTOR shapes
# ("lin").  A true covariate recovered in the WRONG shape scores FN on its true
# shape row + FP on the wrong shape row.  Shape is folded into the y-axis effect
# label ("CL~BW.power", "CL~BW.lin", ...) with the distractor shape placed
# directly BELOW its true shape so a flip reads vertically.
#
# Per (scenario x effect) cell there are two mutually-exclusive regimes:
#   * TRUE effect  (is_true = TRUE) : error = FALSE NEGATIVE
#                                     FN rate = n_FN / n_datasets
#   * FALSE effect (is_true = FALSE): error = FALSE POSITIVE
#                                     FP rate = n_FP / n_datasets
#
# DENOMINATOR.  `n_datasets` in vae_covsel_by_covar.csv is now the cell's TRUE
# dataset total (aggregate_vae_covsel.R derives it from the per-cell fit count in
# diag_long, NOT the per-pair appearance count), so n_FP / n_datasets and
# n_FN / n_datasets are correct FP/FN rates for both regimes with no diag-rates
# work-around.  (Historically n_datasets was the appearance count, forcing a
# join to vae_diag_rates.csv$n_total; that quirk is fixed at source.)
# Both are "how often the selector got THIS effect wrong", so a single
# error-rate fill (0-100%, higher = worse) carries them together. TRUE-effect
# cells (the FN regime) get a bold black outline so the two regimes are
# visually separable; FALSE-effect cells (FP regime) are un-outlined.
#
# LAYOUT (redesigned around the BW<->BMI collinearity story)
#   x     = simulation scenario (1..16)
#   y     = covariate effect "PARAM~COVAR" (CL block then VC block; within each
#           block BW & BMI are ADJACENT so the confusion pair reads vertically)
#   fill  = error rate  (FN rate on outlined true cells, FP rate on plain cells)
#   facet = sample_N (columns)
#           A true "param~BW" (outlined = FN) sits directly ABOVE its
#           "param~BMI" thief (= FP): when BW leaks you SEE a dark FN tile
#           stacked on a dark FP tile.  Structure is collapsed to ONE panel set
#           (linCmt ~ ode are near-identical); override via `structure=`.
#
# KNOBS (all optional):
#   sample_N   : NULL = ALL (40, 80, 300) | subset e.g. c(80, 300) | 300
#   structure  : ONE structure to show (default "linCmt"; use "ode" to swap)
#   labels     : TRUE  = print the error % inside each tile (default)
#                FALSE = fill only (cleaner for the full 3x2 facet grid)
#   min_fp     : FP-regime cells with error rate below this are blanked to
#                grey (de-clutter near-zero noise). Default 0 = show all.
#   save       : FALSE (return only) | TRUE (also ggsave PNG + PDF)
#
# Returns the ggplot object.
#
# USAGE (interactive):
#   source("script/viz/fig_covsel_heatmap.R")
#   fig_covsel_heatmap()                         # all N x both structures
#   fig_covsel_heatmap(sample_N = 300)           # focus large cohort
#   fig_covsel_heatmap(structure = "linCmt", labels = FALSE)
#   fig_covsel_heatmap(min_fp = 0.05)            # hide trivial FP noise
#   fig_covsel_heatmap(save = TRUE)              # write PNG + PDF
# ============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# Recover the correct per-cell FP/FN denominator for a by-covar table.
#
# `n_datasets` is meant to be the number of fitted datasets in a
# (structure [x estimator x outer_opt] x sample_N) cell -- the SHARED
# denominator for true-effect FN rates (n_FN / N) and distractor FP rates
# (n_FP / N).  Some aggregators instead record it PER (var, covar, shape) row,
# so a distractor selected k times gets n_datasets = k, giving a spurious 100%
# FP rate that saturates every nonzero tile red.
#
# The authoritative total is carried by the true-effect rows and is constant
# across scenarios within a cell, so we broadcast max(n_datasets[is_true]) to
# every row of the cell (falling back to the row-wise max when a cell has no
# true-effect rows at all).  When n_datasets is already the correct constant
# (e.g. the VAE aggregator) this is a no-op, so it is safe for the runSCM, PsN
# and VAE workflows alike.
.fix_covsel_denom <- function(dat) {
  if (!all(c("n_datasets", "is_true") %in% names(dat))) return(dat)
  cell <- intersect(c("structure", "estimator", "outer_opt", "sample_N"),
                    names(dat))
  # cell-wide authoritative total (max over the cell's true-effect rows); used
  # as the fallback for the null scenario, which has NO true-effect row.
  dat <- dat |>
    dplyr::group_by(dplyr::across(dplyr::all_of(cell))) |>
    dplyr::mutate(.denom_cell = {
      tv <- suppressWarnings(max(n_datasets[is_true], na.rm = TRUE))
      if (!is.finite(tv)) suppressWarnings(max(n_datasets, na.rm = TRUE)) else tv
    }) |>
    dplyr::ungroup()
  # prefer the per-scenario total (the fitted-dataset count can vary by a few
  # datasets across scenarios, e.g. failed fits), falling back to the cell-wide
  # constant when a scenario carries no true-effect row.
  if ("scenario" %in% names(dat)) {
    dat <- dat |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c(cell, "scenario")))) |>
      dplyr::mutate(.denom_scn =
        suppressWarnings(max(n_datasets[is_true], na.rm = TRUE))) |>
      dplyr::ungroup() |>
      dplyr::mutate(.denom = dplyr::if_else(
        is.finite(.denom_scn) & .denom_scn > 0, .denom_scn, .denom_cell))
  } else {
    dat <- dplyr::mutate(dat, .denom = .denom_cell)
  }
  dat |>
    dplyr::mutate(n_datasets = dplyr::if_else(
      is.finite(.denom) & .denom > 0, as.double(.denom),
      as.double(n_datasets))) |>
    dplyr::select(-tidyselect::any_of(c(".denom", ".denom_cell", ".denom_scn")))
}

# advan4 FIRST -> leftmost structure column (analytic-first convention).
.STRUCT_LAB <- c(advan4 = "ADVAN4 (analytic)", linCmt = "linCmt (analytic)",
                 ode = "ADVAN13 (ODE)")

# canonical effect ordering: structural param first, then covariate
.VAR_LAB   <- c(cl = "CL", vc = "VC", q = "Q", vp = "VP", ka = "KA")

# TOP-DOWN covariate order on the y-axis.  The KEY message is BW<->BMI
# collinearity, so BW and BMI are placed ADJACENT: within each parameter panel a
# true "~BW" effect (outlined = FN) sits directly above its "~BMI" thief
# (= FP).  CrCL (the other continuous true effect) follows, then the
# categoricals SEX / RACE.
.COVAR_ORD <- c("BW", "BMI", "CrCL", "SEX", "RACE")

theme_scm <- function(base_size = 15) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      legend.position  = "top",
      strip.text       = ggplot2::element_text(face = "bold"),
      plot.title       = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(size = base_size - 2),
      axis.text.y      = ggplot2::element_text(size = base_size - 3)
    )
}

# ---- main ------------------------------------------------------------------
fig_covsel_heatmap <- function(
    agg_dir      = "output/vae_covsel_aggregated",
    sample_N     = NULL,       # NULL = all (40, 80, 300)
    structure    = "both",     # "both" = linCmt + ode side by side (default) | "linCmt" | "ode"
    labels       = TRUE,
    min_fp       = 0,      # blank FP cells below this rate
    layout       = c("slide", "wide"),  # "slide" = N stacked (16:9-friendly)
    save         = FALSE,
    out_dir      = "output/figures/vae_covsel",
    csv_name     = "vae_covsel_by_covar.csv",
    title_prefix = "VAE covariate selection") {

  csv <- file.path(trimws(agg_dir), csv_name)
  if (!file.exists(csv)) stop("by-covar CSV not found: ", csv)

  dat <- readr::read_csv(csv, show_col_types = FALSE)

  if (!is.null(sample_N))  dat <- dplyr::filter(dat, sample_N  %in% !!sample_N)
  # one structure per figure (default "linCmt"); "both" keeps ode + linCmt as
  # side-by-side facet columns so they share one slide.
  structure <- structure[1]
  both_struct <- identical(structure, "both")
  if (!both_struct) dat <- dplyr::filter(dat, structure == !!structure)
  if (nrow(dat) == 0L) stop("no rows after filtering (check sample_N / structure)")

  # ---- recover the correct per-cell FP/FN denominator ----------------------
  # `n_datasets` should be the number of fitted datasets in a
  # (structure[x estimator x outer_opt] x sample_N) cell -- the common
  # denominator for BOTH true-effect FN rates and distractor FP rates.  Some
  # aggregators instead store it PER (var,covar,shape) row, so a distractor
  # selected k times gets n_datasets = k and a spurious 100% FP rate (every
  # nonzero tile saturates).  The authoritative total is carried by the
  # true-effect rows and is constant across scenarios within a cell, so
  # broadcast it to every row.  When n_datasets is already the correct constant
  # (e.g. aggregate_vae_covsel.R) this is a no-op.
  dat <- .fix_covsel_denom(dat)

  # per-cell error rate + regime label.
  # y-axis is the full effect label "PARAM~COVAR" (unambiguous), ordered as a CL
  # block then a VC block; WITHIN each block BW & BMI are adjacent so the
  # collinearity pair reads vertically.
  plot_df <- dat |>
    dplyr::mutate(
      err_type = ifelse(is_true, "FN", "FP"),
      err_rate = ifelse(is_true, n_FN / n_datasets, n_FP / n_datasets),
      # blank trivial FP noise if requested
      err_show = dplyr::if_else(!is_true & err_rate < min_fp,
                                NA_real_, err_rate),
      # signed severity: FP -> positive (red), FN -> negative (blue).
      # magnitude is the error rate in %, sign encodes the regime so the two
      # error kinds get contrasting colours on ONE diverging scale.
      err_signed = dplyr::if_else(is_true, -err_show * 100, err_show * 100),
      param    = dplyr::coalesce(.VAR_LAB[var], toupper(var)),
      # effect label folds shape in so a power->lin flip reads vertically
      effect   = paste0(param, "~", covar, ".",
                        dplyr::coalesce(.SHAPE_LAB[tolower(shape)], shape)),
      v_ord    = match(var, names(.VAR_LAB)),
      c_ord    = match(covar, .COVAR_ORD),
      s_ord    = match(tolower(shape), .SHAPE_ORD),
      scenario = factor(scenario, levels = sort(unique(scenario))),
      struct_f = dplyr::recode(structure,
                               linCmt = "linCmt (analytic)", ode = "ODE"),
      sample_N = factor(paste0("N = ", sample_N),
                        levels = paste0("N = ", c(40, 80, 300)))
    ) |>
    dplyr::filter(!is.na(v_ord), !is.na(c_ord), !is.na(s_ord)) |>
    dplyr::mutate(
      struct_f = factor(struct_f,
        levels = intersect(c("linCmt (analytic)", "ODE"), unique(struct_f)))
    )

  # Complete the (sample_N x scenario x effect) grid.  A FALSE-effect distractor
  # only gets a by-covar row when it was falsely selected in >=1 dataset, so
  # never-selected distractors (genuine 0% FP, common at large N) are simply
  # absent and would leave holes in the heatmap.  Fill those as explicit white
  # zeros.  NOTE: `is_true` varies BY SCENARIO for a given (var, covar), so it
  # must NOT be part of the nesting key (that would fabricate a duplicate TRUE
  # + FALSE tile per cell).  Any cell we have to fill is by definition one that
  # was never selected -> a distractor -> is_true = FALSE, FP = 0.
  plot_df <- plot_df |>
    tidyr::complete(
      struct_f, sample_N, scenario,
      tidyr::nesting(effect, param, var, covar, shape, v_ord, c_ord, s_ord),
      fill = list(is_true = FALSE, err_rate = 0, err_show = 0, err_signed = 0)
    )

  # TOP-DOWN effect order: CL block first, then VC; BW,BMI adjacent within each;
  # within a (var,covar) the TRUE shape sits directly above its distractor shape.
  eff_levels <- plot_df |>
    dplyr::distinct(effect, v_ord, c_ord, s_ord) |>
    dplyr::arrange(v_ord, c_ord, s_ord) |>
    dplyr::pull(effect)
  # ggplot draws first factor level at BOTTOM -> reverse so CL~BW sits at top.
  plot_df <- dplyr::mutate(plot_df,
    effect = factor(effect, levels = rev(eff_levels)))

  # LAYOUT: "slide" stacks the three N panels vertically (one per row) so the
  # figure is roughly landscape and drops onto a 16:9 slide; tiles fill the
  # panel width instead of being forced square.  "wide" keeps the original
  # single-row, square-tile arrangement (good for a full-page landscape).
  layout <- match.arg(layout)
  facet_nrow <- if (layout == "slide") 3L else 1L

  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(scenario, effect, fill = err_signed)) +
    ggplot2::geom_tile(colour = "grey85", linewidth = 0.3) +
    # Bold black frame around the TRUE covariate effect of each scenario (the
    # FN regime). Drawn as a separate top layer with fill=NA so the outline is
    # never overdrawn by adjacent tiles; makes the true-vs-distractor regimes
    # visually separable at a glance.
    ggplot2::geom_tile(
      data = dplyr::filter(plot_df, is_true),
      fill = NA, colour = "black", linewidth = 0.2) +
    # diverging scale: blue = FN (true effect missed), red = FP (null selected),
    # white = no error.  Legend re-labelled so both arms read as 0-100%.
    ggplot2::scale_fill_gradient2(
      low = "#08519C", mid = "white", high = "#B30000",
      midpoint = 0, na.value = "grey92", name = NULL,
      limits = c(-100, 100), breaks = seq(-100, 100, 50),
      labels = c("FN 100%", "FN 50%", "0", "FP 50%", "FP 100%")) +
    ggplot2::labs(
      x = "Simulation scenario", y = "Covariate effect (param ~ covariate . shape)"
    ) +
    theme_scm() +
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = 14))

  # FACETING: with both structures, use a N(rows) x structure(cols) grid so ode
  # and linCmt share one slide; otherwise keep the single-structure N layout.
  if (both_struct) {
    p <- p + ggplot2::facet_grid(sample_N ~ struct_f)
  } else {
    p <- p + ggplot2::facet_wrap(~ sample_N, nrow = facet_nrow)
  }

  # square tiles only in the "wide" layout; "slide" lets tiles fill the width
  if (layout == "wide") p <- p + ggplot2::coord_equal()

  if (isTRUE(labels)) {
    lab_df <- dplyr::filter(plot_df, !is.na(err_show))
    p <- p + ggplot2::geom_text(
      data = lab_df,
      ggplot2::aes(label = round(err_show * 100)),
      size = if (both_struct) 4.0 else 3.0,
      colour = ifelse(lab_df$err_show > 0.55, "white", "grey20"))
  }

  if (isTRUE(save)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    n_tag <- if (is.null(sample_N)) "allN" else paste0("N", paste(sample_N, collapse = "-"))
    s_tag <- structure
    stub  <- file.path(out_dir, sprintf("fig_covsel_heatmap_%s_%s", n_tag, s_tag))
    dims  <- if (both_struct) c(12, 10) else if (layout == "slide") c(13, 9) else c(22, 12)
    ggplot2::ggsave(paste0(stub, ".png"), p, width = dims[1], height = dims[2], dpi = 200)
    ggplot2::ggsave(paste0(stub, ".pdf"), p, width = dims[1], height = dims[2])
    message("saved: ", stub, ".{png,pdf}")
  }

  p
}

# ============================================================================
# fig_covsel_heatmap_scm()  --  SCM covariate-selection error-pattern heatmap
# ----------------------------------------------------------------------------
# SCM differs from VAE in TWO ways this function accommodates:
#
#   (1) SHAPE is a selectable dimension.  Each (var, covar) has a TRUE shape
#       (continuous BW/BMI/CrCL -> "power"; categorical SEX/RACE -> "cat") plus
#       DISTRACTOR shapes (e.g. "lin" for a continuous covariate).  A true
#       covariate can therefore be recovered in the WRONG shape -- the
#       lin-vs-power shape-flip that drives the strict-power collapse.  We fold
#       shape into the y-axis effect label ("CL~BW.power", "CL~BW.lin", ...) and
#       place the distractor shape DIRECTLY BELOW its true shape, so a shape
#       flip reads vertically: a dark FN tile on "CL~BW.power" stacked on a dark
#       FP tile on "CL~BW.lin" is the signature of a flip.
#
#   (2) There are ESTIMATOR x OUTER_OPT cells.  The by-covar CSV carries both
#       columns, so a single figure MUST fix one estimator/outer_opt cell
#       (defaults focei/bobyqa).  A guard errors if >1 cell survives filtering.
#
# Reads output/scm_bench_rescue_winner_aggregated/scm_covsel_by_covar.csv
#   cols: sample_N, scenario, structure, estimator, outer_opt, var, covar,
#         shape, n_datasets, is_true, n_detected, detection_rate,
#         n_TP, n_FN, n_FP
# Denominator: `n_datasets` in the by-covar CSV is the cell's true fit count
# (derived in aggregate_scm_estimator2.1.R), so FP/FN rates divide by it
# directly -- no scm_diag_rates.csv join needed.
#
# USAGE:
#   source("script/viz/fig_covsel_heatmap.R")
#   fig_covsel_heatmap_scm()                                  # focei/bobyqa, linCmt, all N
#   fig_covsel_heatmap_scm(sample_N = 300)
#   fig_covsel_heatmap_scm(estimator = "focei", outer_opt = "bobyqa",
#                          structure = "ode")
#   fig_covsel_heatmap_scm(labels = FALSE, save = TRUE)
# ============================================================================

# within-pair shape order: TRUE shape first (top), distractors below.
# `exp` is included for PsN/NONMEM SCM, which tests BOTH an exponential and a
# power form for each continuous covariate ([code] section); truth is `power`,
# so a selected `exp` is a legitimate shape-mismatch distractor.
.SHAPE_ORD <- c("power", "exp", "lin", "cat")
.SHAPE_LAB <- c(power = "power", exp = "exp", lin = "lin", cat = "cat")

fig_covsel_heatmap_scm <- function(
    agg_dir      = "output/scm_bench_rescue_winner_aggregated",
    sample_N     = NULL,        # NULL = all (40, 80, 300)
    structure    = "linCmt",    # ONE structure per figure ("linCmt" | "ode")
    estimator    = "focei",     # fix ONE estimator cell
    outer_opt    = "bobyqa",    # fix ONE outer optimiser cell
    labels       = TRUE,
    min_fp       = 0,           # blank FP cells below this rate
    layout       = c("slide", "wide"),
    show_title   = TRUE,        # FALSE = drop the plot title (tighter figure)
    save         = FALSE,
    out_dir      = "output/figures/scm_covsel",
    csv_name     = "scm_covsel_by_covar.csv",
    title_prefix = "SCM covariate selection") {

  csv <- file.path(trimws(agg_dir), csv_name)
  if (!file.exists(csv)) stop("by-covar CSV not found: ", csv)
  dat <- readr::read_csv(csv, show_col_types = FALSE)

  # ---- fix the estimator x outer_opt cell (SCM-specific) -------------------
  if (!is.null(estimator) && "estimator" %in% names(dat))
    dat <- dplyr::filter(dat, estimator == !!estimator)
  if (!is.null(outer_opt) && "outer_opt" %in% names(dat))
    dat <- dplyr::filter(dat, outer_opt == !!outer_opt)

  if (!is.null(sample_N)) dat <- dplyr::filter(dat, sample_N %in% !!sample_N)
  structure <- structure[1]
  both_struct <- identical(structure, "both")
  if (!both_struct) dat <- dplyr::filter(dat, structure == !!structure)
  if (nrow(dat) == 0L) stop("no rows after filtering (check sample_N / structure / estimator / outer_opt)")

  # guard: exactly ONE estimator x outer_opt cell must remain, else the fill
  # would silently blend distinct estimators onto one tile.
  cells <- dat |> dplyr::distinct(estimator, outer_opt)
  if (nrow(cells) > 1L) {
    msg <- paste(sprintf("%s/%s", cells$estimator, cells$outer_opt), collapse = ", ")
    stop("more than one estimator/outer_opt cell after filtering: ", msg,
         "\n  -> set estimator= and outer_opt= to pick exactly one.")
  }
  est_tag <- paste(cells$estimator[1], cells$outer_opt[1], sep = "_")

  # ---- recover the correct per-cell FP/FN denominator ----------------------
  # See .fix_covsel_denom(): some aggregators store n_datasets PER covariate row
  # (distractor selected k times -> n_datasets = k -> spurious 100% FP rate).
  # The true per-cell total lives on the true-effect rows and is constant across
  # scenarios, so broadcast it.  No-op when n_datasets is already correct.
  dat <- .fix_covsel_denom(dat)

  plot_df <- dat |>
    dplyr::mutate(
      err_type = ifelse(is_true, "FN", "FP"),
      err_rate = ifelse(is_true, n_FN / n_datasets, n_FP / n_datasets),
      err_show = dplyr::if_else(!is_true & err_rate < min_fp, NA_real_, err_rate),
      err_signed = dplyr::if_else(is_true, -err_show * 100, err_show * 100),
      # var/covar casing varies by source (runSCM lowercase 'cl'/'CrCL' vs PsN
      # canonicalised 'CL'/'CRCL'); match case-insensitively so both render.
      param    = dplyr::coalesce(.VAR_LAB[tolower(var)], toupper(var)),
      # effect label folds shape in so flips read vertically
      effect   = paste0(param, "~", toupper(covar), ".",
                        .SHAPE_LAB[tolower(shape)] %||% shape),
      v_ord    = match(tolower(var), names(.VAR_LAB)),
      c_ord    = match(toupper(covar), toupper(.COVAR_ORD)),
      s_ord    = match(tolower(shape), .SHAPE_ORD),
      scenario = factor(scenario, levels = sort(unique(scenario))),
      struct_f = dplyr::recode(structure,
                               linCmt = "linCmt (analytic)", ode = "ADVAN13 (ODE)",
                               advan4 = "ADVAN4 (analytic)"),
      sample_N = factor(paste0("N = ", sample_N),
                        levels = paste0("N = ", c(40, 80, 300)))
    ) |>
    dplyr::filter(!is.na(v_ord), !is.na(c_ord), !is.na(s_ord)) |>
    dplyr::mutate(
      struct_f = factor(struct_f,
        levels = intersect(c("ADVAN4 (analytic)", "linCmt (analytic)", "ADVAN13 (ODE)"),
                           unique(struct_f)))
    )

  # drop unused N / structure levels so tidyr::complete() below does not
  # re-fabricate empty panels for filtered-out sample sizes or structures.
  plot_df <- droplevels(plot_df)

  # complete the (sample_N x scenario x effect) grid.  As with VAE, a distractor
  # (var,covar,shape) only gets a row when falsely selected >=1 time, so
  # never-selected distractors are absent and must be filled as explicit FP=0.
  # is_true varies BY SCENARIO (and by shape), so it is NOT part of the nesting
  # key -- any filled cell is by definition a never-selected distractor
  # (is_true = FALSE, FP = 0).
  plot_df <- plot_df |>
    tidyr::complete(
      struct_f, sample_N, scenario,
      tidyr::nesting(effect, param, var, covar, shape, v_ord, c_ord, s_ord),
      fill = list(is_true = FALSE, err_rate = 0, err_show = 0, err_signed = 0,
                  n_FN = 0, n_FP = 0)
    ) |>
    # per-tile error COUNT: FN count on true cells, FP count on distractor cells
    dplyr::mutate(err_count = ifelse(is_true, n_FN, n_FP))

  # TOP-DOWN order: CL block then VC; within each, covar order, then shape order
  # (true shape first so its distractor sits directly below it).
  eff_levels <- plot_df |>
    dplyr::distinct(effect, v_ord, c_ord, s_ord) |>
    dplyr::arrange(v_ord, c_ord, s_ord) |>
    dplyr::pull(effect)
  plot_df <- dplyr::mutate(plot_df,
    effect = factor(effect, levels = rev(eff_levels)))

  layout <- match.arg(layout)
  facet_nrow <- if (layout == "slide") 3L else 1L

  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(scenario, effect, fill = err_signed)) +
    ggplot2::geom_tile(colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_tile(
      data = dplyr::filter(plot_df, is_true),
      fill = NA, colour = "black", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(
      low = "#08519C", mid = "white", high = "#B30000",
      midpoint = 0, na.value = "grey92", name = NULL,
      limits = c(-100, 100), breaks = seq(-100, 100, 50),
      labels = c("FN 100%", "FN 50%", "0", "FP 50%", "FP 100%")) +
    ggplot2::labs(
      title = if (isTRUE(show_title))
                sprintf("%s  (%s, %s)", title_prefix, est_tag, structure) else NULL,
      x = "Simulation scenario",
      y = "Covariate effect (param ~ covariate . shape)"
    ) +
    theme_scm() +
    ggplot2::theme(
      panel.spacing = grid::unit(0.4, "lines"),
      plot.margin   = grid::unit(c(2, 2, 2, 2), "pt"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.spacing = grid::unit(2, "pt")
    ) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = 14))

  if (both_struct) {
    p <- p + ggplot2::facet_grid(sample_N ~ struct_f)
  } else {
    p <- p + ggplot2::facet_wrap(~ sample_N, nrow = facet_nrow)
  }
  # Pin tiles square in BOTH layouts so they cannot stretch to fill a wide plot
  # pane (the "too wide" symptom).  In "slide" the 3 N panels stack vertically
  # with square tiles -> a compact block that drops onto a 16:9 slide; in "wide"
  # they sit side-by-side.
  p <- p + ggplot2::coord_equal()

  if (isTRUE(labels)) {
    # print the raw FN/FP dataset COUNT in each tile.  Never-selected covariate
    # relationships were filled to 0 by tidyr::complete() above, so (like VAE)
    # they show an explicit "0" rather than a blank tile.
    lab_df <- dplyr::filter(plot_df, !is.na(err_count))
    p <- p + ggplot2::geom_text(
      data = lab_df,
      ggplot2::aes(label = err_count),
      size = if (both_struct) 3.2 else 3.0,
      fontface = "bold",
      colour = ifelse(lab_df$err_show > 0.55, "white", "grey15"))
  }

  if (isTRUE(save)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    n_tag <- if (is.null(sample_N)) "allN" else paste0("N", paste(sample_N, collapse = "-"))
    stub  <- file.path(out_dir,
      sprintf("fig_covsel_heatmap_scm_%s_%s_%s", est_tag, n_tag, structure))
    dims  <- if (both_struct) c(14, 15) else if (layout == "slide") c(9, 15) else c(22, 14)
    ggplot2::ggsave(paste0(stub, ".png"), p, width = dims[1], height = dims[2], dpi = 200)
    ggplot2::ggsave(paste0(stub, ".pdf"), p, width = dims[1], height = dims[2])
    message("saved: ", stub, ".{png,pdf}")
  }

  p
}
