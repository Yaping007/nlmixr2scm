# ============================================================================
# fig_covsel_heatmap.R  --  VAE covariate-selection error-pattern heatmap
# ----------------------------------------------------------------------------
# Answers: "under each scenario, which covariate EFFECTS tend to become
#           false positive (FP) or false negative (FN)?"
#
# Reads the per-covariate roll-up written by aggregate_vae_covsel.R:
#   output/vae_covsel_aggregated/vae_covsel_by_covar.csv
#   cols: sample_N, scenario, structure, var, covar, n_datasets,
#         is_true, n_detected, detection_rate, n_TP, n_FN, n_FP
#
# Per (scenario x effect) cell there are two mutually-exclusive regimes:
#   * TRUE effect  (is_true = TRUE) : error = FALSE NEGATIVE
#                                     FN rate = n_FN / N_cell
#   * FALSE effect (is_true = FALSE): error = FALSE POSITIVE
#                                     FP rate = n_FP / N_cell
#
# CRITICAL denominator note.  In vae_covsel_by_covar.csv the `n_datasets` column
# is NOT the cell's dataset count -- it is the number of datasets in which THAT
# (var, covar) pair appeared in covsel_long.  A distractor only appears when it
# is falsely selected, so for FALSE effects n_datasets == n_FP and n_FP/n_datasets
# is ALWAYS ~1 (a bug if used as an FP rate).  We therefore divide n_FP and n_FN
# by the true per-cell dataset total N_cell (from vae_diag_rates.csv n_total /
# equivalently vae_power.csv N), NOT by n_datasets.  True effects appear in every
# dataset so their n_datasets already equals N_cell, but we use N_cell for both
# regimes for correctness and symmetry.
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

.STRUCT_LAB <- c(linCmt = "linCmt (analytic)", ode = "ODE")

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
      axis.text.x      = ggplot2::element_text(size = base_size - 1),
      axis.text.y      = ggplot2::element_text(size = base_size - 1)
    )
}

# ---- main ------------------------------------------------------------------
fig_covsel_heatmap <- function(
    agg_dir      = "output/vae_covsel_aggregated",
    sample_N     = NULL,       # NULL = all (40, 80, 300)
    structure    = "linCmt",   # ONE structure per figure ("linCmt" | "ode")
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

  # per-cell dataset total N_cell -- the correct FP/FN denominator (see header).
  # n_datasets in the by-covar CSV is per-pair appearance count, not the cell N.
  diag_csv <- file.path(trimws(agg_dir), "vae_diag_rates.csv")
  if (!file.exists(diag_csv)) stop("diag-rates CSV not found (needed for N_cell): ", diag_csv)
  n_cell <- readr::read_csv(diag_csv, show_col_types = FALSE) |>
    dplyr::select(sample_N, scenario, structure, N_cell = n_total)

  # per-cell error rate + regime label.
  # y-axis is the full effect label "PARAM~COVAR" (unambiguous), ordered as a CL
  # block then a VC block; WITHIN each block BW & BMI are adjacent so the
  # collinearity pair reads vertically.
  plot_df <- dat |>
    dplyr::left_join(n_cell, by = c("sample_N", "scenario", "structure")) |>
    dplyr::mutate(
      err_type = ifelse(is_true, "FN", "FP"),
      err_rate = ifelse(is_true, n_FN / N_cell, n_FP / N_cell),
      # blank trivial FP noise if requested
      err_show = dplyr::if_else(!is_true & err_rate < min_fp,
                                NA_real_, err_rate),
      # signed severity: FP -> positive (red), FN -> negative (blue).
      # magnitude is the error rate in %, sign encodes the regime so the two
      # error kinds get contrasting colours on ONE diverging scale.
      err_signed = dplyr::if_else(is_true, -err_show * 100, err_show * 100),
      param    = dplyr::coalesce(.VAR_LAB[var], toupper(var)),
      effect   = paste0(param, "~", covar),
      v_ord    = match(var, names(.VAR_LAB)),
      c_ord    = match(covar, .COVAR_ORD),
      scenario = factor(scenario, levels = sort(unique(scenario))),
      struct_f = dplyr::recode(structure,
                               linCmt = "linCmt (analytic)", ode = "ODE"),
      sample_N = factor(paste0("N = ", sample_N),
                        levels = paste0("N = ", c(40, 80, 300)))
    ) |>
    dplyr::filter(!is.na(v_ord), !is.na(c_ord)) |>
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
      tidyr::nesting(effect, param, var, covar, v_ord, c_ord),
      fill = list(is_true = FALSE, err_rate = 0, err_show = 0, err_signed = 0)
    )

  # TOP-DOWN effect order: CL block first, then VC; BW,BMI adjacent within each.
  eff_levels <- plot_df |>
    dplyr::distinct(effect, v_ord, c_ord) |>
    dplyr::arrange(v_ord, c_ord) |>
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
    # diverging scale: blue = FN (true effect missed), red = FP (null selected),
    # white = no error.  Legend re-labelled so both arms read as 0-100%.
    ggplot2::scale_fill_gradient2(
      low = "#08519C", mid = "white", high = "#B30000",
      midpoint = 0, na.value = "grey92", name = NULL,
      limits = c(-100, 100), breaks = seq(-100, 100, 50),
      labels = c("FN 100%", "FN 50%", "0", "FP 50%", "FP 100%")) +
    ggplot2::labs(
      x = "Simulation scenario", y = "Covariate effect (param ~ covariate)"
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
      size = if (both_struct) 4.5 else 3.4,
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
