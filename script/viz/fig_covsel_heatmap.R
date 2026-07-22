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

theme_scm <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      legend.position  = "top",
      strip.text       = ggplot2::element_text(face = "bold"),
      plot.title       = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(size = base_size - 3),
      axis.text.y      = ggplot2::element_text(size = base_size - 3)
    )
}

# ---- main ------------------------------------------------------------------
fig_covsel_heatmap <- function(
    agg_dir      = "output/vae_covsel_aggregated",
    sample_N     = NULL,       # NULL = all (40, 80, 300)
    structure    = "linCmt",   # ONE structure per figure ("linCmt" | "ode")
    labels       = TRUE,
    min_fp       = 0,      # blank FP cells below this rate
    save         = FALSE,
    out_dir      = "output/figures/vae_covsel",
    csv_name     = "vae_covsel_by_covar.csv",
    title_prefix = "VAE covariate selection") {

  csv <- file.path(trimws(agg_dir), csv_name)
  if (!file.exists(csv)) stop("by-covar CSV not found: ", csv)

  dat <- readr::read_csv(csv, show_col_types = FALSE)

  if (!is.null(sample_N))  dat <- dplyr::filter(dat, sample_N  %in% !!sample_N)
  # one structure per figure (default "linCmt"); keep first if several passed
  structure <- structure[1]
  dat <- dplyr::filter(dat, structure == !!structure)
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
      param    = dplyr::coalesce(.VAR_LAB[var], toupper(var)),
      effect   = paste0(param, "~", covar),
      v_ord    = match(var, names(.VAR_LAB)),
      c_ord    = match(covar, .COVAR_ORD),
      scenario = factor(scenario, levels = sort(unique(scenario))),
      sample_N = factor(paste0("N = ", sample_N),
                        levels = paste0("N = ", c(40, 80, 300)))
    ) |>
    dplyr::filter(!is.na(v_ord), !is.na(c_ord))

  # TOP-DOWN effect order: CL block first, then VC; BW,BMI adjacent within each.
  eff_levels <- plot_df |>
    dplyr::distinct(effect, v_ord, c_ord) |>
    dplyr::arrange(v_ord, c_ord) |>
    dplyr::pull(effect)
  # ggplot draws first factor level at BOTTOM -> reverse so CL~BW sits at top.
  plot_df <- dplyr::mutate(plot_df,
    effect = factor(effect, levels = rev(eff_levels)))

  # bold outline only on TRUE-effect (FN) cells
  true_df <- dplyr::filter(plot_df, is_true)

  struct_lab <- .STRUCT_LAB[[structure]] %||% structure

  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(scenario, effect, fill = err_show * 100)) +
    ggplot2::geom_tile(colour = "grey85", linewidth = 0.3) +
    # re-draw TRUE cells with a bold black border to flag the FN regime
    ggplot2::geom_tile(data = true_df, colour = "black", linewidth = 0.8,
                       fill = NA) +
    ggplot2::scale_fill_gradient(
      low = "#FFF5EB", high = "#B30000", na.value = "grey92",
      name = "Error rate", limits = c(0, 100),
      breaks = seq(0, 100, 25), labels = function(x) paste0(x, "%")) +
    ggplot2::facet_wrap(~ sample_N, nrow = 1) +
    ggplot2::labs(
      title    = sprintf("%s: where do FP / FN errors concentrate?",
                         title_prefix),
      subtitle = paste0(struct_lab,
                        " | Bold outline = TRUE effect (FN rate).  ",
                        "Plain tile = null effect (FP rate).  ",
                        "BW & BMI adjacent to expose their confusion."),
      x = "Simulation scenario", y = "Covariate effect (param ~ covariate)"
    ) +
    theme_scm() +
    ggplot2::coord_equal()

  if (isTRUE(labels)) {
    lab_df <- dplyr::filter(plot_df, !is.na(err_show))
    p <- p + ggplot2::geom_text(
      data = lab_df,
      ggplot2::aes(label = round(err_show * 100)),
      size = 2.6,
      colour = ifelse(lab_df$err_show > 0.55, "white", "grey20"))
  }

  if (isTRUE(save)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    n_tag <- if (is.null(sample_N)) "allN" else paste0("N", paste(sample_N, collapse = "-"))
    s_tag <- structure
    stub  <- file.path(out_dir, sprintf("fig_covsel_heatmap_%s_%s", n_tag, s_tag))
    ggplot2::ggsave(paste0(stub, ".png"), p, width = 11, height = 6, dpi = 150)
    ggplot2::ggsave(paste0(stub, ".pdf"), p, width = 11, height = 6)
    message("saved: ", stub, ".{png,pdf}")
  }

  p
}
