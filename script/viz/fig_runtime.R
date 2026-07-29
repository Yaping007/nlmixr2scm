# ============================================================================
# fig_runtime.R  --  wall-time comparison across estimator / optimiser combos
# ----------------------------------------------------------------------------
# Answers: "how expensive is each estimator x outer-optimiser combination, per
#           fit and in TOTAL compute, and how does that scale with sample size?"
#
# Reads the PER-DATASET diagnostic long tables (one row = one fitted dataset):
#   VAE : <agg>/vae_diag_long.csv   wall col = fit_runtime_sec
#         (single fit; estimator = "vae", outer_opt = NA)
#   SCM : <agg>/scm_diag_long.csv   wall col = wall_total_sec
#         (= wall_base + wall_scm + wall_refit; carries estimator/outer_opt)
#
# The two schemas differ ONLY in the wall-time column name and in whether the
# estimator/outer_opt columns vary, so this harmonises them onto a common
# tibble: (combo, estimator, outer_opt, sample_N, scenario, structure,
#          dataset_id, wall_sec) and stacks any number of sources.
#
# IMPORTANT comparability caveat.  SCM wall_total_sec covers the WHOLE pipeline
# (base fit + stepwise search + tight-tol refit), whereas VAE fit_runtime_sec is
# a SINGLE likelihood-evaluation fit with no search.  They are therefore NOT the
# same workload -- VAE will look far cheaper.  To compare only the base fit set
# the SCM source's wall = "wall_base_sec".
#
# OUTPUTS
#   * plot  : per-fit wall-time distribution (box or violin) by sample size,
#             coloured by combo, faceted by structure.  y in minutes (default)
#             on an optional log scale (runtimes span orders of magnitude).
#   * table : returned as attr(p, "summary") -- per (combo, sample_N, structure)
#             n, median, mean, p90 wall AND total compute (hours = sum of all
#             fits / 3600).  return_data = TRUE returns list(tidy, summary).
#
# KNOBS
#   sources     : NAMED list; each element list(agg_dir=, csv=, wall=, estimator=,
#                 outer_opt=).  `wall` auto-detects if omitted (wall_total_sec
#                 else fit_runtime_sec).  estimator/outer_opt optionally subset a
#                 multi-estimator SCM file.  Default = VAE + focei/bobyqa +
#                 irlsfocei/bobyqa from the rescue-winner tree.
#   sample_N    : NULL = all | subset e.g. c(80, 300)
#   structure   : NULL = all | "linCmt" | "ode"
#   scenario    : NULL = all | subset e.g. c(1, 8, 16)
#   by_scenario : TRUE (default) resolve runtime PER scenario | FALSE pool
#                 across scenarios (x = sample_N, facet structure)
#   style       : "trend" (default) median-vs-scenario line + p90 band, facet
#                 sample_N x structure -- scan all 16 scenarios horizontally |
#                 "dist" per-fit box/violin distribution
#   unit        : "min" (default) | "sec" | "hour"
#   geom        : "box" (default) | "violin"
#   log_y       : TRUE (default) log10 wall axis | FALSE linear
#   return_data : FALSE = ggplot (summary in attr) | TRUE = list(tidy, summary)
#   save        : FALSE | TRUE (ggsave PNG + PDF + _summary.csv)
#
# USAGE:
#   source("script/viz/fig_runtime.R")
#   p   <- fig_runtime()                       # PER-scenario grid (16 rows)
#   attr(p, "summary")                         # per (combo,scenario,N,struct)
#   fig_runtime(by_scenario = FALSE)           # old pooled-over-scenarios view
#   fig_runtime(scenario = c(1, 8, 16))        # just a few scenarios
#   fig_runtime(structure = "linCmt", geom = "violin")
#   tab <- fig_runtime(return_data = TRUE)$summary
#   fig_runtime(save = TRUE)
# ============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(purrr)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

.STRUCT_LAB <- c(linCmt = "linCmt (analytic)", ode = "ODE")

# qualitative palette for combos (extended lazily if >6 sources supplied)
.PAL_COMBO <- c("#7F7F7F", "#1F77B4", "#D62728", "#2CA02C", "#9467BD", "#E8820C")

.UNIT_DIV <- c(sec = 1, min = 60, hour = 3600)

.theme_rt <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "top",
      strip.text       = ggplot2::element_text(face = "bold"),
      plot.title       = ggplot2::element_text(face = "bold")
    )
}

# ---- load + harmonise one source -------------------------------------------
.load_runtime_source <- function(label, spec, sample_N, structure) {
  csv <- file.path(trimws(spec$agg_dir), spec$csv %||%
                     stop("source '", label, "' needs a csv= name"))
  if (!file.exists(csv)) stop("runtime CSV not found for '", label, "': ", csv)
  d <- readr::read_csv(csv, show_col_types = FALSE)

  # wall column: explicit, else auto (SCM total, else VAE single fit)
  wall_col <- spec$wall %||%
    (if ("wall_total_sec" %in% names(d)) "wall_total_sec"
     else if ("fit_runtime_sec" %in% names(d)) "fit_runtime_sec"
     else stop("no wall-time column in '", label, "' (", csv,
               "); set wall= explicitly"))
  if (!wall_col %in% names(d))
    stop("wall column '", wall_col, "' absent in '", label, "'")

  # optional estimator/outer_opt subset (multi-estimator SCM files)
  if (!is.null(spec$estimator) && "estimator" %in% names(d))
    d <- dplyr::filter(d, estimator %in% spec$estimator)
  if (!is.null(spec$outer_opt) && "outer_opt" %in% names(d))
    d <- dplyr::filter(d, outer_opt %in% spec$outer_opt)
  if (!is.null(sample_N))  d <- dplyr::filter(d, sample_N %in% !!sample_N)
  if (!is.null(structure)) d <- dplyr::filter(d, structure %in% !!structure)
  if (nrow(d) == 0L)
    stop("source '", label, "' has no rows after filtering (check estimator / ",
         "outer_opt / sample_N / structure)")

  d |>
    dplyr::transmute(
      combo      = label,
      estimator  = if ("estimator" %in% names(d)) estimator else NA_character_,
      outer_opt  = if ("outer_opt" %in% names(d)) outer_opt else NA_character_,
      sample_N, scenario, structure,
      dataset_id = if ("dataset_id" %in% names(d)) dataset_id else NA_integer_,
      wall_sec   = .data[[wall_col]]
    ) |>
    dplyr::filter(is.finite(wall_sec))
}

# ---- main ------------------------------------------------------------------
fig_runtime <- function(
    sources = list(
      "VAE" = list(
        agg_dir = "output/vae_covsel_aggregated",
        csv     = "vae_diag_long.csv"),
      "FOCEi / bobyqa" = list(
        agg_dir = "output/scm_bench_rescue_winner_aggregated",
        csv     = "scm_diag_long.csv",
        estimator = "focei", outer_opt = "bobyqa"),
      "IRLS-FOCEi / bobyqa" = list(
        agg_dir = "output/scm_bench_rescue_winner_aggregated",
        csv     = "scm_diag_long.csv",
        estimator = "irlsfocei", outer_opt = "bobyqa")
    ),
    sample_N    = NULL,
    structure   = NULL,
    scenario    = NULL,
    by_scenario = TRUE,
    style       = c("trend", "dist"),
    unit        = c("min", "sec", "hour"),
    geom        = c("box", "violin"),
    log_y       = TRUE,
    return_data = FALSE,
    save        = FALSE,
    out_dir     = "output/figures/scm_covsel",
    title       = "Wall-time by estimator / optimiser") {

  unit <- match.arg(unit)
  geom <- match.arg(geom)
  style <- match.arg(style)
  if (is.null(names(sources)) || any(names(sources) == ""))
    stop("`sources` must be a NAMED list (names become the combo labels)")

  # stack all sources, keeping label order for the x/legend
  tidy <- purrr::map2_dfr(names(sources), sources,
                          ~ .load_runtime_source(.x, .y, sample_N, structure))
  if (!is.null(scenario)) {
    tidy <- dplyr::filter(tidy, scenario %in% !!scenario)
    if (nrow(tidy) == 0L)
      stop("no rows after scenario subset: ", paste(scenario, collapse = ", "))
  }
  combo_levels <- names(sources)[names(sources) %in% unique(tidy$combo)]

  div <- .UNIT_DIV[[unit]]
  tidy <- tidy |>
    dplyr::mutate(
      wall      = wall_sec / div,
      combo     = factor(combo, levels = combo_levels),
      sample_N  = factor(sample_N, levels = sort(unique(sample_N))),
      scenario  = factor(scenario, levels = sort(unique(scenario))),
      structure = factor(structure, levels = names(.STRUCT_LAB))
    )

  # ---- summary: per-cell stats + TOTAL compute (hours) ---------------------
  grp <- if (isTRUE(by_scenario))
    dplyr::vars(combo, scenario, sample_N, structure)
  else
    dplyr::vars(combo, sample_N, structure)

  summary_tbl <- tidy |>
    dplyr::group_by(!!!grp) |>
    dplyr::summarise(
      n_fits       = dplyr::n(),
      med_wall     = stats::median(wall),
      mean_wall    = mean(wall),
      p90_wall     = stats::quantile(wall, 0.90, names = FALSE),
      total_hours  = sum(wall_sec) / 3600,
      .groups = "drop"
    ) |>
    dplyr::rename_with(~ paste0(., "_", unit),
                       c(med_wall, mean_wall, p90_wall)) |>
    dplyr::arrange(combo, structure, sample_N)

  if (isTRUE(return_data))
    return(list(tidy = tidy, summary = summary_tbl))

  # ---- plot ----------------------------------------------------------------
  pal <- setNames(.PAL_COMBO[seq_along(combo_levels)], combo_levels)
  y_lab <- sprintf("Wall time per fit (%s)%s",
                   unit, if (log_y) ", log scale" else "")

  if (style == "trend" && isTRUE(by_scenario)) {
    # scenario on x-axis: median line + p90 ribbon, facet sample_N x structure
    med_col <- paste0("med_wall_", unit)
    p90_col <- paste0("p90_wall_", unit)
    plt <- summary_tbl |>
      dplyr::mutate(
        scen_i = as.integer(as.character(scenario)),
        med    = .data[[med_col]],
        p90    = .data[[p90_col]])

    p <- ggplot2::ggplot(
        plt, ggplot2::aes(x = scen_i, y = med,
                          colour = combo, fill = combo, group = combo)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = med, ymax = p90),
                           alpha = 0.12, colour = NA) +
      ggplot2::geom_line(linewidth = 0.6) +
      ggplot2::geom_point(size = 1.4) +
      ggplot2::facet_grid(sample_N ~ structure,
                          labeller = ggplot2::labeller(
                            structure = .STRUCT_LAB,
                            sample_N  = function(x) paste0("N = ", x)),
                          scales = "free_y") +
      ggplot2::scale_x_continuous(
        breaks = sort(unique(plt$scen_i)),
        expand = ggplot2::expansion(mult = 0.02)) +
      ggplot2::scale_fill_manual(values = pal, name = NULL) +
      ggplot2::scale_colour_manual(values = pal, name = NULL) +
      ggplot2::labs(
        title = title, x = "Scenario",
        y = sprintf("Median wall time per fit (%s)%s (band \u2192 p90)",
                    unit, if (log_y) ", log scale" else "")) +
      .theme_rt()
    if (log_y)
      p <- p + ggplot2::scale_y_log10(labels = scales::label_number())

    attr(p, "summary") <- summary_tbl

    if (isTRUE(save)) {
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      combos_tag <- gsub("[^A-Za-z0-9]+", "", paste(combo_levels, collapse = "-"))
      stub <- file.path(out_dir,
                        sprintf("fig_runtime_%s_%s_trend", unit, combos_tag))
      ggplot2::ggsave(paste0(stub, ".png"), p, width = 11, height = 8,
                      dpi = 150, limitsize = FALSE)
      ggplot2::ggsave(paste0(stub, ".pdf"), p, width = 11, height = 8,
                      limitsize = FALSE)
      readr::write_csv(summary_tbl, paste0(stub, "_summary.csv"))
      message("saved: ", stub, ".{png,pdf} (+ _summary.csv)")
    }
    return(p)
  }

  # ---- plot: per-fit wall distribution (style = 'dist') --------------------
  p <- ggplot2::ggplot(
    tidy, ggplot2::aes(x = sample_N, y = wall, fill = combo, colour = combo))

  if (geom == "box") {
    p <- p + ggplot2::geom_boxplot(
      position = ggplot2::position_dodge(width = 0.8),
      alpha = 0.35, outlier.size = 0.5, outlier.alpha = 0.25,
      linewidth = 0.4)
  } else {
    p <- p + ggplot2::geom_violin(
      position = ggplot2::position_dodge(width = 0.8),
      alpha = 0.30, linewidth = 0.4, scale = "width")
  }

  p <- p +
    (if (isTRUE(by_scenario))
       ggplot2::facet_grid(scenario ~ structure,
                           labeller = ggplot2::labeller(structure = .STRUCT_LAB))
     else
       ggplot2::facet_wrap(~ structure, ncol = 1,
                          labeller = ggplot2::labeller(structure = .STRUCT_LAB))) +
    ggplot2::scale_fill_manual(values = pal, name = NULL) +
    ggplot2::scale_colour_manual(values = pal, name = NULL) +
    ggplot2::labs(title = title, x = "Sample size (subjects)", y = y_lab) +
    .theme_rt()

  if (log_y)
    p <- p + ggplot2::scale_y_log10(labels = scales::label_number())

  # attach the summary so it travels with the returned plot
  attr(p, "summary") <- summary_tbl

  if (isTRUE(save)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    combos_tag <- gsub("[^A-Za-z0-9]+", "", paste(combo_levels, collapse = "-"))
    scn_tag <- if (isTRUE(by_scenario)) "_byScn" else ""
    stub <- file.path(out_dir,
                      sprintf("fig_runtime_%s_%s%s", unit, combos_tag, scn_tag))
    n_row <- if (isTRUE(by_scenario)) nlevels(tidy$scenario) else 2L
    h     <- if (isTRUE(by_scenario)) max(8, 1.3 * n_row) else 8
    ggplot2::ggsave(paste0(stub, ".png"), p, width = 9, height = h,
                    dpi = 150, limitsize = FALSE)
    ggplot2::ggsave(paste0(stub, ".pdf"), p, width = 9, height = h,
                    limitsize = FALSE)
    readr::write_csv(summary_tbl, paste0(stub, "_summary.csv"))
    message("saved: ", stub, ".{png,pdf} (+ _summary.csv)")
  }

  p
}
