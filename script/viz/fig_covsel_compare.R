## fig_covsel_compare.R
## Cross-method comparison of covariate-selection PATTERN (shape-resolved).
##
## The per-method heatmap (fig_covsel_heatmap_scm) is information-dense; to
## compare THREE methods we plot the comparison directly instead of three raw
## heatmaps:
##
##   fig_covsel_diff()        signed DIFFERENCE heatmap, method - reference
##                            (default reference = PsN-SCM).  White = matches
##                            PsN, red = method worse (more error), blue =
##                            method better.  Shape rows are KEPT (not
##                            collapsed): power / cat align across methods;
##                            each method's continuous distractor (nlmixr2 =
##                            "lin", PsN = "exp") has no counterpart in the
##                            other vocabulary and is shown as grey (no diff).
##
##   fig_covsel_agreement()   per-scenario scalar agreement vs the reference:
##                            100 - mean |error-rate difference| over all
##                            shared effect rows.  One line per method.
##
## Requires helpers from fig_covsel_heatmap.R (.fix_covsel_denom, .VAR_LAB,
## .COVAR_ORD, .SHAPE_ORD, .SHAPE_LAB, %||%); it is sourced automatically.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

if (!exists(".fix_covsel_denom")) source("script/viz/fig_covsel_heatmap.R")

## ---- method sources (one fixed estimator/outer_opt cell each) --------------

.CMP_SOURCES <- list(
  "nlmixr2-SCM" = list(
    path      = "output/scm_bench_rescue_winner_aggregated/scm_covsel_by_covar.csv",
    estimator = "focei",
    outer_opt = "bobyqa"
  ),
  "PsN-SCM" = list(
    path      = "output/psn_scm_combined_aggregated/scm_covsel_by_covar.csv",
    estimator = "nonmem_scm",
    outer_opt = "focei"
  ),
  "nlmixr2-VAE" = list(
    path      = "output/vae_covsel_full0729_est702_aggregated/vae_covsel_by_covar.csv",
    estimator = NULL,
    outer_opt = NULL
  )
)

.PAL_METHOD_C <- c(
  "nlmixr2-SCM"  = "#0F8B8D",   # teal
  "PsN-SCM"      = "#2B2B2B",   # near-black
  "nlmixr2-VAE"  = "#A64AC9"    # purple
)

## ---- loader: one method -> shape-resolved per-cell error rate --------------
## structure/estimator/outer_opt fixed to ONE cell; returns per (scenario,
## effect) error rate (FN on true rows, FP on distractor rows) plus ordering
## keys and the folded effect label "PARAM~COVAR.shape".

.cmp_load_one <- function(src, structure = "linCmt", sample_N = 80) {
  d <- readr::read_csv(src$path, show_col_types = FALSE)

  if (!is.null(src$estimator) && "estimator" %in% names(d))
    d <- dplyr::filter(d, .data$estimator == src$estimator)
  if (!is.null(src$outer_opt) && "outer_opt" %in% names(d))
    d <- dplyr::filter(d, .data$outer_opt == src$outer_opt)

  ## structure label differs by platform: nlmixr2 uses linCmt; PsN uses advan4.
  struct_keep <- if (identical(structure, "linCmt"))
                   c("linCmt", "advan4") else structure
  d <- dplyr::filter(d, .data$structure %in% struct_keep,
                        .data$sample_N == !!sample_N)
  if (nrow(d) == 0L) stop("no rows for ", src$path,
                          " (structure/sample_N filter)")

  d <- .fix_covsel_denom(d)

  d |>
    dplyr::mutate(
      err_rate = ifelse(is_true, n_FN / n_datasets, n_FP / n_datasets) * 100,
      param    = dplyr::coalesce(.VAR_LAB[tolower(var)], toupper(var)),
      covar_u  = toupper(covar),
      shape_l  = tolower(shape),
      ## canonicalise the continuous-distractor shape: PsN tests "exp", nlmixr2
      ## tests "lin"; both are the SAME "wrong continuous shape" distractor, so
      ## fold them onto one row ("lin") so the 16 var-covar-shape combos align
      ## across methods (otherwise CL~BW.exp vs CL~BW.lin never join).
      shape_c  = dplyr::if_else(shape_l == "exp", "lin", shape_l),
      effect   = paste0(param, "~", covar_u, ".",
                        .SHAPE_LAB[shape_c] %||% shape_c),
      v_ord    = match(tolower(var), names(.VAR_LAB)),
      c_ord    = match(covar_u, toupper(.COVAR_ORD)),
      s_ord    = match(shape_c, .SHAPE_ORD)
    ) |>
    dplyr::filter(!is.na(v_ord), !is.na(c_ord), !is.na(s_ord)) |>
    dplyr::transmute(scenario = as.integer(scenario), effect,
                     v_ord, c_ord, s_ord, is_true, err_rate)
}

## effect-row ordering (CL block then VC; BW/BMI adjacent; true shape above its
## distractor) built from the UNION of all methods so rows are shared.
.cmp_effect_levels <- function(df) {
  df |>
    dplyr::distinct(effect, v_ord, c_ord, s_ord) |>
    dplyr::arrange(v_ord, c_ord, s_ord) |>
    dplyr::pull(effect)
}

## ---- Figure 1 : signed difference heatmaps (method - reference) ------------

fig_covsel_diff <- function(sources   = .CMP_SOURCES,
                            reference = "PsN-SCM",
                            methods   = c("nlmixr2-SCM", "nlmixr2-VAE"),
                            structure = "linCmt",
                            sample_N  = 80,
                            labels    = TRUE,
                            save      = FALSE,
                            out_dir   = "output/figures/covsel_compare",
                            base_size = 14) {
  ref <- .cmp_load_one(sources[[reference]], structure, sample_N) |>
    dplyr::rename(err_ref = err_rate) |>
    dplyr::select(scenario, effect, v_ord, c_ord, s_ord, err_ref)

  diffs <- lapply(methods, function(m) {
    .cmp_load_one(sources[[m]], structure, sample_N) |>
      dplyr::rename(err_m = err_rate, is_true_m = is_true) |>
      ## full-join on the canonical (exp->lin folded) effect so ALL 16
      ## var-covar-shape rows survive; a distractor absent in one method (never
      ## falsely selected) is a genuine 0% error, so coalesce NA -> 0 rather
      ## than dropping the row.
      dplyr::full_join(ref, by = c("scenario", "effect",
                                   "v_ord", "c_ord", "s_ord")) |>
      dplyr::mutate(
        err_m   = dplyr::coalesce(err_m, 0),
        err_ref = dplyr::coalesce(err_ref, 0),
        method  = m, delta = err_m - err_ref)
  })
  dat <- dplyr::bind_rows(diffs)

  ## complete the shared (scenario x effect) grid within each method panel.
  ## A cell absent from BOTH method and reference is a distractor never selected
  ## by either -> 0% error on both sides -> delta = 0 (white "matches"), NOT
  ## missing data.  Filling 0 removes the scattered grey (na.value) boxes.
  dat <- dat |>
    dplyr::mutate(method = factor(method, levels = methods)) |>
    tidyr::complete(method, scenario,
      tidyr::nesting(effect, v_ord, c_ord, s_ord),
      fill = list(delta = 0))

  eff_lvls <- .cmp_effect_levels(dat)
  dat <- dplyr_eff(dat, eff_lvls)

  cap <- sprintf(
    "Difference in error rate vs %s  (%s, N = %d).  Red = method worse, blue = better, white = matches %s.",
    reference, structure, sample_N, reference)

  p <- ggplot(dat, aes(factor(scenario), effect, fill = delta)) +
    geom_tile(colour = "grey85", linewidth = 0.3) +
    facet_wrap(~ method, nrow = 1) +
    scale_fill_gradient2(
      low = "#08519C", mid = "white", high = "#B30000",
      midpoint = 0, na.value = "grey92",
      limits = c(-100, 100), breaks = seq(-100, 100, 50),
      name = "delta error rate (pp)") +
    labs(x = "Simulation scenario",
         y = "Covariate effect (param ~ covariate . shape)",
         caption = cap) +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid       = element_blank(),
      legend.position  = "top",
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text       = element_text(size = base_size),
      axis.text.x      = element_text(size = base_size - 3),
      axis.text.y      = element_text(size = base_size - 4),
      plot.caption     = element_text(size = base_size - 5, hjust = 0)
    ) +
    guides(fill = guide_colourbar(barwidth = 14))

  if (isTRUE(labels)) {
    lab <- dplyr::filter(dat, !is.na(delta))
    p <- p + geom_text(data = lab, aes(label = round(delta)),
                       size = 3.2, fontface = "bold",
                       colour = ifelse(abs(lab$delta) > 55, "white", "grey20"))
  }

  if (isTRUE(save)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    stub <- file.path(out_dir,
      sprintf("fig_covsel_diff_%s_N%d", structure, sample_N))
    ggsave(paste0(stub, ".png"), p, width = 12, height = 6.5, dpi = 200)
    ggsave(paste0(stub, ".pdf"), p, width = 12, height = 6.5)
    message("saved: ", stub, ".{png,pdf}")
  }
  p
}

## small helper: apply reversed effect factor so CL~BW sits at the top
dplyr_eff <- function(df, eff_lvls) {
  dplyr::mutate(df, effect = factor(effect, levels = rev(eff_lvls)))
}

## ---- Figure 2 : per-scenario agreement vs reference ------------------------
## agreement = 100 - mean |error-rate difference| over shared effect rows.

fig_covsel_agreement <- function(sources   = .CMP_SOURCES,
                                 reference = "PsN-SCM",
                                 methods   = c("nlmixr2-SCM", "nlmixr2-VAE"),
                                 structure = "linCmt",
                                 sample_N  = 80,
                                 save      = FALSE,
                                 out_dir   = "output/figures/covsel_compare",
                                 base_size = 15) {
  ref <- .cmp_load_one(sources[[reference]], structure, sample_N) |>
    dplyr::select(scenario, effect, err_ref = err_rate)

  ag <- lapply(methods, function(m) {
    .cmp_load_one(sources[[m]], structure, sample_N) |>
      dplyr::select(scenario, effect, err_m = err_rate) |>
      dplyr::full_join(ref, by = c("scenario", "effect")) |>
      dplyr::mutate(err_m   = dplyr::coalesce(err_m, 0),
                    err_ref = dplyr::coalesce(err_ref, 0)) |>
      dplyr::group_by(scenario) |>
      dplyr::summarise(agreement = 100 - mean(abs(err_m - err_ref)),
                       .groups = "drop") |>
      dplyr::mutate(method = m)
  }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(method = factor(method, levels = methods))

  p <- ggplot(ag, aes(scenario, agreement, colour = method, group = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_colour_manual(values = .PAL_METHOD_C[methods], name = "Method") +
    scale_x_continuous(breaks = 1:16) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    coord_cartesian(ylim = c(min(60, floor(min(ag$agreement) / 5) * 5), 100)) +
    labs(x = "Simulation scenario",
         y = sprintf("Pattern agreement vs %s", reference),
         caption = sprintf(
           "Agreement = 100 - mean |delta error rate| over shared effect rows  (%s, N = %d).",
           structure, sample_N)) +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "top",
      plot.caption     = element_text(size = base_size - 5, hjust = 0)
    )

  if (isTRUE(save)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    stub <- file.path(out_dir,
      sprintf("fig_covsel_agreement_%s_N%d", structure, sample_N))
    ggsave(paste0(stub, ".png"), p, width = 9, height = 4.5, dpi = 200)
    ggsave(paste0(stub, ".pdf"), p, width = 9, height = 4.5)
    message("saved: ", stub, ".{png,pdf}")
  }
  p
}
