## fig_error_compare.R
## Cross-method comparison of ESTIMATION ACCURACY (relative-error metrics),
## shown for BOTH conditioning regimes (Unconditioned vs True selection).
##
##   fig_error_threeway()      Cleveland dot plot: y = parameter (class-grouped),
##                             x = error metric (%), colour = METHOD, faceted by
##                             conditioning (columns) so both regimes show side
##                             by side.  Direct 3-method accuracy comparison.
##
##   fig_error_threeway_diff() signed DIFFERENCE vs reference (default PsN-SCM):
##                             Delta metric = method - PsN per parameter, dot +
##                             segment from 0, vertical 0 reference.  Faceted by
##                             conditioning (columns).  Left of 0 = method MORE
##                             accurate (smaller error) than PsN.
##
## Both fix ONE scenario (default 16 = all-effects-active), ONE structure
## (linCmt for nlmixr2 / advan4 for PsN), ONE sample size (default N = 80), and
## ONE metric (default MARE = robust median |rel err|).  Reads the schema-2.1
## estimation roll-ups (scm_/vae_estim_all.csv + _cond.csv).

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

## ---- method sources (all + cond CSV per method, one fixed cell each) --------

.ERR_SOURCES <- list(
  "nlmixr2-SCM" = list(
    dir       = "output/scm_bench_rescue_winner_aggregated",
    csv_all   = "scm_estim_all.csv",
    csv_cond  = "scm_estim_cond.csv",
    estimator = "focei",
    outer_opt = "bobyqa"
  ),
  "PsN-SCM" = list(
    dir       = "output/psn_scm_combined_aggregated",
    csv_all   = "scm_estim_all.csv",
    csv_cond  = "scm_estim_cond.csv",
    estimator = "nonmem_scm",
    outer_opt = "focei"
  ),
  "nlmixr2-VAE" = list(
    dir       = "output/vae_covsel_full0729_est702_aggregated",
    csv_all   = "vae_estim_all.csv",
    csv_cond  = "vae_estim_cond.csv",
    estimator = NULL,
    outer_opt = NULL
  )
)

.PAL_METHOD_E <- c(
  "nlmixr2-SCM"  = "#0F8B8D",   # teal
  "PsN-SCM"      = "#2B2B2B",   # near-black
  "nlmixr2-VAE"  = "#A64AC9"    # purple
)

## pretty parameter labels + canonical top-down order (mirrors fig_error_metrics)
.PARAM_LAB_E <- c(
  CLBW = "CL~BW", CLcrCL = "CL~CrCL", VcBW = "Vc~BW", VcSEX = "Vc~SEX",
  TVCL = "TVCL", TVVc = "TVVc",
  TVQ = "TVQ", TVVp = "TVVp", TVKA = "TVKA",
  var_CL = "var(CL)", var_Vc = "var(Vc)", cov_VcCL = "cov(Vc,CL)",
  ResErr = "ResErr"
)
.PARAM_ORD_E <- c(
  "CLBW", "CLcrCL", "VcBW", "VcSEX",          # covariate_beta (headline)
  "TVCL", "TVVc",                             # structural_intercept
  "TVQ", "TVVp", "var_CL", "var_Vc", "cov_VcCL", "ResErr", "TVKA")
.CLASS_LAB_E <- c(
  covariate_beta       = "Covariate effects",
  structural_intercept = "Structural intercepts",
  structural_other     = "Structural / variance")

.METRIC_COL_E <- c(RMRSE = "RMRSE_pct", MARE = "MARE_pct")

## ---- loader: one method, both regimes -> long per-parameter error ----------

.err_load_one <- function(src, scenario = 16, structure = "linCmt",
                          sample_N = 80, metric = "MARE") {
  mcol <- unname(.METRIC_COL_E[metric])
  struct_keep <- if (identical(structure, "linCmt"))
                   c("linCmt", "advan4") else structure

  read_one <- function(file, cond_label) {
    p <- file.path(src$dir, file)
    if (!file.exists(p)) stop("estim CSV not found: ", p)
    d <- suppressMessages(readr::read_csv(p, show_col_types = FALSE))
    if (!is.null(src$estimator) && "estimator" %in% names(d))
      d <- dplyr::filter(d, .data$estimator == src$estimator)
    if (!is.null(src$outer_opt) && "outer_opt" %in% names(d))
      d <- dplyr::filter(d, .data$outer_opt == src$outer_opt)
    d |>
      dplyr::filter(.data$scenario == !!scenario,
                    .data$structure %in% struct_keep,
                    .data$sample_N == !!sample_N) |>
      dplyr::transmute(parameter, param_class,
                       value = .data[[mcol]],
                       conditioning = cond_label)
  }

  dplyr::bind_rows(
    read_one(src$csv_all,  "Unconditioned"),
    read_one(src$csv_cond, "True selection")
  )
}

## stack all methods, drop fixed TVKA, attach ordering/labels
.err_load_all <- function(sources, methods, scenario, structure,
                          sample_N, metric, drop_fixed = TRUE) {
  dat <- lapply(methods, function(m) {
    .err_load_one(sources[[m]], scenario, structure, sample_N, metric) |>
      dplyr::mutate(method = m)
  }) |> dplyr::bind_rows()

  if (drop_fixed) dat <- dplyr::filter(dat, parameter != "TVKA")

  present <- .PARAM_ORD_E[.PARAM_ORD_E %in% unique(dat$parameter)]
  lab_present <- dplyr::coalesce(.PARAM_LAB_E[present], present)

  dat |>
    dplyr::mutate(
      method    = factor(method, levels = methods),
      param_lab = dplyr::coalesce(.PARAM_LAB_E[parameter], parameter),
      param_lab = factor(param_lab, levels = rev(lab_present)),
      class_lab = factor(dplyr::coalesce(.CLASS_LAB_E[param_class], param_class),
                         levels = unname(.CLASS_LAB_E)),
      conditioning = factor(conditioning,
                            levels = c("Unconditioned", "True selection"))
    )
}

.err_theme <- function(base_size = 15) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "top",
      strip.text         = element_text(size = base_size),
      strip.background   = element_rect(fill = "grey92", colour = NA),
      strip.text.y       = element_blank(),
      strip.background.y = element_blank(),
      axis.text.y        = element_text(size = base_size - 2),
      axis.text.x        = element_text(size = base_size - 2),
      # widen the gap between the Unconditioned | True selection facet columns
      # so their 0% / 60% axis ticks don't collide.
      panel.spacing.x    = grid::unit(2, "lines"),
      plot.caption       = element_text(size = base_size - 3.5, hjust = 0,
                                        colour = "grey35")
    )
}

## ---- Figure 1 : three-method Cleveland dots, both regimes ------------------

fig_error_threeway <- function(sources    = .ERR_SOURCES,
                               methods    = c("nlmixr2-SCM", "PsN-SCM",
                                              "nlmixr2-VAE"),
                               scenario   = 16,
                               structure  = "linCmt",
                               sample_N   = 80,
                               metric     = "MARE",
                               x_cap      = 100,
                               save       = FALSE,
                               out_dir    = "output/figures/error_compare",
                               base_size  = 15) {
  metric <- match.arg(metric, names(.METRIC_COL_E))
  dat <- .err_load_all(sources, methods, scenario, structure, sample_N, metric)

  x_up  <- max(20, x_cap)
  dat <- dplyr::mutate(dat,
    value_plot = pmin(value, x_up),
    off_scale  = !is.na(value) & value > x_up)

  dodge <- position_dodge(width = 0.6)
  p <- ggplot(dat, aes(value_plot, param_lab,
                       colour = method, group = method)) +
    geom_linerange(aes(xmin = 0, xmax = value_plot),
                   position = dodge, linewidth = 0.8, alpha = 0.45) +
    geom_point(position = dodge, size = 2.8) +
    facet_grid(class_lab ~ conditioning, scales = "free_y", space = "free_y") +
    scale_colour_manual(values = .PAL_METHOD_E[methods], name = "Method") +
    scale_x_continuous(labels = function(x) paste0(x, "%"),
                       limits = c(0, x_up), oob = scales::oob_squish,
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = sprintf("%s (%%)", metric), y = "Parameter",
         caption = sprintf(
           "Estimation accuracy, scenario %d  (%s, N = %d).  Lower = more accurate.  Values > %d%% squished to edge (italic label).",
           scenario, structure, sample_N, x_up)) +
    .err_theme(base_size)

  if (any(dat$off_scale)) {
    p <- p + geom_text(data = dplyr::filter(dat, off_scale),
      aes(label = paste0(round(value), "%")),
      position = dodge, hjust = 1.1, size = 3.0, fontface = "italic",
      show.legend = FALSE)
  }

  if (isTRUE(save)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    stub <- file.path(out_dir,
      sprintf("fig_error_threeway_scn%d_%s_N%d_%s",
              scenario, structure, sample_N, metric))
    ggsave(paste0(stub, ".png"), p, width = 7, height = 7, dpi = 200)
    ggsave(paste0(stub, ".pdf"), p, width = 7, height = 7)
    message("saved: ", stub, ".{png,pdf}")
  }
  p
}

## ---- Figure 2 : difference vs reference, both regimes ----------------------

fig_error_threeway_diff <- function(sources    = .ERR_SOURCES,
                                    reference  = "PsN-SCM",
                                    methods    = c("nlmixr2-SCM", "nlmixr2-VAE"),
                                    scenario   = 16,
                                    structure  = "linCmt",
                                    sample_N   = 80,
                                    metric     = "MARE",
                                    save       = FALSE,
                                    out_dir    = "output/figures/error_compare",
                                    base_size  = 15) {
  metric <- match.arg(metric, names(.METRIC_COL_E))
  all_m  <- unique(c(reference, methods))
  dat <- .err_load_all(sources, all_m, scenario, structure, sample_N, metric)

  ref <- dat |>
    dplyr::filter(method == reference) |>
    dplyr::select(parameter, param_lab, class_lab, conditioning,
                  value_ref = value)

  diff <- dat |>
    dplyr::filter(method %in% methods) |>
    dplyr::inner_join(ref, by = c("parameter", "param_lab", "class_lab",
                                  "conditioning")) |>
    dplyr::mutate(delta = value - value_ref,
                  method = factor(method, levels = methods))

  dodge <- position_dodge(width = 0.6)
  lim <- max(abs(diff$delta), na.rm = TRUE)
  p <- ggplot(diff, aes(delta, param_lab, colour = method, group = method)) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.5) +
    geom_linerange(aes(xmin = 0, xmax = delta),
                   position = dodge, linewidth = 0.8, alpha = 0.45) +
    geom_point(position = dodge, size = 2.8) +
    facet_grid(class_lab ~ conditioning, scales = "free_y", space = "free_y") +
    scale_colour_manual(values = .PAL_METHOD_E[methods], name = "Method") +
    scale_x_continuous(labels = function(x) paste0(x, "pp"),
                       limits = c(-lim, lim),
                       expand = expansion(mult = 0.05)) +
    labs(x = sprintf("Delta %s vs %s (pp)", metric, reference),
         y = "Parameter",
         caption = sprintf(
           "Difference in %s vs %s, scenario %d  (%s, N = %d).  Left of 0 = more accurate than %s; right = less accurate.",
           metric, reference, scenario, structure, sample_N, reference)) +
    .err_theme(base_size)

  if (isTRUE(save)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    stub <- file.path(out_dir,
      sprintf("fig_error_threeway_diff_scn%d_%s_N%d_%s",
              scenario, structure, sample_N, metric))
    ggsave(paste0(stub, ".png"), p, width = 11, height = 7, dpi = 200)
    ggsave(paste0(stub, ".pdf"), p, width = 11, height = 7)
    message("saved: ", stub, ".{png,pdf}")
  }
  p
}
