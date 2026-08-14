# ==============================================================================
# Covariate correlation figure (Fig 2 style) across three virtual-population
# cohorts (N = 40, 80, 300), restricted to the first 100 replicate datasets.
# ------------------------------------------------------------------------------
#   Inputs (relative to workspace root):
#     simulated_virtual_dataset_eta_filtered_N40/sim_obs_all_scenarios.rds
#     simulated_virtual_dataset_eta_filtered_N80/sim_obs_all_scenarios.rds
#     simulated_virtual_dataset_eta_filtered/sim_obs_all_scenarios.rds     (N=300)
#     C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/NHANES_dataset/
#         {DEMO_L, BMX_L, BIOPRO_L}.XPT
#
#   Outputs:
#     simulated_virtual_dataset/ref_cor.rds
#     output/figures/correlationMatrix/covariate_cor_summary.csv
#     output/figures/correlationMatrix/fig2_covariate_correlation_{N40,N80,N300}.png
#
#   Design notes:
#     * Only datasets 1-100 are used (DATASET <= DATASET_MAX); set DATASET_MAX
#       below to change the range.
#     * Covariates in each sim dataset are scenario-invariant (iCov is shared
#       across the 16 scenarios in simulate_scenario_v2()), so we dedupe on
#       (DATASET, SUBJECT) within SCENARIO == 1.
#     * Variable order matches the reference Fig 2: BMI, BW, CrCL, RACE, SEX.
#     * Upper triangle = filled circles (size ~ |r|) via corrplot.
#     * Lower triangle = two-line label "r \n +/- SD" from the 250 per-
#       dataset correlations, overlaid manually with base text().
#     * Diagonal = variable name in steelblue via tl.pos = "d".
#     * NHANES XPT files are required; the script stops with an explicit
#       error pointing at the expected folder if any are missing.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(haven)
  library(corrplot)
})

# ---- Paths -----------------------------------------------------------------
NHANES_DIR <- "C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/NHANES_dataset"

cohort_paths <- c(
  N40  = "Inputdataset/sim_obs_N40/sim_obs_scenario_01.rds",
  N80  = "Inputdataset/sim_obs_N80/sim_obs_scenario_01.rds",
  N300 = "Inputdataset/sim_obs_N300/sim_obs_scenario_01.rds"
)

out_fig_dir <- "output/figures/correlationMatrix"
out_ref_dir <- "simulated_virtual_dataset"
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_ref_dir, showWarnings = FALSE, recursive = TRUE)

VAR_ORDER <- c("BMI", "BW", "CrCL", "RACE", "SEX")

# Focus the analysis on the first 100 replicate datasets (1-100) only.
DATASET_MAX <- 100L

# ==============================================================================
# 1. Rebuild NHANES reference correlation
# ------------------------------------------------------------------------------
#   Mirrors the block in scripts/Refitmodel_robustness.R (lines 45-66):
#     - Non-Hispanic White (3) + Non-Hispanic Asian (6), age > 17
#     - Cockcroft-Gault CrCL, physiologic filter (0 < CrCL < 250)
# ==============================================================================
xpt_paths <- c(
  demo   = file.path(NHANES_DIR, "DEMO_L.XPT"),
  bmx    = file.path(NHANES_DIR, "BMX_L.XPT"),
  biopro = file.path(NHANES_DIR, "BIOPRO_L.XPT")
)

missing_xpt <- xpt_paths[!file.exists(xpt_paths)]
if (length(missing_xpt) > 0L) {
  stop(
    "Missing NHANES XPT file(s):\n  ",
    paste(missing_xpt, collapse = "\n  "),
    "\nExpected under: ", NHANES_DIR
  )
}

demo <- haven::read_xpt(xpt_paths[["demo"]]) |>
  dplyr::select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH3)
bmx <- haven::read_xpt(xpt_paths[["bmx"]]) |>
  dplyr::select(SEQN, BMXWT, BMXHT, BMXBMI)
biopro <- haven::read_xpt(xpt_paths[["biopro"]]) |>
  dplyr::select(SEQN, LBXSCR)

pop <- demo |>
  dplyr::inner_join(bmx,    by = "SEQN") |>
  dplyr::inner_join(biopro, by = "SEQN") |>
  dplyr::filter(RIDAGEYR > 17, RIDRETH3 %in% c(3, 6)) |>
  dplyr::mutate(
    SEX  = ifelse(RIAGENDR == 2, 0, 1),
    RACE = ifelse(RIDRETH3  == 6, 0, 1),
    BW   = BMXWT,
    BMI  = BMXBMI,
    CrCL = ((140 - RIDAGEYR) * BW / (72 * LBXSCR)) * ifelse(SEX == 0, 0.85, 1)
  ) |>
  dplyr::select(dplyr::all_of(VAR_ORDER)) |>
  tidyr::drop_na() |>
  dplyr::filter(CrCL > 0, CrCL < 250)

ref_cor <- stats::cor(pop, method = "pearson")
saveRDS(ref_cor, file.path(out_ref_dir, "ref_cor.rds"))

# ==============================================================================
# 2. Load sim datasets + dedupe to subject-level covariates
# ------------------------------------------------------------------------------
#   Filter to SCENARIO == 1 (covariates are shared across all 16 scenarios)
#   then distinct(DATASET, SUBJECT, <covariates>). Result: 250 * N_subj rows.
# ==============================================================================
load_covariates <- function(path, var_order = VAR_ORDER,
                            dataset_max = DATASET_MAX) {
  if (!file.exists(path)) {
    stop("Missing sim dataset: ", path)
  }
  d <- readRDS(path)
  # Covariates are scenario-invariant; these files hold scenario 1 only. Keep
  # the SCENARIO filter defensively in case a combined file is supplied.
  if ("SCENARIO" %in% names(d)) d <- dplyr::filter(d, SCENARIO == 1L)
  d |>
    dplyr::filter(DATASET <= dataset_max) |>
    dplyr::distinct(DATASET, SUBJECT,
                    dplyr::across(dplyr::all_of(var_order))) |>
    dplyr::select(DATASET, SUBJECT, dplyr::all_of(var_order))
}

cohort_cov <- purrr::map(cohort_paths, load_covariates)

# Sanity: expect DATASET_MAX datasets per cohort
cohort_n_ds <- vapply(cohort_cov, function(d) dplyr::n_distinct(d$DATASET), integer(1))
if (any(cohort_n_ds != DATASET_MAX)) {
  warning("Unexpected dataset counts: ",
          paste(sprintf("%s=%d", names(cohort_n_ds), cohort_n_ds), collapse = ", "))
}

# ==============================================================================
# 3. Pooled and per-DATASET Pearson correlation matrices per cohort
# ==============================================================================
pooled_cor <- function(cov_df, var_order = VAR_ORDER) {
  stats::cor(as.matrix(cov_df[, var_order, drop = FALSE]),
             method = "pearson")
}

per_ds_cor_long <- function(cov_df, var_order = VAR_ORDER) {
  cov_df |>
    dplyr::group_by(DATASET) |>
    dplyr::group_modify(function(sub, key) {
      m <- stats::cor(as.matrix(sub[, var_order, drop = FALSE]),
                      method = "pearson")
      tibble::as_tibble(m, rownames = "var1") |>
        tidyr::pivot_longer(-var1, names_to = "var2", values_to = "r")
    }) |>
    dplyr::ungroup()
}

summarise_cor_long <- function(per_ds_long) {
  per_ds_long |>
    dplyr::group_by(var1, var2) |>
    dplyr::summarise(
      mean_r   = mean(r,   na.rm = TRUE),
      sd_r     = stats::sd(r, na.rm = TRUE),
      median_r = stats::median(r, na.rm = TRUE),
      q05_r    = stats::quantile(r, 0.05, na.rm = TRUE, names = FALSE),
      q95_r    = stats::quantile(r, 0.95, na.rm = TRUE, names = FALSE),
      .groups  = "drop"
    )
}

pooled_mats <- purrr::map(cohort_cov, pooled_cor)
per_ds_long <- purrr::map(cohort_cov, per_ds_cor_long)
per_ds_sum  <- purrr::map(per_ds_long, summarise_cor_long)

# Build SD matrix per cohort (same row/col order as pooled)
build_sd_matrix <- function(sum_df, var_order = VAR_ORDER) {
  m <- matrix(NA_real_, nrow = length(var_order), ncol = length(var_order),
              dimnames = list(var_order, var_order))
  for (i in seq_along(var_order)) {
    for (j in seq_along(var_order)) {
      hit <- sum_df |>
        dplyr::filter(var1 == var_order[i], var2 == var_order[j])
      if (nrow(hit) == 1L) m[i, j] <- hit$sd_r
    }
  }
  m
}
sd_mats <- purrr::map(per_ds_sum, build_sd_matrix)

# ==============================================================================
# 4. Tidy summary table incl. ref_cor preservation flag
# ==============================================================================
ref_long <- tibble::as_tibble(ref_cor, rownames = "var1") |>
  tidyr::pivot_longer(-var1, names_to = "var2", values_to = "ref_r")

cohort_summary <- purrr::imap_dfr(per_ds_sum, function(sum_df, cohort) {
  pooled <- pooled_mats[[cohort]]
  pooled_long <- tibble::as_tibble(pooled, rownames = "var1") |>
    tidyr::pivot_longer(-var1, names_to = "var2", values_to = "pooled_r")
  sum_df |>
    dplyr::left_join(pooled_long, by = c("var1", "var2")) |>
    dplyr::left_join(ref_long,    by = c("var1", "var2")) |>
    dplyr::mutate(
      cohort        = cohort,
      abs_delta_ref = abs(pooled_r - ref_r),
      within_0.05   = abs_delta_ref <= 0.05
    ) |>
    dplyr::select(cohort, var1, var2, pooled_r, sd_r, mean_r,
                  median_r, q05_r, q95_r, ref_r, abs_delta_ref, within_0.05)
})

readr::write_csv(cohort_summary,
                 file.path(out_fig_dir, "covariate_cor_summary.csv"))

# Warn on any off-diagonal preservation breach
breach <- cohort_summary |>
  dplyr::filter(var1 != var2, !within_0.05)
if (nrow(breach) > 0L) {
  warning("Bootstrap preservation breached (|pooled - ref| > 0.05) in ",
          nrow(breach), " off-diagonal cells:\n",
          paste(sprintf("  %s: %s ~ %s  delta=%.3f",
                        breach$cohort, breach$var1, breach$var2,
                        breach$abs_delta_ref),
                collapse = "\n"))
}

# ==============================================================================
# 5. Fig 2-style renderer: upper circles + lower two-line "r \n +/- SD"
# ------------------------------------------------------------------------------
#   corrplot places matrix cell (i, j) at plot coordinate (j, n + 1 - i)
#   for an n x n matrix drawn top-left origin. We overlay the strict lower
#   triangle manually with text() to get the two-line "r \n +/- SD" labels
#   (corrplot's built-in addCoef.col only supports single-line values).
# ==============================================================================
render_fig2_panel <- function(pooled, sd_mat = NULL, title, out_png,
                              var_order = VAR_ORDER,
                              width_px = 2400, height_px = 2400, res_dpi = 300) {
  stopifnot(identical(rownames(pooled), var_order),
            identical(colnames(pooled), var_order))
  if (!is.null(sd_mat)) {
    stopifnot(identical(dim(pooled), dim(sd_mat)))
  }
  n <- length(var_order)

  grDevices::png(out_png, width = width_px, height = height_px, res = res_dpi)
  on.exit(grDevices::dev.off(), add = TRUE)

  # Mask lower triangle + diagonal to NA and draw the FULL 5x5 layout so the
  # corrplot plot coordinate system stays a straightforward (col, n+1-row)
  # grid. NA cells render blank, leaving the diagonal and lower triangle
  # free for our manual overlay. Circles in the upper triangle are drawn by
  # corrplot with radius proportional to |r|.
  upper_only <- pooled
  upper_only[lower.tri(upper_only, diag = TRUE)] <- NA_real_

  corrplot::corrplot(
    upper_only,
    method       = "circle",
    type         = "full",
    diag         = TRUE,
    col          = "black",
    bg           = "white",
    tl.pos       = "n",                   # suppress corrplot's own labels
    cl.pos       = "n",                   # suppress colour legend
    addgrid.col  = "grey60",
    na.label     = " ",                   # blank NAs (lower tri + diag)
    na.label.col = "white",
    mar          = c(0, 0, 2, 0),
    title        = title
  )

  # Diagonal: variable name in steelblue
  for (k in seq_len(n)) {
    graphics::text(
      x      = k,
      y      = n + 1 - k,
      labels = var_order[k],
      col    = "steelblue",
      cex    = 1.4,
      font   = 2
    )
  }

  # Strict lower triangle overlay:
  #   sd_mat supplied -> two-line  "r \n +/- SD"  (cohort panels: 250 datasets)
  #   sd_mat = NULL   -> single-line "r"          (reference panel: 1 dataset)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i > j) {
        lbl <- if (is.null(sd_mat)) {
          sprintf("%.2f", pooled[i, j])
        } else {
          sprintf("%.2f\n\u00b1%.2f", pooled[i, j], sd_mat[i, j])
        }
        graphics::text(
          x      = j,
          y      = n + 1 - i,
          labels = lbl,
          cex    = if (is.null(sd_mat)) 1.2 else 1.0,
          font   = 2
        )
      }
    }
  }
}

cohort_titles <- c(
  N40  = sprintf("N = 40 per dataset (%d datasets)", DATASET_MAX),
  N80  = sprintf("N = 80 per dataset (%d datasets)", DATASET_MAX),
  N300 = sprintf("N = 300 per dataset (%d datasets)", DATASET_MAX)
)

for (cohort in names(cohort_paths)) {
  render_fig2_panel(
    pooled  = pooled_mats[[cohort]],
    sd_mat  = sd_mats[[cohort]],
    title   = cohort_titles[[cohort]],
    out_png = file.path(out_fig_dir,
                        sprintf("fig2_covariate_correlation_%s.png", cohort))
  )
}

# ---- NHANES reference panel ------------------------------------------------
#   Single dataset (no 250-replicate distribution) so the lower triangle
#   shows only r, matching the original Fig 2 reference layout. Placed
#   side-by-side with the three cohort panels this visually demonstrates
#   that the stratified bootstrap preserves the NHANES correlation
#   structure at every N (quantified in outputs/covariate_cor_summary.csv).
render_fig2_panel(
  pooled  = ref_cor,
  sd_mat  = NULL,
  title   = sprintf("NHANES 2021-2023 reference (N = %d)", nrow(pop)),
  out_png = file.path(out_fig_dir, "fig2_covariate_correlation_ref.png")
)

# ==============================================================================
# 6. Return object for interactive inspection
# ==============================================================================
list(
  ref_cor         = ref_cor,
  pooled_matrices = pooled_mats,
  sd_matrices     = sd_mats,
  per_dataset     = per_ds_sum,
  summary_df      = cohort_summary
)
