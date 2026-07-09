dir.create(file.path(tempdir(), "vscode-R"), showWarnings = FALSE)
library(testthat)
library(nlmixr2scm)
library(nlmixr2utils)
#devtools::load_all("C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm")
library(nlmixr2)
library(rxode2)
library(tidyverse)

library(devtools)
library(haven)
test_ds02_m02 <- readRDS("simulated_virtual_dataset_eta_filtered_N80/stage1_pilot_scn16_N80_m02/test_full_fast_ds02.rds")
test_ds02 <- readRDS("simulated_virtual_dataset_eta_filtered_N80/stage1_pilot_scn16_N80/test_full_fast_ds02.rds")
test_ds03 <- readRDS("simulated_virtual_dataset_eta_filtered_N80/stage1_pilot_scn16_N80/test_full_fast_ds03.rds")
test_ds13 <- readRDS("output_N80/test_ds13.rds")
res_refit_true_N300sc16ds02 <- readRDS("output/refit_true_N300/scn16_none/res_ds002.rds")


true_params <- readRDS("Inputdataset/true_params_long.rds")
## Pop params for RMRSE (TVKA dropped: fixed in model)
pop_params_rmrse <- c("TVCL", "TVVc", "TVQ", "TVVp",
                      "var_CL", "var_Vc", "cov_VcCL", "ResErr")
cov_params_rmrse <- c("CLBW", "CLcrCL", "VcBW", "VcSEX")

## ============================================================================
## HPCE output_N80 aggregation from saved test_ds*.rds
## ----------------------------------------------------------------------------
## Reads packaged test objects from output directories and computes:
##   - power-related metrics (via existing helper if available; fallback: selected_any rate)
##   - mean root squared relative error (MRSRE; via helper if available; fallback from parFixed)
## ============================================================================
## ============================================================================
# Part 3: Operating-characteristics pilot (4 datasets, scenario 16, N=80)-----------

##   Wraps the fast-route full-SCM pipeline used above (single dataset) into
##   a per-dataset driver, loops 4 replicates of scenario 16, then computes
##   six OC artefacts:
##     - oc_pilot_power       : Power, PowerCN, PowerMinSuc
##     - oc_pilot_relpower    : relative power at k = 1..Ntrue
##     - oc_pilot_rmrse_uncond: unconditional RMRSE %
##     - oc_pilot_rmrse_cond  : conditional RMRSE % (exact-match runs only)
##     - oc_pilot_timing      : per-dataset + aggregate wall-clock
##     - oc_pilot_per_ds      : audit log (one row per dataset)
##
##   Helpers below are intentionally scenario-agnostic so the same code drops
##   onto the 100 x 16 scale-up.  See the SCALE-UP STUB at the end.
## ============================================================================


# ---- (P3.1) OC helpers (reusable across scenarios) --------------------------

##  Canonical-shape normalizer.  runSCM emits per-LEVEL theta names for
##  categorical covariates (e.g. cov_SEX_1_vc, cov_RACE_2_vc); the parser
##  in package_scm_result() leaves the level digit in `shape`.  For
##  semantic matching against truth (where the relation is one entity)
##  we collapse any all-digit shape to "cat".  This also folds a 3-level
##  cat that selected multiple levels into ONE relation row.
.canon_shape_tbl <- function(rel) {
  if (nrow(rel) == 0L) {
    return(tibble::tibble(var   = character(),
                          covar = character(),
                          shape = character()))
  }
  rel %>%
    dplyr::mutate(shape = dplyr::if_else(
      grepl("^[0-9]+$", shape), "cat", as.character(shape)
    )) %>%
    dplyr::distinct(var, covar, shape)
}

##  Build the truth tibble of "what runSCM should have selected" for a given
##  scenario.  Reads true_long, keeps rows with non-zero true_value AND a
##  parameter name that maps to a (var, covar, shape) triple.  Defaults
##  cover the article scenarios (scn 1-16); pass `shape_map` if you add
##  custom covariate effects later.
extract_true_relations <- function(
    true_long, scenario_id,
    shape_map = list(
      CLBW   = list(var = "cl", covar = "BW",   shape = "power"),
      CLcrCL = list(var = "cl", covar = "CrCL", shape = "power"),
      VcBW   = list(var = "vc", covar = "BW",   shape = "power"),
      VcSEX  = list(var = "vc", covar = "SEX",  shape = "cat")
    )) {
  hits <- true_long %>%
    dplyr::filter(scenario == scenario_id,
                  parameter %in% names(shape_map),
                  !is.na(true_value), true_value != 0) %>%
    dplyr::pull(parameter)
  if (length(hits) == 0L) {
    return(tibble::tibble(var = character(), covar = character(),
                          shape = character()))
  }
  purrr::map_dfr(hits, function(p) {
    m <- shape_map[[p]]
    tibble::tibble(var = m$var, covar = m$covar, shape = m$shape)
  })
}

##  Compare a single run's `selected` tibble against truth.  Returns a
##  list with: n_true_hit, n_false_pos, exact_match (set equality on
##  (var, covar, shape) after canonicalization).
match_selected_to_truth <- function(selected, true_rel) {
  if (is.null(selected) || nrow(selected) == 0L) {
    return(list(n_true_hit  = 0L,
                n_false_pos = 0L,
                exact_match = nrow(true_rel) == 0L))
  }
  sel <- .canon_shape_tbl(selected)
  tru <- .canon_shape_tbl(true_rel)
  true_hit  <- dplyr::inner_join(sel, tru, by = c("var", "covar", "shape"))
  false_pos <- dplyr::anti_join (sel, tru, by = c("var", "covar", "shape"))
  miss      <- dplyr::anti_join (tru, sel, by = c("var", "covar", "shape"))
  list(
    n_true_hit  = nrow(true_hit),
    n_false_pos = nrow(false_pos),
    exact_match = (nrow(false_pos) == 0L) && (nrow(miss) == 0L)
  )
}

##  Compute Power, PowerCN (cond_num_cor < cn_cor_cut), PowerMinSuc
##  (converged == TRUE), plus a relative-power tibble (k = 1..n_true).
##  cond_num_cor is the correlation-matrix condition number (lambda_max /
##  lambda_min) from fit$conditionNumberCor.  Convention: < 1000 is "well
##  conditioned" (FDA / pharma pop-PK threshold).  We do NOT take a sqrt
##  any more -- the old cond_num_sqrt was sqrt(cond_num_cov), and Cor and
##  Cov condition numbers are on different scales (Cor is more sensitive
##  to high parameter correlations, Cov is dominated by scale differences).
compute_power_block <- function(per_ds, n_true, cn_cor_cut = 1000) {
  N      <- nrow(per_ds)
  ok_cn  <- per_ds$cond_num_cor < cn_cor_cut & !is.na(per_ds$cond_num_cor)
  ok_min <- per_ds$converged
  ok_min[is.na(ok_min)] <- FALSE

  power_main <- tibble::tibble(
    metric = c("Power", "PowerCN", "PowerMinSuc"),
    num    = c(
      sum(per_ds$exact_match, na.rm = TRUE),
      sum(per_ds$exact_match & ok_cn,  na.rm = TRUE),
      sum(per_ds$exact_match & ok_min, na.rm = TRUE)
    ),
    denom  = c(N, sum(ok_cn), sum(ok_min))
  ) %>%
    dplyr::mutate(value = num / denom)

  rel_power <- tibble::tibble(
    k                   = seq_len(n_true),
    n_at_least_k        = vapply(seq_len(n_true), function(k)
      sum(per_ds$n_true_hit >= k, na.rm = TRUE), integer(1)),
    fraction_at_least_k = vapply(seq_len(n_true), function(k)
      mean(per_ds$n_true_hit >= k, na.rm = TRUE), numeric(1))
  )

  list(power = power_main, rel_power = rel_power)
}

##  Compute unconditional + conditional RMRSE / MARE for a parameter set.
##    pop_params: vector of population-parameter names (always estimated).
##                Unconditional denom = runs with finite estimate.
##    cov_params: vector of true covariate-effect parameter names.
##                Unconditional denom = runs where that cov was SELECTED
##                (and therefore estimated, => non-NA estimate).
##    Conditional denom for BOTH groups = runs with exact_match == TRUE.
##
##  Metrics returned per parameter:
##    RMRSE_pct: 100 * sqrt(mean(rel_err^2))   -- mean accuracy + variance
##                                                penalty, sensitive to outliers
##    MARE_pct:  100 * median(|rel_err|)       -- robust central tendency,
##                                                outlier-resistant
##    n_used:    # of finite contributions
##
##  Rationale: with N = 4 (pilot) or even N = 250 (scale-up) draws, a single
##  numerically-unstable estimate (e.g. cov_VcCL near the noise floor) can
##  dominate RMRSE.  MARE gives a complementary view that's stable to that
##  single bad draw.  Report both side-by-side.
compute_rmrse_block <- function(per_ds, true_long, scenario_id,
                                pop_params, cov_params) {
  truth <- true_long %>%
    dplyr::filter(scenario == scenario_id,
                  parameter %in% c(pop_params, cov_params)) %>%
    dplyr::select(parameter, true_value)

  long <- purrr::imap_dfr(per_ds$final_est, function(est, i) {
    if (is.null(est)) return(NULL)
    dplyr::mutate(est,
                  dataset_id  = per_ds$dataset_id[i],
                  exact_match = per_ds$exact_match[i])
  }) %>%
    dplyr::inner_join(truth, by = "parameter") %>%
    dplyr::mutate(
      rel_err = (estimate - true_value) / true_value,
      sq_rel  = rel_err^2,
      abs_rel = abs(rel_err)
    )

  rmrse_one <- function(df) {
    df <- dplyr::filter(df, is.finite(sq_rel))
    if (nrow(df) == 0L) {
      tibble::tibble(RMRSE_pct = NA_real_, MARE_pct = NA_real_, n_used = 0L)
    } else {
      tibble::tibble(
        RMRSE_pct = 100 * sqrt(mean(df$sq_rel)),
        MARE_pct  = 100 * stats::median(df$abs_rel),
        n_used    = nrow(df)
      )
    }
  }

  uncond <- long %>%
    dplyr::group_by(parameter) %>%
    dplyr::group_modify(~ rmrse_one(.x)) %>%
    dplyr::ungroup()

  cond <- long %>%
    dplyr::filter(exact_match) %>%
    dplyr::group_by(parameter) %>%
    dplyr::group_modify(~ rmrse_one(.x)) %>%
    dplyr::ungroup()

  list(rmrse_uncond = uncond, rmrse_cond = cond)
}



compute_rmrse_block <- function(per_ds, true_long, scenario_id,
                                pop_params, cov_params) {
  truth <- true_long %>%
    dplyr::filter(scenario == scenario_id,
                  parameter %in% c(pop_params, cov_params)) %>%
    dplyr::select(parameter, true_value)

  empty <- tibble::tibble(parameter = character(),
                          RMRSE_pct = numeric(),
                          MARE_pct  = numeric(),
                          n_used    = integer())

  # Guard: no final_est column, or all NULL
  if (!"final_est" %in% names(per_ds) ||
      all(vapply(per_ds$final_est, is.null, logical(1)))) {
    warning("compute_rmrse_block: no final_est available; returning empty tibbles")
    return(list(rmrse_uncond = empty, rmrse_cond = empty))
  }

  long <- purrr::imap_dfr(per_ds$final_est, function(est, i) {
    if (is.null(est) || !"parameter" %in% names(est)) return(NULL)
    dplyr::mutate(est,
                  dataset_id  = per_ds$dataset_id[i],
                  exact_match = per_ds$exact_match[i])
  })

  if (nrow(long) == 0L) {
    warning("compute_rmrse_block: no rows to join; returning empty tibbles")
    return(list(rmrse_uncond = empty, rmrse_cond = empty))
  }

  long <- long %>%
    dplyr::inner_join(truth, by = "parameter") %>%
    dplyr::mutate(
      rel_err = (estimate - true_value) / true_value,
      sq_rel  = rel_err^2,
      abs_rel = abs(rel_err)
    )

  rmrse_one <- function(df) {
    df <- dplyr::filter(df, is.finite(sq_rel))
    if (nrow(df) == 0L) {
      tibble::tibble(RMRSE_pct = NA_real_, MARE_pct = NA_real_, n_used = 0L)
    } else {
      tibble::tibble(
        RMRSE_pct = 100 * sqrt(mean(df$sq_rel)),
        MARE_pct  = 100 * stats::median(df$abs_rel),
        n_used    = nrow(df)
      )
    }
  }

  uncond <- long %>% dplyr::group_by(parameter) %>%
    dplyr::group_modify(~ rmrse_one(.x)) %>% dplyr::ungroup()
  cond   <- long %>% dplyr::filter(exact_match) %>% dplyr::group_by(parameter) %>%
    dplyr::group_modify(~ rmrse_one(.x)) %>% dplyr::ungroup()

  list(rmrse_uncond = uncond, rmrse_cond = cond)
}

## ---------------------------------------------------------------------------
## extract_individual_rse()
##
## Returns a long tibble with ONE ROW PER (dataset_id, parameter), containing
## the raw estimate, the true value, the relative error, its square, and its
## absolute value.  This is the "atomic" table that compute_rmrse_block()
## aggregates -- exposing it lets you find the datasets that dominate RMRSE.
##
## Args:
##   per_ds        : oc_full_per_ds (must have $dataset_id, $final_est,
##                   $exact_match, and ideally $cond_num_cor, $converged)
##   true_long     : true_params (long-format truth table)
##   scenario_id   : scenario integer (e.g. 16)
##   params        : character vector of parameter names to extract.
##                   Default = c("var_CL", "cov_VcCL").  Pass NULL to get
##                   every parameter present in `truth` for that scenario.
##
## Returns: tibble sorted by |rel_err| descending, so the worst offenders
##          are on top.
## ---------------------------------------------------------------------------
extract_individual_rse <- function(per_ds, true_long, scenario_id,
                                   params = c("var_CL", "cov_VcCL")) {

  truth <- true_long %>%
    dplyr::filter(scenario == scenario_id) %>%
    dplyr::select(parameter, true_value)

  if (!is.null(params)) {
    truth <- dplyr::filter(truth, parameter %in% params)
  }
  if (nrow(truth) == 0L) {
    warning("extract_individual_rse: no matching parameters in truth")
    return(tibble::tibble())
  }

  long <- purrr::imap_dfr(per_ds$final_est, function(est, i) {
    if (is.null(est) || !"parameter" %in% names(est)) return(NULL)
    dplyr::mutate(
      est,
      dataset_id   = per_ds$dataset_id[i],
      exact_match  = per_ds$exact_match[i],
      converged    = if ("converged"    %in% names(per_ds)) per_ds$converged[i]    else NA,
      cond_num_cor = if ("cond_num_cor" %in% names(per_ds)) per_ds$cond_num_cor[i] else NA_real_
    )
  })

  if (nrow(long) == 0L) {
    warning("extract_individual_rse: no final_est rows")
    return(tibble::tibble())
  }

  long %>%
    dplyr::inner_join(truth, by = "parameter") %>%
    dplyr::mutate(
      rel_err    = (estimate - true_value) / true_value,
      sq_rel     = rel_err^2,
      abs_rel    = abs(rel_err),
      rel_err_pct = 100 * rel_err,
      abs_rel_pct = 100 * abs_rel
    ) %>%
    dplyr::select(dataset_id, parameter, estimate, true_value,
                  rel_err, rel_err_pct, abs_rel_pct, sq_rel,
                  exact_match, converged, cond_num_cor) %>%
    dplyr::arrange(parameter, dplyr::desc(abs_rel_pct))
}


## ---------------------------------------------------------------------------
## summarise_rse_contributions()
##
## Given the per-dataset table from extract_individual_rse(), report the
## top-k contributors to sum(sq_rel) for each parameter.  If one or two
## datasets account for > 50% of the sum, that's your outlier problem.
## ---------------------------------------------------------------------------
summarise_rse_contributions <- function(indiv_rse, top_k = 10) {
  indiv_rse %>%
    dplyr::group_by(parameter) %>%
    dplyr::mutate(
      total_sq   = sum(sq_rel, na.rm = TRUE),
      pct_of_sum = 100 * sq_rel / total_sq
    ) %>%
    dplyr::slice_max(order_by = sq_rel, n = top_k, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(parameter, dataset_id, estimate, true_value,
                  rel_err_pct, pct_of_sum, exact_match,
                  converged, cond_num_cor)
}





# ============================================================================
# Part 4: Aggregate HPCE output_N80 (scenario 16, N=80, full run) -------------
# ----------------------------------------------------------------------------
# Each output_N80/test_dsXXX.rds IS the bare `test` object (fields at top level:
# $label, $selected, $step_hist, $final_fit, $final_est, $rel_err, $diag,
# $parFixed, $cov_done, $runtime_sec). We wrap each into the pilot's per-ds
# shape so P3's OC helpers apply unchanged.
# ============================================================================

out_dir_N80 <- "output/output_N80/output_N80_sc16"   
ds_files    <- list.files(out_dir_N80, pattern = "^test_ds\\d+\\.rds$",
                          full.names = TRUE)

## ---- (P4.1) Load + wrap as pilot-shape list --------------------------------
scm_full_scn16_N80 <- purrr::map(ds_files, function(f) {
  tt    <- readRDS(f)
  ds_id <- as.integer(sub("^test_ds(\\d+)\\.rds$", "\\1", basename(f)))
  list(
    ds_id       = ds_id,
    ds_tag      = tt$label %||% sprintf("scn16_ds%d_N80", ds_id),
    t_base_sec  = NA_real_,             # not split in HPCE output
    t_scm_sec   = NA_real_,
    t_total_sec = tt$runtime_sec %||% NA_real_,
    test        = tt
  )
})

## ---- (P4.2) Truth for scenario 16 (reuse pilot's true_params) --------------
true_params <- readRDS("Inputdataset/true_params_long.rds")
true_rel_scn16 <- extract_true_relations(true_params, scenario_id = 16)
n_true_scn16   <- nrow(true_rel_scn16)

## ---- (P4.3) Per-dataset audit (mirrors P3.4) -------------------------------
oc_full_per_ds <- purrr::map_dfr(scm_full_scn16_N80, function(p) {
  if (is.null(p) || is.null(p$test)) return(NULL)
  tt  <- p$test
  mat <- match_selected_to_truth(tt$selected, true_rel_scn16)
  tibble::tibble(
    dataset_id   = p$ds_id,
    ds_tag       = p$ds_tag,
    t_base_min   = p$t_base_sec  / 60,
    t_scm_min    = p$t_scm_sec   / 60,
    t_total_min  = p$t_total_sec / 60,
    converged    = isTRUE(tt$diag$converged),
    cov_ok       = isTRUE(tt$diag$cov_ok),
    cond_num_cor = tt$diag$cond_num_cor,
    n_selected   = if (is.null(tt$selected)) NA_integer_ else nrow(tt$selected),
    n_true_hit   = mat$n_true_hit,
    n_false_pos  = mat$n_false_pos,
    exact_match  = mat$exact_match,
    selected     = list(tt$selected),
    final_est    = list(tt$final_est)
  )
})


## ---- (P4.4) Six OC artefacts (mirrors P3.5) --------------------------------
power_out_full <- compute_power_block(oc_full_per_ds,
                                      n_true     = n_true_scn16,
                                      cn_cor_cut = 1000)
oc_full_power    <- power_out_full$power
oc_full_relpower <- power_out_full$rel_power

## Pop params for RMRSE (TVKA dropped: fixed in model)
pop_params_rmrse <- c("TVCL", "TVVc", "TVQ", "TVVp",
                      "var_CL", "var_Vc", "cov_VcCL", "ResErr")
cov_params_rmrse <- c("CLBW", "CLcrCL", "VcBW", "VcSEX")

rmrse_out_full <- compute_rmrse_block(
  per_ds      = oc_full_per_ds,
  true_long   = true_params,
  scenario_id = 16,
  pop_params  = pop_params_rmrse,
  cov_params  = cov_params_rmrse
)
oc_full_rmrse_uncond <- rmrse_out_full$rmrse_uncond
oc_full_rmrse_cond   <- rmrse_out_full$rmrse_cond

oc_full_timing <- tibble::tibble(
  metric = c("n_datasets", "median_min_per_ds", "max_min_per_ds",
             "total_wall_min"),
  value  = c(
    nrow(oc_full_per_ds),
    stats::median(oc_full_per_ds$t_total_min, na.rm = TRUE),
    max(oc_full_per_ds$t_total_min,           na.rm = TRUE),
    sum(oc_full_per_ds$t_total_min,           na.rm = TRUE)
  )
)

##inspect RSE for each dataset in scenario 16
indiv_rse250_sc16 <- extract_individual_rse(
  per_ds      = oc_full_per_ds,
  true_long   = true_params,
  scenario_id = 16,
  params      = c("var_CL", "cov_VcCL")
)

indiv_rse_tail <- indiv_rse250_sc16  %>%
  dplyr::filter(abs_rel_pct >= 60)

ggplot(indiv_rse_tail, aes(x = abs_rel_pct, fill = parameter)) +
  geom_histogram(bins = 40, color = "white", alpha = 0.85) +
  facet_wrap(~ parameter, scales = "free", ncol = 1) +
  labs(
    title    = "Distribution of |Relative Error| per Dataset (Scenario 16, N = 80)",
    subtitle = "Unconditional; long right tail = a few datasets dominate RMRSE",
    x        = "|Relative Error| (%)",
    y        = "Number of datasets"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))


test_sc_16_ds203 <- readRDS("output/output_N80/output_N80_sc16/test_ds203.rds")



# ============================================================================
# Part 5: Full-grid aggregation — per (sample_size, scenario) OC summaries ----
# ----------------------------------------------------------------------------
# Directory layout:
#   output/output_N{NN}/output_N{NN}_sc{XX}/test_ds{DDD}.rds
# For every (N, scenario) cell we reuse the P3 helpers to compute:
#   - power block  (Power, PowerCN, PowerMinSuc)
#   - rel-power    (k = 1..n_true)
#   - RMRSE uncond + cond
#   - timing
# Results are stacked into long tibbles keyed by (sample_size, scenario).
# ============================================================================

root_dir <- "output"                       # parent of output_N40, output_N80, ...
sample_size_dirs <- list.dirs(root_dir, recursive = FALSE)
sample_size_dirs <- sample_size_dirs[grepl("output_N\\d+$", sample_size_dirs)]

## ---- (P5.1) One-cell driver: load + audit + OC for a single scenario dir ---
aggregate_one_cell <- function(sc_dir, N, scenario_id, true_long,
                               pop_params, cov_params, cn_cor_cut = 1000) {
  ds_files <- list.files(sc_dir, pattern = "^test_ds\\d+\\.rds$",
                         full.names = TRUE)
  if (length(ds_files) == 0L) return(NULL)

  true_rel <- extract_true_relations(true_long, scenario_id = scenario_id)
  n_true   <- nrow(true_rel)

  per_ds <- purrr::map_dfr(ds_files, function(f) {
    tt    <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(tt)) return(NULL)
    ds_id <- as.integer(sub("^test_ds(\\d+)\\.rds$", "\\1", basename(f)))
    mat   <- match_selected_to_truth(tt$selected, true_rel)
    tibble::tibble(
      dataset_id   = ds_id,
      ds_tag       = tt$label %||% sprintf("scn%02d_ds%d_N%d",
                                           scenario_id, ds_id, N),
      t_total_min  = (tt$runtime_sec %||% NA_real_) / 60,
      converged    = isTRUE(tt$diag$converged),
      cov_ok       = isTRUE(tt$diag$cov_ok),
      cond_num_cor = tt$diag$cond_num_cor,
      n_selected   = if (is.null(tt$selected)) NA_integer_ else nrow(tt$selected),
      n_true_hit   = mat$n_true_hit,
      n_false_pos  = mat$n_false_pos,
      exact_match  = mat$exact_match,
      selected     = list(tt$selected),
      final_est    = list(tt$final_est)
    )
  })
  if (nrow(per_ds) == 0L) return(NULL)

  pw <- compute_power_block(per_ds, n_true = n_true, cn_cor_cut = cn_cor_cut)
  rm <- compute_rmrse_block(per_ds, true_long, scenario_id,
                            pop_params, cov_params)

  timing <- tibble::tibble(
    metric = c("n_datasets", "median_min_per_ds",
               "max_min_per_ds", "total_wall_min"),
    value  = c(nrow(per_ds),
               stats::median(per_ds$t_total_min, na.rm = TRUE),
               max(per_ds$t_total_min,           na.rm = TRUE),
               sum(per_ds$t_total_min,           na.rm = TRUE))
  )

  list(per_ds       = per_ds,
       power        = pw$power,
       rel_power    = pw$rel_power,
       rmrse_uncond = rm$rmrse_uncond,
       rmrse_cond   = rm$rmrse_cond,
       timing       = timing,
       n_true       = n_true)
}

## ---- (P5.2) Walk the full grid ---------------------------------------------
oc_grid <- purrr::map_dfr(sample_size_dirs, function(N_dir) {
  N       <- as.integer(sub(".*output_N(\\d+)$", "\\1", N_dir))
  sc_dirs <- list.dirs(N_dir, recursive = FALSE)
  sc_dirs <- sc_dirs[grepl(sprintf("output_N%d_sc\\d+$", N), sc_dirs)]

  if (length(sc_dirs) == 0L) {
    message(sprintf("[N=%d] no scenario dirs yet — skipping", N))
    return(NULL)
  }

  cell_N <- purrr::map_dfr(sc_dirs, function(sc_dir) {
    scenario_id <- as.integer(sub(".*_sc(\\d+)$", "\\1", sc_dir))
    message(sprintf("[N=%d scn=%02d] %s", N, scenario_id, sc_dir))
    cell <- aggregate_one_cell(sc_dir, N, scenario_id,
                               true_params, pop_params_rmrse, cov_params_rmrse)
    if (is.null(cell)) {
      message(sprintf("[N=%d scn=%02d] empty — skipping", N, scenario_id))
      return(NULL)
    }
    tibble::tibble(
      sample_size  = N, scenario = scenario_id,
      n_datasets   = nrow(cell$per_ds), n_true = cell$n_true,
      power        = list(cell$power),
      rel_power    = list(cell$rel_power),
      rmrse_uncond = list(cell$rmrse_uncond),
      rmrse_cond   = list(cell$rmrse_cond),
      timing       = list(cell$timing),
      per_ds       = list(cell$per_ds)
    )
  })
  saveRDS(cell_N, file.path(root_dir, sprintf("oc_grid_N%d.rds", N)))
  cell_N
})
## ---- (P5.3) Flat long tibbles for plotting / export ------------------------
oc_power_long <- oc_grid %>%
  dplyr::select(sample_size, scenario, power) %>%
  tidyr::unnest(power)

oc_relpower_long <- oc_grid %>%
  dplyr::select(sample_size, scenario, rel_power) %>%
  tidyr::unnest(rel_power)

oc_rmrse_uncond_long <- oc_grid %>%
  dplyr::select(sample_size, scenario, rmrse_uncond) %>%
  tidyr::unnest(rmrse_uncond)

oc_rmrse_cond_long <- oc_grid %>%
  dplyr::select(sample_size, scenario, rmrse_cond) %>%
  tidyr::unnest(rmrse_cond)

oc_timing_long <- oc_grid %>%
  dplyr::select(sample_size, scenario, timing) %>%
  tidyr::unnest(timing)

## ---- (P5.4) Persist --------------------------------------------------------
saveRDS(oc_grid,               file.path(root_dir, "oc_grid_nested.rds"))
saveRDS(oc_power_long,         file.path(root_dir, "oc_power_long.rds"))
saveRDS(oc_relpower_long,      file.path(root_dir, "oc_relpower_long.rds"))
saveRDS(oc_rmrse_uncond_long,  file.path(root_dir, "oc_rmrse_uncond_long.rds"))
saveRDS(oc_rmrse_cond_long,    file.path(root_dir, "oc_rmrse_cond_long.rds"))
saveRDS(oc_timing_long,        file.path(root_dir, "oc_timing_long.rds"))

message(sprintf("[done] %d cells aggregated across %d sample sizes",
                nrow(oc_grid), dplyr::n_distinct(oc_grid$sample_size)))
print(oc_power_long, n = Inf)







#### Visualization:1.1  [Power] ----------------------
library(ggplot2)
##Aggregated graph with N=40/80/300 in the same graph
library(scales)
## ---- Prep: filter to the three sample sizes and pick a metric --------------
## The reference figure (Fig. 8) shows a SINGLE power definition per line.
## We'll default to "Power" (exact-match, unconditional).  Change `which_metric`
## to "PowerCN" or "PowerMinSuc" if you want the conditional versions.
which_metric <- "Power"

plot_df_all <- oc_power_long %>%
  dplyr::filter(sample_size %in% c(40, 80, 300),
                metric      == which_metric) %>%
  dplyr::mutate(
    scenario    = as.integer(scenario),
    sample_size = factor(sample_size,
                         levels = c(300, 80, 40),
                         labels = c("300 subj", "80 subj", "40 subj"))
  )

## ---- Plot ------------------------------------------------------------------
p_power_by_N <- ggplot(plot_df_all,
                       aes(x = scenario, y = value,
                           color = sample_size,
                           shape = sample_size,
                           group = sample_size)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey40") +
  scale_x_continuous(breaks = 1:16) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_color_manual(
    values = c("300 subj" = "#1f77b4",   # blue
               "80 subj"  = "#ff7f0e",   # orange
               "40 subj"  = "#7f7f7f")   # grey
  ) +
  scale_shape_manual(
    values = c("300 subj" = 16, "80 subj" = 16, "40 subj" = 16)
  ) +
  labs(
    title    = sprintf("Power by Scenario and Sample Size (metric = %s)",
                       which_metric),
    subtitle = "Dashed line = 80% power threshold",
    x        = "Simulation Scenario",
    y        = "Power (%)",
    color    = NULL,
    shape    = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "top",
    legend.direction = "horizontal",
    panel.grid.minor = element_blank()
  )

print(p_power_by_N)
ggsave("figs/power_by_scenario_by_N.png", p_power_by_N,
       width = 9, height = 5, dpi = 150)













#### Visualization: 1.2 [RMRSE tables] ----------------------
## ---------------------------------------------------------------------------
## Config: which parameters to show, and in what order + display name
## Mirrors Table 2 of the reference paper.  KA is dropped because TVKA is
## fixed in the model (not estimated).
## ---------------------------------------------------------------------------
pop_param_order <- c(
  "TVVc"     = "Vc (%)",
  "TVCL"     = "CL (%)",
  "TVVp"     = "Vp (%)",
  "TVQ"      = "Q (%)",
  "ResErr"   = "Err (%)",
  "var_Vc"   = "\u03C91-Vc (%)",     # ω1-Vc
  "cov_VcCL" = "cov (Vc,CL) (%)",
  "var_CL"   = "\u03C92-CL (%)"      # ω2-CL
)

## ---------------------------------------------------------------------------
## Helper: pivot one RMRSE long tibble (uncond or cond) to Table-2 wide form
## ---------------------------------------------------------------------------
make_rmrse_table <- function(rmrse_long, oc_grid, N_target = 300,
                             param_map   = pop_param_order,
                             digits      = 1,
                             per_param_n = FALSE) {   # <- NEW

  ## n_relations per scenario (from oc_grid$n_true)
  n_rel_tbl <- oc_grid %>%
    dplyr::filter(sample_size == N_target) %>%
    dplyr::select(scenario, n_relations = n_true)

  ## --- pull rows for this N and the parameters we display -------------
  base <- rmrse_long %>%
    dplyr::filter(sample_size == N_target,
                  parameter   %in% names(param_map)) %>%
    dplyr::mutate(RMRSE_pct = round(RMRSE_pct, digits))

  ## --- wide RMRSE table (same as before) ------------------------------
  wide_rmrse <- base %>%
    dplyr::select(scenario, parameter, RMRSE_pct) %>%
    tidyr::pivot_wider(names_from = parameter, values_from = RMRSE_pct)

  for (p in names(param_map)) {
    if (!p %in% names(wide_rmrse)) wide_rmrse[[p]] <- NA_real_
  }
  wide_rmrse <- wide_rmrse[, c("scenario", names(param_map))]

  ## --- n_used --------------------------------------------------------
  if (per_param_n) {
    ## one n_used column per parameter (…_n suffix)
    wide_n <- base %>%
      dplyr::select(scenario, parameter, n_used) %>%
      tidyr::pivot_wider(names_from  = parameter,
                         values_from = n_used,
                         names_glue  = "{parameter}_n")
    for (p in names(param_map)) {
      col <- paste0(p, "_n")
      if (!col %in% names(wide_n)) wide_n[[col]] <- 0L
    }
    wide_n <- wide_n[, c("scenario", paste0(names(param_map), "_n"))]

    out <- n_rel_tbl %>%
      dplyr::left_join(wide_rmrse, by = "scenario") %>%
      dplyr::left_join(wide_n,     by = "scenario") %>%
      dplyr::arrange(scenario)

    ## rename RMRSE cols to display headers; keep _n cols as-is
    names(out)[match(names(param_map), names(out))] <- unname(param_map)
    names(out)[1:2] <- c("Scenario", "n relations")

  } else {
    ## single scenario-level "n datasets" column
    ## = max n_used across the displayed parameters (ignores structural NAs)
    n_used_tbl <- base %>%
      dplyr::group_by(scenario) %>%
      dplyr::summarise(`n datasets` = suppressWarnings(max(n_used, na.rm = TRUE)),
                       .groups = "drop") %>%
      dplyr::mutate(`n datasets` = ifelse(is.finite(`n datasets`),
                                          `n datasets`, NA_integer_))

    out <- n_rel_tbl %>%
      dplyr::left_join(n_used_tbl, by = "scenario") %>%
      dplyr::left_join(wide_rmrse, by = "scenario") %>%
      dplyr::arrange(scenario)

    names(out)[match(names(param_map), names(out))] <- unname(param_map)
    names(out)[1:3] <- c("Scenario", "n relations", "n datasets")
  }

  out
}

## ---------------------------------------------------------------------------
## Build the six tables (uncond + cond) x (N = 40, 80, 300)
## ---------------------------------------------------------------------------
table2_uncond_N300 <- make_rmrse_table(oc_rmrse_uncond_long, oc_grid, N_target = 300)
table2_cond_N300   <- make_rmrse_table(oc_rmrse_cond_long,   oc_grid, N_target = 300)
table2_uncond_N80  <- make_rmrse_table(oc_rmrse_uncond_long, oc_grid, N_target = 80)
table2_cond_N80    <- make_rmrse_table(oc_rmrse_cond_long,   oc_grid, N_target = 80)
table2_uncond_N40  <- make_rmrse_table(oc_rmrse_uncond_long, oc_grid, N_target = 40)
table2_cond_N40    <- make_rmrse_table(oc_rmrse_cond_long,   oc_grid, N_target = 40)

## Bundle into a named list -- makes downstream export & printing DRY
table2_all <- list(
  "RMRSE_uncond_N40"  = table2_uncond_N40,
  "RMRSE_cond_N40"    = table2_cond_N40,
  "RMRSE_uncond_N80"  = table2_uncond_N80,
  "RMRSE_cond_N80"    = table2_cond_N80,
  "RMRSE_uncond_N300" = table2_uncond_N300,
  "RMRSE_cond_N300"   = table2_cond_N300
)

## ---------------------------------------------------------------------------
## Save each table to its own CSV
## ---------------------------------------------------------------------------
dir.create("tables", showWarnings = FALSE, recursive = TRUE)
purrr::iwalk(table2_all, function(tbl, nm) {
  readr::write_csv(tbl, file.path("tables", sprintf("table2_%s.csv", nm)))
})

## ---------------------------------------------------------------------------
## Save all six as sheets of a single XLSX
## ---------------------------------------------------------------------------
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(
    table2_all,
    path = "tables/table2_RMRSE_all_N.xlsx"
  )
  message("Wrote tables/table2_RMRSE_all_N.xlsx with ",
          length(table2_all), " sheets.")
}

## ---------------------------------------------------------------------------
## Pretty-print all six with knitr::kable
## ---------------------------------------------------------------------------
if (requireNamespace("knitr", quietly = TRUE)) {

  ## Human-friendly caption per key
  caption_lookup <- c(
    "RMRSE_uncond_N40"  = "Table 2a. RMRSE (unconditional), N = 40",
    "RMRSE_cond_N40"    = "Table 2b. RMRSE (conditional on exact_match), N = 40",
    "RMRSE_uncond_N80"  = "Table 2c. RMRSE (unconditional), N = 80",
    "RMRSE_cond_N80"    = "Table 2d. RMRSE (conditional on exact_match), N = 80",
    "RMRSE_uncond_N300" = "Table 2e. RMRSE (unconditional), N = 300",
    "RMRSE_cond_N300"   = "Table 2f. RMRSE (conditional on exact_match), N = 300"
  )

  purrr::iwalk(table2_all, function(tbl, nm) {
    cat(sprintf("\n=== %s ===\n", caption_lookup[[nm]]))
    print(knitr::kable(tbl, format = "simple",
                       caption = caption_lookup[[nm]]))
  })
}

## ---------------------------------------------------------------------------
## (Optional) One combined "long" table with an N column,
## useful for a single-sheet export or ggplot heatmaps across N.
## ---------------------------------------------------------------------------
table2_combined <- purrr::imap_dfr(table2_all, function(tbl, nm) {
  tbl %>%
    dplyr::mutate(
      sample_size = as.integer(sub(".*_N(\\d+)$", "\\1", nm)),
      condition   = ifelse(grepl("uncond", nm), "unconditional", "conditional"),
      .before     = 1
    )
})
readr::write_csv(table2_combined, "tables/table2_RMRSE_combined_long.csv")
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(
    list("RMRSE_all"        = table2_combined,
         "RMRSE_uncond_N40" = table2_uncond_N40,
         "RMRSE_cond_N40"   = table2_cond_N40,
         "RMRSE_uncond_N80" = table2_uncond_N80,
         "RMRSE_cond_N80"   = table2_cond_N80,
         "RMRSE_uncond_N300"= table2_uncond_N300,
         "RMRSE_cond_N300"  = table2_cond_N300),
    path = "tables/table2_RMRSE_all_N_plus_combined.xlsx"
  )
}






## visualization convergenve rate----------------------- 
library(dplyr)
library(tidyr)
library(purrr)

## ---------------------------------------------------------------------------
## Helper: summarise one per_ds tibble
## ---------------------------------------------------------------------------
summarise_per_ds <- function(per_ds) {

  if (is.null(per_ds) || !nrow(per_ds))
    return(tibble(n_ds = 0L, n_conv = 0L, conv_rate = NA_real_,
                  cn_median = NA_real_, cn_q1 = NA_real_,
                  cn_q3 = NA_real_, cn_iqr = NA_real_,
                  cn_n = 0L))

  ## be robust to alternate column names
  cn_col <- intersect(c("cond_num_cor"),
                      names(per_ds))[1]
  cn_vec <- if (!is.na(cn_col)) per_ds[[cn_col]] else rep(NA_real_, nrow(per_ds))

  ## condition number IQR computed on converged runs only
  conv_lgl <- as.logical(per_ds$converged)
  cn_conv  <- cn_vec[conv_lgl & is.finite(cn_vec)]

  qs <- if (length(cn_conv))
          stats::quantile(cn_conv, c(.5, .75, 1), na.rm = TRUE, names = FALSE)
        else c(NA_real_, NA_real_, NA_real_)

  tibble(
    n_ds      = nrow(per_ds),
    n_conv    = sum(conv_lgl, na.rm = TRUE),
    conv_rate = mean(conv_lgl, na.rm = TRUE),
    cn_n      = length(cn_conv),
    cn_median = qs[1],
    cn_q3     = qs[2],
    cn_max    = qs[3]
  )
}

## ---------------------------------------------------------------------------
## Apply across oc_grid
## ---------------------------------------------------------------------------
oc_convergence_cn <- oc_grid %>%
  mutate(.summ = purrr::map(per_ds, summarise_per_ds)) %>%
  select(sample_size, scenario, .summ) %>%
  tidyr::unnest(.summ) %>%
  arrange(sample_size, scenario)

