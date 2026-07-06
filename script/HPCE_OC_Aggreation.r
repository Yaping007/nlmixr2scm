library(testthat)
library(nlmixr2scm)
library(nlmixr2utils)
devtools::load_all("C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm")
library(nlmixr2)
library(rxode2)
library(tidyverse)

library(devtools)
library(haven)
test_ds02_m02 <- readRDS("simulated_virtual_dataset_eta_filtered_N80/stage1_pilot_scn16_N80_m02/test_full_fast_ds02.rds")
test_ds02 <- readRDS("simulated_virtual_dataset_eta_filtered_N80/stage1_pilot_scn16_N80/test_full_fast_ds02.rds")
test_ds03 <- readRDS("simulated_virtual_dataset_eta_filtered_N80/stage1_pilot_scn16_N80/test_full_fast_ds03.rds")
test_ds13 <- readRDS("output_N80/test_ds13.rds")


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
  rel %>%
    dplyr::mutate(shape = ifelse(grepl("^[0-9]+$", shape), "cat", shape)) %>%
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


# ============================================================================
# Part 4: Aggregate HPCE output_N80 (scenario 16, N=80, full run) -------------
# ----------------------------------------------------------------------------
# Each output_N80/test_dsXXX.rds IS the bare `test` object (fields at top level:
# $label, $selected, $step_hist, $final_fit, $final_est, $rel_err, $diag,
# $parFixed, $cov_done, $runtime_sec). We wrap each into the pilot's per-ds
# shape so P3's OC helpers apply unchanged.
# ============================================================================

out_dir_N80 <- "output_N80"
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
true_params <- readRDS("simulated_virtual_dataset/true_params_long.rds")
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

## ---- (P4.5) Quick console summary ------------------------------------------
message(sprintf("[output_N80] loaded %d datasets", nrow(oc_full_per_ds)))
print(oc_full_power)
print(oc_full_relpower)
print(oc_full_rmrse_uncond)
print(oc_full_rmrse_cond)
print(oc_full_timing)





power_out <- compute_power_block(oc_pilot_per_ds,
                                 n_true      = n_true_scn16,
                                 cn_cor_cut  = 1000)
oc_pilot_power    <- power_out$power
oc_pilot_relpower <- power_out$rel_power

rmrse_out <- compute_rmrse_block(
  per_ds      = oc_pilot_per_ds,
  true_long   = true_params,
  scenario_id = 16,
  pop_params  = pop_params_rmrse,
  cov_params  = cov_params_rmrse
)
oc_pilot_rmrse_uncond <- rmrse_out$rmrse_uncond
oc_pilot_rmrse_cond   <- rmrse_out$rmrse_cond

## Timing summary + scale-up projection
oc_pilot_timing <- tibble::tibble(
  metric = c("median_min_per_ds", "max_min_per_ds", "total_pilot_min",
             "projected_hours_100x16_seq",
             "projected_hours_100x16_workers4"),
  value  = c(
    stats::median(oc_pilot_per_ds$t_total_min, na.rm = TRUE),
    max(oc_pilot_per_ds$t_total_min,           na.rm = TRUE),
    sum(oc_pilot_per_ds$t_total_min,           na.rm = TRUE),
    stats::median(oc_pilot_per_ds$t_total_min, na.rm = TRUE) * 100 * 16 / 60,
    stats::median(oc_pilot_per_ds$t_total_min, na.rm = TRUE) * 100 * 16 / 60 / 4
  )
)
