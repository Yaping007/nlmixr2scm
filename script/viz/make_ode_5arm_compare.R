# ============================================================================
# make_ode_5arm_compare.R
# Five-arm comparison for the ODE structural model, N = 80, competing
# (power + lin/exp + cat) covariate search space:
#
#   1. bobyqa   warm-up ON   output/scm_focei_bobyqa_703est_0805_aggregated
#   2. bobyqa   warm-up OFF  output/scm_profile_off_bobyqa_aggregated
#   3. lbfgsb3c warm-up ON   output/lbfgsb3c_est703_08132026_warmon_aggregated
#   4. lbfgsb3c warm-up OFF  output/lbfgsb3c_est703_08132026_warmoff_aggregated
#   5. PsN-SCM  (NONMEM)     output/psn_scm_combined_aggregated
#
# Figures written to output/figures/ode_5arm/:
#   fig_5arm_convergence_N80.{png,pdf}   Converged% vs scenario (nlmixr2 arms)
#   fig_5arm_cn1000_N80.{png,pdf}        % runs with CN < 1000 (nlmixr2 arms)
#   fig_5arm_power_N80.{png,pdf}         Selection power (all 5 arms)
#   fig_5arm_error_N80.{png,pdf}         MARE% per parameter, scenario 16
#   fig_5arm_covsel_N80.{png,pdf}        FN/FP selection-pattern heatmap (x5)
#
# DATA CAVEATS handled here:
#   * lbfgsb3c warm-OFF scm_covsel_by_covar.csv is EMPTY -> the covsel panel
#     for that arm is derived from scm_covsel_long.csv.
#   * PsN ODE rows have no covariance step: Converged_pct / CNBelowCutoff_pct
#     are 0/NA (missing, not real). PsN is therefore EXCLUDED from the
#     convergence and CN panels (noted in the subtitle).
#
# Usage:  source("script/viz/make_ode_5arm_compare.R")
# ============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

SN     <- 80
STRUCT <- "ode"
OUT    <- "output/figures/ode_5arm"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- arm registry ----------------------------------------------------------
ARMS <- tibble::tribble(
  ~cond,          ~dir,                                                 ~estimator,   ~outer_opt,
  "bobyqa ON",    "output/scm_focei_bobyqa_703est_0805_aggregated",     "focei",      "bobyqa",
  "bobyqa OFF",   "output/scm_profile_off_bobyqa_aggregated",           "focei",      "bobyqa",
  "lbfgsb3c ON",  "output/lbfgsb3c_est703_08132026_warmon_aggregated",  "focei",      "lbfgsb3c",
  "lbfgsb3c OFF", "output/lbfgsb3c_est703_08132026_warmoff_aggregated", "focei",      "lbfgsb3c",
  "PsN-SCM",      "output/psn_scm_combined_aggregated",                 "nonmem_scm", "focei"
)
COND_LVL <- ARMS$cond

## colour: teal = bobyqa, salmon = lbfgsb3c, black = PsN; light = warm OFF
PAL <- c("bobyqa ON"   = "#0F8B8D", "bobyqa OFF"   = "#8EC1C2",
         "lbfgsb3c ON" = "#C1666B", "lbfgsb3c OFF" = "#E3A9AC",
         "PsN-SCM"     = "#2B2B2B")
## linetype: warm ON solid, warm OFF dashed
LTY <- c("bobyqa ON" = "solid", "bobyqa OFF" = "22",
         "lbfgsb3c ON" = "solid", "lbfgsb3c OFF" = "22", "PsN-SCM" = "solid")

.read0 <- function(p) {
  if (!file.exists(p) || file.info(p)$size == 0) return(NULL)
  suppressWarnings(readr::read_csv(p, show_col_types = FALSE))
}

## generic per-arm metric loader (filters ode / N=80 / arm's est+opt) ---------
load_metric <- function(csv) {
  bind_rows(lapply(seq_len(nrow(ARMS)), function(i) {
    a <- ARMS[i, ]
    d <- .read0(file.path(a$dir, csv))
    if (is.null(d) || !"estimator" %in% names(d)) return(NULL)
    d |>
      filter(estimator == a$estimator, outer_opt == a$outer_opt,
             sample_N == SN, structure == STRUCT) |>
      mutate(cond = a$cond)
  })) |>
    mutate(cond = factor(cond, levels = COND_LVL))
}

theme_5 <- function() {
  theme_bw(base_size = 14) +
    theme(legend.position = "top", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))
}

line5 <- function(df, ycol, ylab, title, sub = NULL) {
  ggplot(df, aes(scenario, .data[[ycol]], colour = cond,
                 linetype = cond, group = cond)) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
    scale_colour_manual(values = PAL, name = NULL) +
    scale_linetype_manual(values = LTY, name = NULL) +
    scale_x_continuous(breaks = 1:16) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(x = "Simulation scenario", y = ylab, title = title, subtitle = sub) +
    theme_5()
}

## ============================================================================
## 1. Convergence rate  (nlmixr2 arms only)
## ============================================================================
diag <- load_metric("scm_diag_rates.csv")

p_conv <- line5(filter(diag, cond != "PsN-SCM"), "Converged_pct",
                "Converged runs", "Convergence rate  (N = 80, ODE)",
                "PsN ODE ($DES) reports no convergence flag - excluded")
ggsave(file.path(OUT, "fig_5arm_convergence_N80.png"), p_conv, width = 9, height = 5.5, dpi = 200)
ggsave(file.path(OUT, "fig_5arm_convergence_N80.pdf"), p_conv, width = 9, height = 5.5)

## ============================================================================
## 2. Condition number < 1000  (nlmixr2 arms only)
## ============================================================================
p_cn <- line5(filter(diag, cond != "PsN-SCM"), "CNBelowCutoff_pct",
              "Runs with condition number < 1000",
              "Well-conditioned covariance step  (N = 80, ODE)",
              "PsN ODE ($DES) reports no condition number - excluded")
ggsave(file.path(OUT, "fig_5arm_cn1000_N80.png"), p_cn, width = 9, height = 5.5, dpi = 200)
ggsave(file.path(OUT, "fig_5arm_cn1000_N80.pdf"), p_cn, width = 9, height = 5.5)

## ============================================================================
## 3. Selection power  (all 5 arms)
## ============================================================================
pw <- load_metric("scm_power.csv") |> mutate(Power = Power * 100)
p_pow <- line5(pw, "Power", "Power (true covariate detected)",
               "Selection power  (N = 80, ODE)")
ggsave(file.path(OUT, "fig_5arm_power_N80.png"), p_pow, width = 9, height = 5.5, dpi = 200)
ggsave(file.path(OUT, "fig_5arm_power_N80.pdf"), p_pow, width = 9, height = 5.5)

## ============================================================================
## 4. Estimation accuracy  (MARE%, scenario 16)
## ============================================================================
.PARAM_LAB <- c(
  CLBW = "CL~BW", CLcrCL = "CL~CrCL", VcBW = "Vc~BW", VcSEX = "Vc~SEX",
  TVCL = "TVCL", TVVc = "TVVc", TVQ = "TVQ", TVVp = "TVVp", TVKA = "TVKA",
  var_CL = "var(CL)", var_Vc = "var(Vc)", cov_VcCL = "cov(Vc,CL)", ResErr = "ResErr")
.PARAM_ORD <- c("CLBW","CLcrCL","VcBW","VcSEX","TVCL","TVVc",
                "TVQ","TVVp","var_CL","var_Vc","cov_VcCL","ResErr","TVKA")

err <- load_metric("scm_estim_all.csv")
if (nrow(err) && "MARE_pct" %in% names(err)) {
  err <- err |>
    filter(scenario == 16, parameter %in% .PARAM_ORD) |>
    mutate(param_lab = factor(.PARAM_LAB[parameter],
                              levels = rev(.PARAM_LAB[.PARAM_ORD])),
           value = pmin(MARE_pct, 70))

  p_err <- ggplot(err, aes(value, param_lab, colour = cond)) +
    geom_point(position = position_dodge(width = 0.7), size = 2.4) +
    scale_colour_manual(values = PAL, name = NULL) +
    scale_x_continuous(labels = function(x) paste0(x, "%")) +
    labs(x = "MARE (%)", y = "Parameter",
         title = "Estimation accuracy, scenario 16  (N = 80, ODE; lower = better)") +
    theme_5()
  ggsave(file.path(OUT, "fig_5arm_error_N80.png"), p_err, width = 9, height = 6.5, dpi = 200)
  ggsave(file.path(OUT, "fig_5arm_error_N80.pdf"), p_err, width = 9, height = 6.5)
} else {
  warning("No MARE data for scenario 16 - estimation-accuracy panel skipped.")
}

## ============================================================================
## 5. Selection-pattern heatmap  (5 arms; lbfgsb3c OFF derived from long form)
## ============================================================================
denom <- diag |> select(cond, scenario, n_datasets = n_total)

build_covsel <- function(i) {
  a <- ARMS[i, ]
  bycov <- .read0(file.path(a$dir, "scm_covsel_by_covar.csv"))
  if (!is.null(bycov) && "estimator" %in% names(bycov) && nrow(bycov)) {
    return(bycov |>
      filter(estimator == a$estimator, outer_opt == a$outer_opt,
             sample_N == SN, structure == STRUCT) |>
      transmute(cond = a$cond, scenario, var, covar, shape,
                is_true, n_TP, n_FN, n_FP, n_datasets))
  }
  ## fallback: derive per-covariate counts from the long form
  lng <- .read0(file.path(a$dir, "scm_covsel_long.csv"))
  if (is.null(lng) || !"verdict" %in% names(lng)) {
    warning("No covsel data for '", a$cond, "' - skipped"); return(NULL)
  }
  lng |>
    filter(estimator == a$estimator, outer_opt == a$outer_opt,
           sample_N == SN, structure == STRUCT) |>
    group_by(scenario, var, covar, shape) |>
    summarise(is_true = any(in_true),
              n_TP = sum(verdict == "TP"),
              n_FP = sum(verdict == "FP"),
              n_FN = sum(verdict == "FN"), .groups = "drop") |>
    mutate(cond = a$cond) |>
    left_join(filter(denom, cond == a$cond) |> select(scenario, n_datasets),
              by = "scenario") |>
    transmute(cond, scenario, var, covar, shape,
              is_true, n_TP, n_FN, n_FP, n_datasets)
}

cov_raw <- bind_rows(lapply(seq_len(nrow(ARMS)), build_covsel))

if (nrow(cov_raw)) {
  ## PsN encodes the linear shape as "exp"; fold to "lin" for a common y-axis.
  ## Normalise var/covar casing so arms agree on one effect label.
  cov_raw <- cov_raw |>
    mutate(shape = ifelse(tolower(shape) == "exp", "lin", tolower(shape)),
           var   = tolower(var), covar = toupper(covar),
           n_datasets = ifelse(is.na(n_datasets) | n_datasets == 0, 100, n_datasets),
           err_rate   = ifelse(is_true, n_FN / n_datasets, n_FP / n_datasets),
           err_signed = ifelse(is_true, -err_rate, err_rate) * 100,
           effect     = paste0(toupper(var), "~", covar, ".", shape),
           scenario   = factor(scenario, levels = sort(unique(scenario))),
           cond       = factor(cond, levels = COND_LVL))

  eff_lvls <- cov_raw |> distinct(effect, var, covar, shape) |>
    arrange(var, covar, shape) |> distinct(effect) |> pull(effect)
  cov_raw <- mutate(cov_raw, effect = factor(effect, levels = rev(eff_lvls)))

  p_cov <- ggplot(cov_raw, aes(scenario, effect, fill = err_signed)) +
    geom_tile(colour = "grey85", linewidth = 0.25) +
    geom_tile(data = filter(cov_raw, is_true), fill = NA, colour = "black",
              linewidth = 0.2) +
    facet_wrap(~ cond, nrow = 1) +
    scale_fill_gradient2(low = "#08519C", mid = "white", high = "#B30000",
      midpoint = 0, limits = c(-100, 100), breaks = seq(-100, 100, 50),
      labels = c("FN 100%","FN 50%","0","FP 50%","FP 100%"), name = NULL) +
    labs(x = "Simulation scenario",
         y = "Covariate effect (param ~ covariate . shape)",
         title = "Covariate-selection error pattern  (N = 80, ODE)",
         subtitle = "Black outline = true effect (FN rate); others = FP rate") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), legend.position = "top",
          plot.title = element_text(face = "bold"),
          axis.text.x = element_text(size = 7)) +
    guides(fill = guide_colourbar(barwidth = 16))
  ggsave(file.path(OUT, "fig_5arm_covsel_N80.png"), p_cov, width = 16, height = 9, dpi = 200)
  ggsave(file.path(OUT, "fig_5arm_covsel_N80.pdf"), p_cov, width = 16, height = 9)
} else {
  warning("No covariate-selection data - covsel heatmap skipped.")
}

message("\nFive-arm ODE (N = 80) figures written under ", OUT)
