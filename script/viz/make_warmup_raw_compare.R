# ============================================================================
# make_warmup_raw_compare.R
# RAW (not difference) side-by-side of the SCM warm-start (Mechanism B,
# profileInitOnStall) ON vs OFF, with linCmt AND ode in ONE figure per metric.
# All panels are N = 80, focei / bobyqa.
#
#   warmup-on  : output/scm_focei_bobyqa_703est_0805_aggregated
#   warmup-off : output/scm_profile_off_aggregated   (now includes ode)
#
# Four combined figures (output/figures/warmup_compare/raw/):
#   fig_raw_convergence_N80.{png,pdf}  Converged% / CN-below% vs scenario
#   fig_raw_power_N80.{png,pdf}        Power / PowerCN / PowerMinSuc vs scenario
#   fig_raw_error_N80.{png,pdf}        MARE% per parameter, scenario 16
#   fig_raw_covsel_N80.{png,pdf}       FN/FP selection-pattern heatmaps (2x2)
#
# Usage:  source("script/viz/make_warmup_raw_compare.R")
# ============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

## helpers for the covsel heatmap (.fix_covsel_denom, .VAR_LAB, .COVAR_ORD,
## .SHAPE_ORD, .SHAPE_LAB, theme_scm, %||%)
source("script/viz/fig_covsel_heatmap.R")

ON_DIR  <- "output/scm_focei_bobyqa_703est_0805_aggregated"
OFF_DIR <- "output/scm_profile_off_aggregated"
PSN_DIR <- "output/psn_scm_combined_aggregated"
OUT     <- "output/figures/warmup_compare/raw"
EST <- "focei"; OUT_OPT <- "bobyqa"; SN <- 80
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

PAL <- c("warm-up ON" = "#0F8B8D", "warm-up OFF" = "#C1666B",
         "PsN-SCM" = "#2B2B2B")
STRUCT_LAB <- c(linCmt = "linCmt (analytic)", ode = "ODE")

## load one metric CSV from both folders, tag condition, keep N=80 focei/bobyqa
load_both <- function(csv) {
  rd <- function(dir, cond) read_csv(file.path(dir, csv), show_col_types = FALSE) |>
    filter(estimator == EST, outer_opt == OUT_OPT, sample_N == SN,
           structure %in% c("linCmt", "ode")) |>
    mutate(cond = cond)
  bind_rows(rd(ON_DIR, "warm-up ON"), rd(OFF_DIR, "warm-up OFF")) |>
    mutate(cond      = factor(cond, levels = names(PAL)),
           struct_f  = factor(STRUCT_LAB[structure],
                              levels = unname(STRUCT_LAB)))
}

## ---- 1. convergence --------------------------------------------------------
diag <- load_both("scm_diag_rates.csv") |>
  select(scenario, cond, struct_f,
         Converged = Converged_pct, `CN below cutoff` = CNBelowCutoff_pct) |>
  pivot_longer(c(Converged, `CN below cutoff`), names_to = "metric",
               values_to = "pct") |>
  mutate(metric = factor(metric, levels = c("Converged", "CN below cutoff")))

p_conv <- ggplot(diag, aes(scenario, pct, colour = cond, group = cond)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  facet_grid(metric ~ struct_f) +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_x_continuous(breaks = 1:16) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Simulation scenario", y = "Run-quality rate",
       title = "Convergence & covariance step  (N = 80)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUT, "fig_raw_convergence_N80.png"), p_conv, width = 11, height = 6.5, dpi = 200)
ggsave(file.path(OUT, "fig_raw_convergence_N80.pdf"), p_conv, width = 11, height = 6.5)

## ---- 2. power (Power only; warm-up ON/OFF + PsN-SCM) -----------------------
## PsN uses estimator="nonmem_scm"/outer_opt="focei" and its analytic structure
## is "advan4" (folded to the linCmt analytic panel); ode maps to ODE.
STRUCT_MAP <- c(linCmt = "linCmt (analytic)", advan4 = "linCmt (analytic)",
                ode = "ODE")

read_power <- function(dir, cond, est, oopt) {
  read_csv(file.path(dir, "scm_power.csv"), show_col_types = FALSE) |>
    filter(estimator == est, outer_opt == oopt, sample_N == SN,
           structure %in% c("linCmt", "advan4", "ode")) |>
    transmute(scenario, Power = Power * 100,
              struct_f = factor(STRUCT_MAP[structure],
                                levels = unname(STRUCT_LAB)),
              cond = cond)
}

pw <- bind_rows(
  read_power(ON_DIR,  "warm-up ON",  EST, OUT_OPT),
  read_power(OFF_DIR, "warm-up OFF", EST, OUT_OPT),
  read_power(PSN_DIR, "PsN-SCM",     "nonmem_scm", "focei")) |>
  mutate(cond = factor(cond, levels = names(PAL)))

p_pow <- ggplot(pw, aes(scenario, Power, colour = cond, group = cond)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  facet_wrap(~ struct_f, nrow = 1) +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_x_continuous(breaks = 1:16) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Simulation scenario", y = "Power (true covariate detected)",
       title = "Selection power  (N = 80)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUT, "fig_raw_power_N80.png"), p_pow, width = 11, height = 5.5, dpi = 200)
ggsave(file.path(OUT, "fig_raw_power_N80.pdf"), p_pow, width = 11, height = 5.5)

## ---- 3. estimation accuracy (MARE%, scenario 16, unconditioned) ------------
.PARAM_LAB <- c(
  CLBW = "CL~BW", CLcrCL = "CL~CrCL", VcBW = "Vc~BW", VcSEX = "Vc~SEX",
  TVCL = "TVCL", TVVc = "TVVc", TVQ = "TVQ", TVVp = "TVVp", TVKA = "TVKA",
  var_CL = "var(CL)", var_Vc = "var(Vc)", cov_VcCL = "cov(Vc,CL)", ResErr = "ResErr")
.PARAM_ORD <- c("CLBW","CLcrCL","VcBW","VcSEX","TVCL","TVVc",
                "TVQ","TVVp","var_CL","var_Vc","cov_VcCL","ResErr","TVKA")

err <- load_both("scm_estim_all.csv") |>
  filter(scenario == 16, parameter %in% .PARAM_ORD) |>
  mutate(param_lab = factor(.PARAM_LAB[parameter],
                            levels = rev(.PARAM_LAB[.PARAM_ORD])),
         value = pmin(MARE_pct, 70))

p_err <- ggplot(err, aes(value, param_lab, colour = cond)) +
  geom_linerange(aes(xmin = 0, xmax = value, group = cond),
                 position = position_dodge(width = 0.6),
                 linewidth = 0.8, alpha = 0.4) +
  geom_point(position = position_dodge(width = 0.6), size = 2.6) +
  facet_wrap(~ struct_f, nrow = 1) +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "MARE (%)", y = "Parameter",
       title = "Estimation accuracy, scenario 16  (N = 80; lower = better)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUT, "fig_raw_error_N80.png"), p_err, width = 10, height = 6.5, dpi = 200)
ggsave(file.path(OUT, "fig_raw_error_N80.pdf"), p_err, width = 10, height = 6.5)

## ---- 4. selection-pattern heatmaps (2x2: structure x condition) ------------
build_covsel <- function(dir, cond) {
  read_csv(file.path(dir, "scm_covsel_by_covar.csv"), show_col_types = FALSE) |>
    filter(estimator == EST, outer_opt == OUT_OPT, sample_N == SN,
           structure %in% c("linCmt", "ode")) |>
    .fix_covsel_denom() |>
    mutate(cond = cond)
}
cov_raw <- bind_rows(build_covsel(ON_DIR, "warm-up ON"),
                     build_covsel(OFF_DIR, "warm-up OFF")) |>
  mutate(
    err_rate   = ifelse(is_true, n_FN / n_datasets, n_FP / n_datasets),
    err_signed = ifelse(is_true, -err_rate, err_rate) * 100,
    err_count  = ifelse(is_true, n_FN, n_FP),
    param      = coalesce(.VAR_LAB[tolower(var)], toupper(var)),
    effect     = paste0(param, "~", toupper(covar), ".",
                        .SHAPE_LAB[tolower(shape)] %||% shape),
    v_ord = match(tolower(var), names(.VAR_LAB)),
    c_ord = match(toupper(covar), toupper(.COVAR_ORD)),
    s_ord = match(tolower(shape), .SHAPE_ORD),
    scenario = factor(scenario, levels = sort(unique(scenario))),
    cond     = factor(cond, levels = names(PAL)),
    struct_f = factor(STRUCT_LAB[structure], levels = unname(STRUCT_LAB))) |>
  filter(!is.na(v_ord), !is.na(c_ord), !is.na(s_ord))

eff_lvls <- cov_raw |> distinct(effect, v_ord, c_ord, s_ord) |>
  arrange(v_ord, c_ord, s_ord) |> pull(effect)
cov_raw <- mutate(cov_raw, effect = factor(effect, levels = rev(eff_lvls)))

p_cov <- ggplot(cov_raw, aes(scenario, effect, fill = err_signed)) +
  geom_tile(colour = "grey85", linewidth = 0.3) +
  geom_tile(data = filter(cov_raw, is_true), fill = NA, colour = "black",
            linewidth = 0.2) +
  geom_text(aes(label = err_count), size = 2.7, fontface = "bold",
            colour = ifelse(cov_raw$err_rate > 0.55, "white", "grey15")) +
  facet_grid(cond ~ struct_f) +
  scale_fill_gradient2(low = "#08519C", mid = "white", high = "#B30000",
    midpoint = 0, limits = c(-100, 100), breaks = seq(-100, 100, 50),
    labels = c("FN 100%","FN 50%","0","FP 50%","FP 100%"), name = NULL) +
  labs(x = "Simulation scenario",
       y = "Covariate effect (param ~ covariate . shape)",
       title = "Covariate-selection error pattern  (N = 80)") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(), legend.position = "top",
        plot.title = element_text(face = "bold"),
        strip.text = element_text(size = 13)) +
  guides(fill = guide_colourbar(barwidth = 16))
ggsave(file.path(OUT, "fig_raw_covsel_N80.png"), p_cov, width = 13, height = 11, dpi = 200)
ggsave(file.path(OUT, "fig_raw_covsel_N80.pdf"), p_cov, width = 13, height = 11)

message("\nRaw warm-up on/off (linCmt + ODE) figures written under ", OUT)
