# ============================================================================
# make_lbfgsb3c_raw_compare.R
# RAW warm-start ON vs OFF comparison for the focei + lbfgsb3c benchmark,
# N = 80, ODE only (the lbfgsb3c warm runs contain ode rows only). Mirrors
# make_warmup_raw_compare.R (bobyqa) with a "_lbfgsb3c" filename suffix.
#
#   warm-up ON  : output/lbfgsb3c_est703_08132026_warmon_aggregated
#   warm-up OFF : output/lbfgsb3c_est703_08132026_warmoff_aggregated
#   PsN-SCM     : output/psn_scm_combined_aggregated   (power overlay only)
#
# Figures (output/figures/warmup_compare/raw/):
#   fig_raw_convergence_N80_lbfgsb3c.{png,pdf}
#   fig_raw_power_N80_lbfgsb3c.{png,pdf}
#   fig_raw_error_N80_lbfgsb3c.{png,pdf}
#   fig_raw_covsel_N80_lbfgsb3c.{png,pdf}
#
# Usage:  source("script/viz/make_lbfgsb3c_raw_compare.R")
# ============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

## helpers for the covsel heatmap (.fix_covsel_denom, .VAR_LAB, .COVAR_ORD,
## .SHAPE_ORD, .SHAPE_LAB, %||%)
source("script/viz/fig_covsel_heatmap.R")

ON_DIR  <- "output/lbfgsb3c_est703_08132026_warmon_aggregated"
OFF_DIR <- "output/lbfgsb3c_est703_08132026_warmoff_aggregated"
PSN_DIR <- "output/psn_scm_combined_aggregated"
OUT     <- "output/figures/warmup_compare/raw"
EST <- "focei"; OUT_OPT <- "lbfgsb3c"; SN <- 80
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

PAL <- c("warm-up ON" = "#0F8B8D", "warm-up OFF" = "#C1666B",
         "PsN-SCM" = "#2B2B2B")
STRUCT_LAB <- c(ode = "ODE")
STRUCT_MAP <- c(ode = "ODE")   # PsN "ode" folds to the same ODE panel

## load one metric CSV from both warm folders, tag condition (N=80, focei/lbfgsb3c, ode)
load_both <- function(csv) {
  rd <- function(dir, cond) read_csv(file.path(dir, csv), show_col_types = FALSE) |>
    filter(estimator == EST, outer_opt == OUT_OPT, sample_N == SN,
           structure == "ode") |>
    mutate(cond = cond)
  bind_rows(rd(ON_DIR, "warm-up ON"), rd(OFF_DIR, "warm-up OFF")) |>
    mutate(cond     = factor(cond, levels = names(PAL)),
           struct_f = factor(STRUCT_LAB[structure], levels = unname(STRUCT_LAB)))
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
       title = "Convergence & covariance step  (N = 80) [lbfgsb3c]") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUT, "fig_raw_convergence_N80_lbfgsb3c.png"), p_conv, width = 7, height = 6.5, dpi = 200)
ggsave(file.path(OUT, "fig_raw_convergence_N80_lbfgsb3c.pdf"), p_conv, width = 7, height = 6.5)

## ---- 2. power (warm-up ON/OFF + PsN-SCM) -----------------------------------
read_power <- function(dir, cond, est, oopt) {
  read_csv(file.path(dir, "scm_power.csv"), show_col_types = FALSE) |>
    filter(estimator == est, outer_opt == oopt, sample_N == SN,
           structure == "ode") |>
    transmute(scenario, Power = Power * 100,
              struct_f = factor("ODE", levels = unname(STRUCT_LAB)),
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
       title = "Selection power  (N = 80) [lbfgsb3c]") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUT, "fig_raw_power_N80_lbfgsb3c.png"), p_pow, width = 7, height = 5.5, dpi = 200)
ggsave(file.path(OUT, "fig_raw_power_N80_lbfgsb3c.pdf"), p_pow, width = 7, height = 5.5)

## ---- 3. estimation accuracy (MARE%, scenario 16) ---------------------------
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
       title = "Estimation accuracy, scenario 16  (N = 80; lower = better) [lbfgsb3c]") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUT, "fig_raw_error_N80_lbfgsb3c.png"), p_err, width = 7, height = 6.5, dpi = 200)
ggsave(file.path(OUT, "fig_raw_error_N80_lbfgsb3c.pdf"), p_err, width = 7, height = 6.5)

## ---- 4. selection-pattern heatmaps (rows: condition) -----------------------
## NB: some lbfgsb3c aggregated dirs ship an EMPTY scm_covsel_by_covar.csv
## (e.g. warm-off); build_covsel skips any dir with no usable rows.
build_covsel <- function(dir, cond) {
  d <- read_csv(file.path(dir, "scm_covsel_by_covar.csv"), show_col_types = FALSE)
  if (nrow(d) == 0 || !"estimator" %in% names(d)) {
    warning("No covsel rows for '", cond, "' (", basename(dir), ") - skipped")
    return(NULL)
  }
  d |>
    filter(estimator == EST, outer_opt == OUT_OPT, sample_N == SN,
           structure == "ode") |>
    .fix_covsel_denom() |>
    mutate(cond = cond)
}
cov_raw <- bind_rows(build_covsel(ON_DIR, "warm-up ON"),
                     build_covsel(OFF_DIR, "warm-up OFF"))
if (nrow(cov_raw) == 0) {
  warning("No covariate-selection data for lbfgsb3c - skipping covsel heatmap.")
} else {
cov_raw <- cov_raw |>
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
       title = "Covariate-selection error pattern  (N = 80) [lbfgsb3c]") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(), legend.position = "top",
        plot.title = element_text(face = "bold"),
        strip.text = element_text(size = 13)) +
  guides(fill = guide_colourbar(barwidth = 16))
ggsave(file.path(OUT, "fig_raw_covsel_N80_lbfgsb3c.png"), p_cov, width = 8, height = 11, dpi = 200)
ggsave(file.path(OUT, "fig_raw_covsel_N80_lbfgsb3c.pdf"), p_cov, width = 8, height = 11)
}

message("\nRaw warm-up on/off (lbfgsb3c, ODE) figures written under ", OUT)
