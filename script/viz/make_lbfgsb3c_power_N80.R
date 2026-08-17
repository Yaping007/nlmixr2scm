# ============================================================================
# make_lbfgsb3c_power_N80.R
# Selection-power overlay for the nlmixr2-SCM lbfgsb3c benchmark (N = 80, ODE).
#
# Three platforms overlaid in ONE panel (Power = true covariate detected):
#   warm-up ON  (teal)   : output/lbfgsb3c_est703_08132026_warmon_aggregated
#   warm-up OFF (salmon) : output/lbfgsb3c_est703_08132026_warmoff_aggregated
#   PsN-SCM     (black)  : output/psn_scm_combined_aggregated
#
# Power comes precomputed from each dir's scm_power.csv (Power = n_exact / N),
# scaled to a percentage. Only the ODE structure is plotted (the lbfgsb3c warm
# runs contain ode rows only); PsN "ode" maps to the same panel.
#
# Output (output/figures/warmup_compare/raw/):
#   fig_lbfgsb3c_power_N80.{png,pdf}
#
# Usage:  source("script/viz/make_lbfgsb3c_power_N80.R")
# ============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
})

ON_DIR  <- "output/lbfgsb3c_est703_08132026_warmon_aggregated"
OFF_DIR <- "output/lbfgsb3c_est703_08132026_warmoff_aggregated"
PSN_DIR <- "output/psn_scm_combined_aggregated"
OUT     <- "output/figures/warmup_compare/raw"
SN      <- 80
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

PAL <- c("warm-up ON" = "#0F8B8D", "warm-up OFF" = "#C1666B",
         "PsN-SCM" = "#2B2B2B")

## Read scm_power.csv from one aggregated dir, keep N=80 ODE rows for the given
## estimator/optimizer, return scenario + Power(%) tagged with the platform.
read_power <- function(dir, cond, est, oopt) {
  read_csv(file.path(dir, "scm_power.csv"), show_col_types = FALSE) |>
    filter(estimator == est, outer_opt == oopt, sample_N == SN,
           structure == "ode") |>
    transmute(scenario, Power = Power * 100, cond = cond)
}

pw <- bind_rows(
  read_power(ON_DIR,  "warm-up ON",  "focei",       "lbfgsb3c"),
  read_power(OFF_DIR, "warm-up OFF", "focei",       "lbfgsb3c"),
  read_power(PSN_DIR, "PsN-SCM",     "nonmem_scm",  "focei")) |>
  mutate(cond = factor(cond, levels = names(PAL)))

## sanity check: warn if any platform failed to contribute rows
missing <- setdiff(names(PAL), as.character(unique(pw$cond)))
if (length(missing))
  warning("No N=80 ODE power rows for: ", paste(missing, collapse = ", "))

p_pow <- ggplot(pw, aes(scenario, Power, colour = cond, group = cond)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_x_continuous(breaks = 1:16) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Simulation scenario", y = "Power (true covariate detected)",
       title = "Selection power \u2014 lbfgsb3c  (N = 80, ODE)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT, "fig_lbfgsb3c_power_N80.png"), p_pow,
       width = 8, height = 5.5, dpi = 200)
ggsave(file.path(OUT, "fig_lbfgsb3c_power_N80.pdf"), p_pow,
       width = 8, height = 5.5)

message("Wrote fig_lbfgsb3c_power_N80.{png,pdf} under ", OUT)
