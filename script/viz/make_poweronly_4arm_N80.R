# ============================================================================
# make_poweronly_4arm_N80.R
# ----------------------------------------------------------------------------
# Selection-power overlay, TWO panels (Analytic | ODE), four arms
# (power-shape-only space, N=80):
#   nlmixr2-SCM bobyqa   (teal)    output/scm_poweronly_est703_08132026_aggregated
#   nlmixr2-SCM lbfgsb3c (orange)  output/lbfgsb3c_poweronly_N80_08162026_warmon_aggregated
#   PsN-SCM  (nonmem)    (black)   advan4 (Analytic) + ode dirs
#   nlmixr2-VAE          (purple)  output/vae_covsel_poweronly_N80_08162026_aggregated
#
# Power = exact recovery (Power = n_exact / N), scaled to %. Structure folded:
#   linCmt / advan4 -> "Analytic"  ;  ode -> "ODE".
#
# Output: output/figures/poweronly_N80/fig_poweronly_4arm_N80.{png,pdf}
# Usage:  source("script/viz/make_poweronly_4arm_N80.R")
# ============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
})

SN     <- 80
METRIC <- "Power"          # exact-recovery power
OUT    <- "output/figures/poweronly_N80"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

PAL <- c("nlmixr2-SCM (bobyqa)"   = "#0F8B8D",
         "nlmixr2-SCM (lbfgsb3c)" = "#E8820C",
         "PsN-SCM"                = "#2B2B2B",
         "nlmixr2-VAE"            = "#A64AC9")

## Read a power CSV, keep N=80 rows (+ est/opt filter), fold structure into
## Analytic/ODE, return scenario + Power(%) tagged with the arm label.
read_power <- function(dir, csv, arm, est = NULL, oopt = NULL) {
  d <- read_csv(file.path(dir, csv), show_col_types = FALSE) |>
    filter(sample_N == SN)
  if (!is.null(est)  && "estimator" %in% names(d)) d <- filter(d, estimator == est)
  if (!is.null(oopt) && "outer_opt" %in% names(d)) d <- filter(d, outer_opt == oopt)
  d |>
    transmute(scenario, Power = .data[[METRIC]] * 100, arm = arm,
              struct = case_when(structure %in% c("linCmt", "advan4") ~ "Analytic",
                                 structure == "ode" ~ "ODE",
                                 TRUE ~ structure))
}

pw <- bind_rows(
  read_power("output/scm_poweronly_est703_08132026_aggregated",
             "scm_power.csv", "nlmixr2-SCM (bobyqa)",   "focei", "bobyqa"),
  read_power("output/lbfgsb3c_poweronly_N80_08162026_warmon_aggregated",
             "scm_power.csv", "nlmixr2-SCM (lbfgsb3c)", "focei", "lbfgsb3c"),
  # PsN power-only ships split by structure: advan4 -> Analytic, ode -> ODE
  read_power("output/psn_scm_poweronly_advan4_aggregated",
             "scm_power.csv", "PsN-SCM",                "nonmem_scm", "focei"),
  read_power("output/psn_scm_poweronly_ode_aggregated",
             "scm_power.csv", "PsN-SCM",                "nonmem_scm", "focei"),
  read_power("output/vae_covsel_poweronly_N80_08162026_aggregated",
             "vae_power.csv", "nlmixr2-VAE")
) |>
  mutate(arm    = factor(arm, levels = names(PAL)),
         struct = factor(struct, levels = c("Analytic", "ODE")))

missing <- setdiff(names(PAL), as.character(unique(pw$arm)))
if (length(missing))
  warning("No N=80 power rows for: ", paste(missing, collapse = ", "))

p_pow <- ggplot(pw, aes(scenario, Power, colour = arm, group = arm)) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "grey70") +
  geom_line(linewidth = 0.9) + geom_point(size = 1.9) +
  facet_wrap(~ struct, nrow = 1) +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_x_continuous(breaks = seq(2, 16, by = 2)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Simulation scenario",
       y = "Power (true covariate + shape detected)",
       title = "Selection power \u2014 power-only covariate space  (N = 80)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold")) +
  guides(colour = guide_legend(nrow = 1))

ggsave(file.path(OUT, "fig_poweronly_4arm_N80.png"), p_pow,
       width = 12, height = 5.5, dpi = 200)
ggsave(file.path(OUT, "fig_poweronly_4arm_N80.pdf"), p_pow,
       width = 12, height = 5.5)

message("Wrote fig_poweronly_4arm_N80.{png,pdf} under ", OUT)
