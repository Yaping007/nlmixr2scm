# ============================================================================
# make_vae_colinear_power_compare.R
# One-panel power comparison of the three PR #921 correlation-aware VAE
# covariate-selection arms (all N80, ODE structure):
#   colON080  = covSelectColinear = TRUE,  cut = 0.80
#   colON095  = covSelectColinear = TRUE,  cut = 0.95
#   colOoff95 = covSelectColinear = FALSE, cut = 0.95  (pre-PR baseline)
# Power = exact-match / N across the 16 scenarios. Dashed line = 80%.
# Usage:  source("script/viz/make_vae_colinear_power_compare.R")
# ============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(purrr)
})

FIG_DIR <- "output/figures/vae_colinear_compare"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

arms <- c(
  "Colinear ON, cut 0.80"         = "output/vae_covsel_colON080_N80_aggregated",
  "Colinear ON, cut 0.95"         = "output/vae_covsel_colON095_N80_aggregated",
  "Colinear OFF (baseline), 0.95" = "output/vae_covsel_colOoff95_N80_aggregated"
)

dat <- purrr::imap_dfr(arms, function(dir, lab) {
  csv <- file.path(dir, "vae_power.csv")
  if (!file.exists(csv)) stop("power CSV not found: ", csv)
  readr::read_csv(csv, show_col_types = FALSE) |>
    dplyr::transmute(scenario, structure, Power, arm = lab)
})

plot_df <- dat |>
  dplyr::mutate(
    scenario  = factor(scenario, levels = sort(unique(scenario))),
    arm       = factor(arm, levels = names(arms)),
    power_pct = Power * 100
  )

PAL_ARM <- c(
  "Colinear ON, cut 0.80"         = "#1F77B4",
  "Colinear ON, cut 0.95"         = "#E8820C",
  "Colinear OFF (baseline), 0.95" = "#7F7F7F"
)

p <- ggplot2::ggplot(plot_df,
    ggplot2::aes(scenario, power_pct, colour = arm, group = arm)) +
  ggplot2::geom_hline(yintercept = 80, linetype = "dashed",
                      colour = "grey40", linewidth = 0.4) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_colour_manual(values = PAL_ARM, name = "Correlation-aware arm") +
  ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                              labels = function(x) paste0(x, "%")) +
  ggplot2::labs(
    title    = "PR #921 correlation-aware VAE covariate selection (N80, ODE)",
    subtitle = "Power = exact recovery / N across the 16 scenarios",
    x = "Simulation scenario", y = "Power"
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "top",
    strip.background = ggplot2::element_rect(fill = "grey92", colour = NA)
  )

stub <- file.path(FIG_DIR, "fig_power_colinear_compare")
ggplot2::ggsave(paste0(stub, ".png"), p, width = 10, height = 5, dpi = 150)
ggplot2::ggsave(paste0(stub, ".pdf"), p, width = 10, height = 5)
message("saved: ", stub, ".{png,pdf}")

p
