# ==============================================================================
# aggregate_scm_timing.R
# ------------------------------------------------------------------------------
# Extract WALL-clock vs CPU timing from the schema-2.1 SCM bench records and
# contrast them to quantify the nlmixr2 fork-parallelism benefit.
#
# Directory scanned:
#   output/scm_bench/N<N>/scn<SS>_<structure>/<est>_<opt>/res_ds<DDD>.rds
# (the `.fit.rds` sidecars are skipped -- they are raw fits, not records.)
#
# TIMING FIELDS PER RECORD (written by package_scm_schema21()):
#   r$runtime$base_sec / scm_sec / refit_sec / total_sec   -- WALL seconds/phase
#   r$cpu$base_sec / scm_sec / refit_sec / total_sec        -- CPU seconds/phase
#   r$cpu$hog_factor = cpu$total_sec / runtime$total_sec    -- realised speedup
#   r$scm_workers, r$rx_threads                             -- parallelism knobs
#   flat mirrors: r$wall_total_sec, r$cpu_total_sec, r$hog_factor
#
# INTERPRETING hog_factor (READ THIS):
#   hog = CPU_total / WALL_total is "how many cores of work per wall-second".
#     hog > 1  -> parallel speedup captured (fork workers' CPU was reaped).
#     hog ~ 1  -> effectively serial.
#     hog < 1  -> forked-worker CPU NOT captured by proc.time()'s child columns
#                 (a MEASUREMENT artefact, NOT a slowdown). In that regime the
#                 in-record CPU under-reports total work and the TRUE benefit
#                 must be read from LSF `bacct -l <jobid>` CPU_T instead.
#   The summary below prints the hog distribution so you can see which regime
#   your data is in before drawing conclusions.
#
# OUTPUTS (to output/scm_timing/):
#   scm_timing_long.csv       one row per fit: phase wall/cpu + hog + knobs
#   scm_timing_by_cell.csv    per (N, scenario, structure, est, opt): medians
#   scm_timing_by_group.csv   per (N, structure): medians + hog summary
#   scm_wall_vs_cpu.png       scatter: wall vs cpu, y=x reference (benefit line)
#   scm_hog_hist.png          hog_factor distribution, faceted by structure
#   scm_wall_by_N.png         wall-time growth with cohort N, by structure
#
# Usage:
#   Rscript script/aggregate_scm_timing.R [--root output/scm_bench] \
#                                         [--out_dir output/scm_timing]
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

# ---- CLI -------------------------------------------------------------------
parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  get <- function(flag, default) {
    i <- match(flag, args)
    if (is.na(i) || i == length(args)) default else args[[i + 1L]]
  }
  list(
    root    = get("--root",    "output/scm_bench"),
    out_dir = get("--out_dir", "output/scm_timing")
  )
}

# ---- record -> one tidy timing row -----------------------------------------
# Path is the backstop for coordinate keys the record may be missing.
read_timing <- function(path) {
  r <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(r)) {
    return(tibble(file = path, status = "unreadable"))
  }

  # path coordinates: .../N<N>/scn<SS>_<structure>/<est>_<opt>/res_ds<DDD>.rds
  m_N     <- str_match(path, "/N(\\d+)/")[, 2]
  m_scn   <- str_match(path, "/scn(\\d+)_([A-Za-z]+)/")[, 2:3]
  m_eo    <- str_match(path, "/([a-z0-9]+)_([a-z0-9]+)/res_ds")[, 2:3]
  m_ds    <- str_match(path, "res_ds(\\d+)\\.rds$")[, 2]

  wall <- r$runtime %||% list()
  cpu  <- r$cpu     %||% list()

  tibble(
    file        = path,
    sample_N    = as.integer(r$sample_N    %||% m_N),
    scenario_id = as.integer(r$scenario_id %||% m_scn[[1]]),
    structure   = r$structure %||% m_scn[[2]] %||% r$model_type %||% NA_character_,
    estimator   = r$estimator %||% m_eo[[1]] %||% NA_character_,
    outer_opt   = r$outer_opt %||% m_eo[[2]] %||% NA_character_,
    dataset_id  = as.integer(r$dataset_id  %||% m_ds),
    status      = r$status %||% NA_character_,
    scm_workers = as.integer(r$scm_workers %||% NA),
    rx_threads  = as.integer(r$rx_threads  %||% NA),
    # WALL (elapsed) seconds
    wall_base   = as.numeric(wall$base_sec  %||% NA),
    wall_scm    = as.numeric(wall$scm_sec   %||% NA),
    wall_refit  = as.numeric(wall$refit_sec %||% NA),
    wall_total  = as.numeric(wall$total_sec %||% r$wall_total_sec %||% NA),
    # CPU (self + reaped child) seconds
    cpu_base    = as.numeric(cpu$base_sec   %||% NA),
    cpu_scm     = as.numeric(cpu$scm_sec    %||% NA),
    cpu_refit   = as.numeric(cpu$refit_sec  %||% NA),
    cpu_total   = as.numeric(cpu$total_sec  %||% r$cpu_total_sec %||% NA),
    hog_factor  = as.numeric(cpu$hog_factor %||% r$hog_factor %||% NA)
  )
}

# ---- main ------------------------------------------------------------------
main <- function() {
  args <- parse_args()
  dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(
    args$root,
    pattern    = "^res_ds\\d+\\.rds$",   # excludes *.fit.rds and *.meta.json
    recursive  = TRUE,
    full.names = TRUE
  )
  message(sprintf("Found %d record files under %s", length(files), args$root))
  if (length(files) == 0L) stop("No res_ds*.rds records found. Check --root.")

  long <- map_dfr(files, read_timing)

  ok <- long %>% filter(!is.na(wall_total))
  message(sprintf("Parsed %d records (%d with usable wall time).",
                  nrow(long), nrow(ok)))

  # ---- hog_factor regime report (the key diagnostic) ----------------------
  hog_ok <- ok %>% filter(is.finite(hog_factor))
  if (nrow(hog_ok) > 0) {
    qs <- quantile(hog_ok$hog_factor, c(.05, .25, .5, .75, .95), na.rm = TRUE)
    frac_gt1 <- mean(hog_ok$hog_factor > 1, na.rm = TRUE)
    message("\n--- hog_factor (CPU_total / WALL_total) ---")
    message(sprintf("  median = %.2f   IQR = [%.2f, %.2f]   p05/p95 = %.2f / %.2f",
                    qs[3], qs[2], qs[4], qs[1], qs[5]))
    message(sprintf("  fraction with hog > 1 (parallel benefit visible): %.1f%%",
                    100 * frac_gt1))
    if (qs[3] < 1) {
      message("  NOTE: median hog < 1 -> forked-worker CPU is NOT captured in")
      message("        the records (proc.time child cols empty). The in-record")
      message("        CPU under-reports work; cross-check true CPU via")
      message("        `bacct -l <jobid>` CPU_T for the parallelism story.")
    }
  }

  # ---- per-cell medians ---------------------------------------------------
  by_cell <- ok %>%
    group_by(sample_N, scenario_id, structure, estimator, outer_opt) %>%
    summarise(
      n          = dplyr::n(),
      wall_med   = median(wall_total, na.rm = TRUE),
      wall_iqr   = IQR(wall_total, na.rm = TRUE),
      cpu_med    = median(cpu_total,  na.rm = TRUE),
      hog_med    = median(hog_factor, na.rm = TRUE),
      wall_scm_med   = median(wall_scm,   na.rm = TRUE),
      cpu_scm_med    = median(cpu_scm,    na.rm = TRUE),
      wall_refit_med = median(wall_refit, na.rm = TRUE),
      workers    = dplyr::first(scm_workers),
      rx_threads = dplyr::first(rx_threads),
      .groups    = "drop"
    ) %>%
    arrange(structure, sample_N, scenario_id, estimator, outer_opt)

  # ---- per (N, structure) rollup ------------------------------------------
  by_group <- ok %>%
    group_by(sample_N, structure) %>%
    summarise(
      n         = dplyr::n(),
      wall_med  = median(wall_total, na.rm = TRUE),
      cpu_med   = median(cpu_total,  na.rm = TRUE),
      hog_med   = median(hog_factor, na.rm = TRUE),
      hog_p95   = quantile(hog_factor, .95, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    arrange(structure, sample_N)

  # ---- write CSVs ---------------------------------------------------------
  write_csv(long,     file.path(args$out_dir, "scm_timing_long.csv"))
  write_csv(by_cell,  file.path(args$out_dir, "scm_timing_by_cell.csv"))
  write_csv(by_group, file.path(args$out_dir, "scm_timing_by_group.csv"))

  # ---- plots --------------------------------------------------------------
  # 1. wall vs cpu: the y=x line is the "no benefit" boundary. Points ABOVE
  #    (cpu > wall) mean parallel work compressed into less wall time.
  p_wc <- ggplot(ok, aes(wall_total, cpu_total, colour = structure)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
    geom_point(alpha = 0.35, size = 1) +
    coord_equal() +
    labs(
      title    = "SCM fit: CPU work vs wall-clock time",
      subtitle = "Points above the dashed y = x line = parallel speedup (CPU > wall)",
      x = "Wall-clock time (s)", y = "CPU time (s)", colour = "Structure"
    ) +
    theme_minimal()
  ggsave(file.path(args$out_dir, "scm_wall_vs_cpu.png"), p_wc,
         width = 7, height = 6, dpi = 150)

  # 2. hog_factor distribution by structure
  p_hog <- ok %>%
    filter(is.finite(hog_factor)) %>%
    ggplot(aes(hog_factor, fill = structure)) +
    geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey30") +
    facet_wrap(~ structure, ncol = 1, scales = "free_y") +
    labs(
      title    = "Realised parallelism: hog_factor = CPU / wall",
      subtitle = "Right of the dashed line (>1) = speedup from fork workers",
      x = "hog_factor", y = "count"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  ggsave(file.path(args$out_dir, "scm_hog_hist.png"), p_hog,
         width = 7, height = 6, dpi = 150)

  # 3. wall-time growth with cohort N, per structure (median +/- IQR band)
  p_N <- by_group %>%
    ggplot(aes(factor(sample_N), wall_med, colour = structure,
               group = structure)) +
    geom_line() + geom_point(size = 2) +
    labs(
      title = "Median SCM wall-time by cohort size",
      x = "Cohort N", y = "Median wall time (s)", colour = "Structure"
    ) +
    theme_minimal()
  ggsave(file.path(args$out_dir, "scm_wall_by_N.png"), p_N,
         width = 7, height = 5, dpi = 150)

  message(sprintf("\nWrote CSVs + 3 plots to %s", args$out_dir))
  by_group
}

if (sys.nframe() == 0L) {
  print(main())
}
