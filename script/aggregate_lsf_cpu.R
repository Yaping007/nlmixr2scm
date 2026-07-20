# ==============================================================================
# aggregate_lsf_cpu.R
# ------------------------------------------------------------------------------
# HONEST parallelism story for the SCM bench, using LSF-accounted CPU time.
#
# WHY THIS EXISTS
# ---------------
# The in-record cpu_total_sec (from R's proc.time()) does NOT see the CPU burned
# by runSCM()'s forked workers -- proc.time()'s user.child/sys.child are only
# filled for children reaped by wait(), which the fork pool does not do. Hence
# in-record hog_factor came out ~0.03-0.25 (a MEASUREMENT artefact, not a real
# slowdown).
#
# LSF, by contrast, accounts CPU time for the WHOLE process tree (parent + all
# forked workers + OpenMP threads). That number is written by LSF into two
# equivalent places:
#   (a) the "Resource usage summary" block appended to each task's -o .out log
#   (b) `bacct -l <jobid>` output on the cluster
# Both report the same `CPU time` and `Run time`. This script parses (a) by
# default (already on disk, self-keying via job name) and (b) as a fallback.
#
#   true_hog = LSF CPU_time / LSF Run_time   (>1 == real fork speedup)
#
# JOB -> COORDINATE MAPPING (from submit_one_array.sh + bench_array.lsf)
#   job name : bench_N<N>_scn<SS>_<structure>_<est>_<opt>
#   log file : logs/bench_N<N>/<jobname>.<jobid>.<idx>.out   (idx == dataset_id)
#
# HOW TO GET THE bacct FALLBACK (run on the cluster if .out logs were rotated):
#   for jid in <jobid1> <jobid2> ...; do bacct -l "$jid"; done > logs/bacct_dump.txt
#   # or capture the array jobs you still have:
#   bacct -l -J 'bench_*' > logs/bacct_dump.txt
#
# OUTPUTS (to output/scm_timing/):
#   lsf_cpu_long.csv        one row per task: cpu_t, run_t, true_hog + coords
#   lsf_cpu_by_group.csv    per (N, structure): median cpu_t / run_t / true_hog
#   lsf_true_hog_hist.png   true_hog distribution, faceted by structure
#   lsf_cpu_vs_wall.png     LSF CPU vs LSF wall, y=x reference
#
# Usage:
#   Rscript script/aggregate_lsf_cpu.R [--logs logs] [--bacct logs/bacct_dump.txt] \
#                                      [--out_dir output/scm_timing]
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

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  get <- function(flag, default) {
    i <- match(flag, args)
    if (is.na(i) || i == length(args)) default else args[[i + 1L]]
  }
  list(
    logs    = get("--logs",    "logs"),
    bacct   = get("--bacct",   NA_character_),   # optional fallback dump
    out_dir = get("--out_dir", "output/scm_timing")
  )
}

# Job name -> coordinates. structure is fixed vocab so est/opt split cleanly.
parse_jobname <- function(x) {
  m <- str_match(
    x, "bench_N(\\d+)_scn(\\d+)_(linCmt|ode)_([A-Za-z]+)_([A-Za-z0-9]+)"
  )
  tibble(
    sample_N    = as.integer(m[, 2]),
    scenario_id = as.integer(m[, 3]),
    structure   = m[, 4],
    estimator   = m[, 5],
    outer_opt   = m[, 6]
  )
}

# --- source (a): LSF -o .out logs -------------------------------------------
# Each log's tail has:
#   Resource usage summary:
#       CPU time :   1234.56 sec.
#       ...
#       Run time :   456 sec.
read_out_log <- function(path) {
  txt <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
  if (length(txt) == 0L) return(NULL)
  blob <- paste(txt, collapse = "\n")

  # LSF may append MULTIPLE "Resource usage summary" blocks (e.g. on requeue);
  # take the LAST CPU time / Run time, which reflects the successful run.
  cpu_all <- str_match_all(blob, "CPU time\\s*:\\s*([0-9.]+)")[[1]][, 2]
  run_all <- str_match_all(blob, "Run time\\s*:\\s*([0-9.]+)")[[1]][, 2]
  cpu_t <- if (length(cpu_all)) as.numeric(tail(cpu_all, 1)) else NA_real_
  run_t <- if (length(run_all)) as.numeric(tail(run_all, 1)) else NA_real_
  if (is.na(cpu_t) && is.na(run_t)) return(NULL)   # job hadn't finished / no summary

  # Coordinates come from the FILENAME, which is unambiguous:
  #   bench_N<N>_scn<SS>_<structure>_<est>_<opt>.<jobid>.<idx>.out
  # NB: do NOT parse the echoed "Job name :" line -- the .out embeds the LSBATCH
  # script text, whose literal `echo "Job name : ${LSB_JOBNAME}"` matches first.
  fn <- basename(path)
  jobid <- str_match(fn, "\\.(\\d+)\\.\\d+\\.out$")[, 2]
  idx   <- str_match(fn, "\\.(\\d+)\\.out$")[, 2]

  coords <- parse_jobname(fn)
  bind_cols(
    tibble(
      file       = path,
      source     = "out_log",
      jobid      = as.integer(jobid),
      dataset_id = as.integer(idx),
      cpu_t_sec  = cpu_t,
      run_t_sec  = run_t
    ),
    coords
  )
}

# --- source (b): bacct -l dump ----------------------------------------------
# Splits on "-----" record separators; each record has a "Job <...>, Job Name
# <name>" header and "CPU_TIME"/"TURNAROUND"/"RUN_TIME" fields.
read_bacct_dump <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  raw <- readLines(path, warn = FALSE)
  blob <- paste(raw, collapse = "\n")
  # bacct separates accounting records with rows of dashes
  recs <- str_split(blob, "\\n-{5,}\\n")[[1]]
  map_dfr(recs, function(r) {
    if (!str_detect(r, "Job Name")) return(NULL)
    jobname <- str_match(r, "Job Name <([^>]+)>")[, 2]
    jobidx  <- str_match(r, "Job <\\d+\\[(\\d+)\\]>")[, 2]     # array element
    jobid   <- str_match(r, "Job <(\\d+)")[, 2]
    cpu_t   <- as.numeric(str_match(r, "CPU_T(?:IME)?\\s*[: ]\\s*([0-9.]+)")[, 2])
    run_t   <- as.numeric(str_match(r, "RUN_T(?:IME)?\\s*[: ]\\s*([0-9.]+)")[, 2])
    if (is.na(cpu_t)) return(NULL)
    bind_cols(
      tibble(
        file       = path,
        source     = "bacct",
        jobid      = as.integer(jobid),
        dataset_id = as.integer(jobidx),
        cpu_t_sec  = cpu_t,
        run_t_sec  = run_t
      ),
      parse_jobname(jobname)
    )
  })
}

main <- function() {
  args <- parse_args()
  dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

  out_files <- list.files(args$logs, pattern = "\\.out$",
                          recursive = TRUE, full.names = TRUE)
  out_files <- out_files[str_detect(basename(out_files), "^bench_N")]
  message(sprintf("Found %d bench .out logs under %s", length(out_files), args$logs))

  from_logs  <- if (length(out_files)) compact(map(out_files, read_out_log)) else list()
  from_logs  <- if (length(from_logs)) bind_rows(from_logs) else NULL
  from_bacct <- read_bacct_dump(args$bacct)

  dat <- bind_rows(from_logs, from_bacct)
  if (is.null(dat) || nrow(dat) == 0L) {
    stop("No LSF CPU records parsed. Point --logs at the logs/ dir on the ",
         "cluster, or supply --bacct <dump.txt>.")
  }

  # dedupe: prefer out_log over bacct for the same (jobid, dataset_id)
  dat <- dat %>%
    filter(!is.na(cpu_t_sec), !is.na(run_t_sec), run_t_sec > 0) %>%
    mutate(true_hog = cpu_t_sec / run_t_sec) %>%
    arrange(desc(source == "out_log")) %>%
    distinct(sample_N, scenario_id, structure, estimator, outer_opt,
             dataset_id, .keep_all = TRUE)

  message(sprintf("Usable LSF task records: %d", nrow(dat)))

  # guard: catch a coordinate-parse regression before it silently collapses
  n_bad <- sum(is.na(dat$structure) | is.na(dat$sample_N))
  if (n_bad > 0) {
    warning(sprintf("%d records have NA coordinates (jobname parse failed).",
                    n_bad))
  }

  # --- headline diagnostic --------------------------------------------------
  qs <- quantile(dat$true_hog, c(.05, .25, .5, .75, .95), na.rm = TRUE)
  message("\n--- TRUE hog (LSF CPU_time / LSF Run_time) ---")
  message(sprintf("  median = %.2f   IQR = [%.2f, %.2f]   p05/p95 = %.2f / %.2f",
                  qs[3], qs[2], qs[4], qs[1], qs[5]))
  message(sprintf("  fraction with true_hog > 1 (real speedup): %.1f%%",
                  100 * mean(dat$true_hog > 1, na.rm = TRUE)))

  by_group <- dat %>%
    group_by(sample_N, structure) %>%
    summarise(
      n            = dplyr::n(),
      cpu_t_med    = median(cpu_t_sec, na.rm = TRUE),
      run_t_med    = median(run_t_sec, na.rm = TRUE),
      true_hog_med = median(true_hog,  na.rm = TRUE),
      true_hog_p95 = quantile(true_hog, .95, na.rm = TRUE),
      .groups      = "drop"
    ) %>%
    arrange(structure, sample_N)

  write_csv(dat,      file.path(args$out_dir, "lsf_cpu_long.csv"))
  write_csv(by_group, file.path(args$out_dir, "lsf_cpu_by_group.csv"))

  # --- plots ----------------------------------------------------------------
  p_hog <- ggplot(dat, aes(true_hog, fill = structure)) +
    geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey30") +
    facet_wrap(~ structure, ncol = 1, scales = "free_y") +
    labs(
      title    = "Honest parallelism: LSF CPU_time / Run_time",
      subtitle = "Right of the dashed line (>1) = fork workers really sped it up",
      x = "true_hog", y = "count"
    ) +
    theme_minimal() + theme(legend.position = "none")
  ggsave(file.path(args$out_dir, "lsf_true_hog_hist.png"), p_hog,
         width = 7, height = 6, dpi = 150)

  p_cw <- ggplot(dat, aes(run_t_sec, cpu_t_sec, colour = structure)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
    geom_point(alpha = 0.35, size = 1) +
    labs(
      title    = "LSF CPU vs wall (Run time)",
      subtitle = "Points above y = x = CPU work compressed into less wall time",
      x = "LSF Run time (s)", y = "LSF CPU time (s)", colour = "Structure"
    ) +
    theme_minimal()
  ggsave(file.path(args$out_dir, "lsf_cpu_vs_wall.png"), p_cw,
         width = 7, height = 6, dpi = 150)

  message(sprintf("\nWrote CSVs + 2 plots to %s", args$out_dir))
  by_group
}

if (sys.nframe() == 0L) {
  print(main())
}
