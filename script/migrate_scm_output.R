# ============================================================================
# migrate_scm_output.R
# ----------------------------------------------------------------------------
# One-time relocation of the legacy SCM output tree into the new bench layout.
#
#   OLD:  output/output_N{NN}/output_N{NN}_sc{SS}/res_ds{DDD}.rds
#   NEW:  output/scm_bench/N{NN}/scn{SS}/focei_bobyqa/res_ds{DDD}.rds
#
# Also backfills 5 coordinate keys inside each RDS so downstream aggregation
# never needs path parsing:
#   sample_N, scenario_id, dataset_id, estimator, outer_opt
#
# Safe to re-run: file.copy(..., overwrite = FALSE) skips existing targets.
#
# Usage (from repo root):
#   Rscript script/migrate_scm_output.R
#   # or:  source("script/migrate_scm_output.R"); migrate_scm_output()
#   # dry run first:
#   source("script/migrate_scm_output.R"); migrate_scm_output(dry_run = TRUE)
# ============================================================================

migrate_scm_output <- function(old_root  = "output",
                               new_root  = "output/scm_bench",
                               estimator = "focei",
                               outer_opt = "bobyqa",
                               dry_run   = FALSE,
                               verbose   = TRUE) {
  cell_tag <- paste(estimator, outer_opt, sep = "_")

  old_dirs <- list.files(old_root, pattern = "^output_N\\d+$",
                         full.names = TRUE)
  if (!length(old_dirs)) {
    message("No legacy output_N* directories found under ", old_root)
    return(invisible(NULL))
  }

  n_copied <- 0L; n_backfilled <- 0L; n_skipped <- 0L

  for (d in old_dirs) {
    N <- as.integer(sub("^output_N", "", basename(d)))
    sc_dirs <- list.files(d, pattern = "^output_N\\d+_sc\\d+$",
                          full.names = TRUE)

    for (sd in sc_dirs) {
      sc <- as.integer(sub("^.*_sc(\\d+)$", "\\1", basename(sd)))
      tgt <- file.path(new_root, sprintf("N%d", N),
                       sprintf("scn%02d", sc), cell_tag)
      if (!dry_run) dir.create(tgt, showWarnings = FALSE, recursive = TRUE)

      rds <- list.files(sd, "^res_ds\\d+\\.rds$", full.names = TRUE)
      for (f in rds) {
        ds_id <- as.integer(sub("^res_ds(\\d+)\\.rds$", "\\1", basename(f)))
        new_f <- file.path(tgt, basename(f))

        if (file.exists(new_f)) { n_skipped <- n_skipped + 1L; next }

        if (dry_run) {
          if (verbose) message(sprintf("[dry] %s -> %s", f, new_f))
          next
        }

        ok <- file.copy(f, new_f, overwrite = FALSE, copy.date = TRUE)
        if (!ok) { warning("copy failed: ", f); next }
        n_copied <- n_copied + 1L

        # Backfill coordinate keys inside the RDS
        r <- tryCatch(readRDS(new_f), error = function(e) NULL)
        if (is.null(r) || !is.list(r)) next
        added <- FALSE
        if (is.null(r$sample_N))    { r$sample_N    <- N;         added <- TRUE }
        if (is.null(r$scenario_id)) { r$scenario_id <- sc;        added <- TRUE }
        if (is.null(r$dataset_id))  { r$dataset_id  <- ds_id;     added <- TRUE }
        if (is.null(r$estimator))   { r$estimator   <- estimator; added <- TRUE }
        if (is.null(r$outer_opt))   { r$outer_opt   <- outer_opt; added <- TRUE }
        if (added) {
          saveRDS(r, new_f)
          n_backfilled <- n_backfilled + 1L
        }
      }
    }
  }

  if (verbose) {
    message(sprintf(
      "migrate_scm_output: copied=%d backfilled=%d skipped=%d (target: %s)",
      n_copied, n_backfilled, n_skipped, new_root
    ))
  }
  invisible(list(copied = n_copied, backfilled = n_backfilled,
                 skipped = n_skipped))
}

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) == 0L) {
  migrate_scm_output()
}
