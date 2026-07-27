# ==============================================================================
# parse_pilot_scn16_N300.R  --  parse the 5 pilot cells into records/
# ------------------------------------------------------------------------------
# Run from the repo root AFTER the LSF jobs finish:
#   Rscript script/psn_scm/parse_pilot_scn16_N300.R
#
# Loops parse_psn_scm.R over scn16 / N300 / ds 1..5, writing
# output/psn_scm/records/N300/scn16/ds00{1..5}/psn_scm_record.{rds,meta.json},
# then prints a one-line summary per cell (selected set, objf, cond#, wall).
# ==============================================================================
suppressPackageStartupMessages(library(jsonlite))

source("script/psn_scm/parse_psn_scm.R")   # defines build_record(), %||%

N        <- 300L
SCEN     <- 16L
DS_LIST  <- 1:5
runs     <- file.path("output", "psn_scm", "runs")
tp       <- readRDS("Inputdataset/true_params_long.rds")

summ <- lapply(DS_LIST, function(ds) {
  cell <- file.path(runs, sprintf("N%d", N), sprintf("scn%02d", SCEN),
                    sprintf("ds%03d", ds))
  if (!dir.exists(cell)) {
    message(sprintf("skip ds%03d (no cell dir)", ds)); return(NULL)
  }
  rec <- tryCatch(build_record(cell, N, SCEN, ds, tp),
                  error = function(e) { message(sprintf("ds%03d FAILED: %s", ds,
                    conditionMessage(e))); NULL })
  if (is.null(rec)) return(NULL)

  out_dir <- sub("([\\\\/])runs([\\\\/])", "\\1records\\2", cell)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(rec, file.path(out_dir, "psn_scm_record.rds"))

  sel <- paste(sprintf("%s~%s(%s)", rec$scm$selected$var,
                       rec$scm$selected$covar, rec$scm$selected$shape),
               collapse = "; ")
  data.frame(
    dataset   = ds,
    n_sel     = nrow(rec$scm$selected),
    objf      = round(rec$objf, 2),
    converged = rec$converged,
    cond_num  = rec$cond_num_cor,
    wall_sec  = rec$runtime$total_sec,
    selected  = sel,
    stringsAsFactors = FALSE
  )
})

res <- do.call(rbind, summ)
cat("\n=== pilot summary (scn16 N300 ds1-5) ===\n")
print(res, row.names = FALSE)
res
