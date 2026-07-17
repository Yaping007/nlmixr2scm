# ==============================================================================
# output_schema.R
# ------------------------------------------------------------------------------
# Single source of truth for the per-fit result record written by every
# benchmark driver (refit_one_dataset.R, bench_refit_estimators.R, and the
# SCM driver). Factoring the record here guarantees the drivers cannot drift
# apart, and stamps each record with a schema_version + model_type so that
# downstream aggregators can branch safely as the schema evolves.
#
# Sourcing order (after refit_helpers.R + true_model_factory.R):
#   source("script/output_schema.R")
#
# Public API:
#   SCHEMA_VERSION                      -- character constant, bump on change
#   fit_model_type(true_mod)            -- "ode" | "linCmt" from model attr
#   assemble_common(fit, true_mod, ...) -- named list, the common record core
#   write_fit_sidecar(result, out_rds)  -- saveRDS + tiny .meta.json sidecar
#
# 2026-07-14: introduced at schema_version "2.0". Relative to the ad-hoc 1.x
# records the drivers built inline, 2.0 adds: schema_version, model_type,
# and a stable ordering of the diagnostic/estimate/error blocks. GQ is never
# part of the record -- fit$objf (FOCEi-approx OBJF) is the only likelihood
# reported, consistent with the SAEM control decision (logLik=FALSE).
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# ---- Schema version --------------------------------------------------------
# Bump when the record shape changes in a way aggregators must know about.
#   2.0 (2026-07-14): unified drivers, model_type + schema_version added.
#   2.1 (2026-07-14): dropped redundant final_est (rel_err carries estimate);
#     promoted objf/converged/cond_num_cor to top-level scalars; collapsed the
#     duplicated runtime_sec into a single fit_runtime_sec; store the raw
#     covariance matrix (cov) and back-fill parFixed$SE/%RSE from it when
#     nlmixr2 computes cov but leaves covMethod empty (the IRLS-FOCEi case);
#     the full nlmixr2 fit object is now persisted to a parallel *.fit.rds.
SCHEMA_VERSION <- "2.1"

# ---- %||% shim (NULL / NA -> fallback) -------------------------------------
`%||%` <- function(x, y) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) y else x
}

# ---- Model-type tag --------------------------------------------------------
# true_model_factory.R stamps attr(fn, "structure") = "ode" | "linCmt".
# Older models (built before the structure arg) carry no attribute; treat
# those as "linCmt" for backward compatibility.
fit_model_type <- function(true_mod) {
  st <- attr(true_mod, "structure")
  if (is.null(st) || is.na(st)) "linCmt" else as.character(st)
}

# ---- Covariance-aware parFixed --------------------------------------------
# nlmixr2 sometimes computes a valid covariance matrix (fit$cov, correctly
# dimnamed and PD) yet leaves fit$env$covMethod empty, so it never propagates
# SE / %RSE into fit$parFixedDf -- every SE comes back NA. This is the
# IRLS-FOCEi (muModel="irls") behaviour on the ODE model: the 110 s covariance
# step DID run, the matrix is good (cond_num_cor recoverable via cov2cor), but
# the SE table is unpopulated. Rather than lose the standard errors, back-fill
# them from sqrt(diag(cov)) by name (fixed thetas like lTVKA are absent from
# cov and correctly stay NA). If nlmixr2 already populated SE, leave it be.
.fill_parFixed_se <- function(fit) {
  pf <- tryCatch(as.data.frame(fit$parFixedDf), error = function(e) NULL)
  if (is.null(pf)) return(NULL)
  cov <- tryCatch(fit$cov, error = function(e) NULL)
  needs_fill <- !is.null(cov) && "SE" %in% names(pf) && all(is.na(pf$SE))
  if (needs_fill) {
    se  <- sqrt(diag(cov))
    hit <- intersect(rownames(pf), names(se))
    if (length(hit) > 0L) {
      pf[hit, "SE"] <- se[hit]
      if ("%RSE" %in% names(pf)) {
        est <- pf[hit, "Estimate"]
        pf[hit, "%RSE"] <- ifelse(is.finite(est) & est != 0,
                                  100 * abs(se[hit] / est), NA_real_)
      }
      attr(pf, "se_source") <- "cov2se"  # SE reconstructed from fit$cov
    }
  }
  pf
}

# ---- Common record core ----------------------------------------------------
# Builds the estimator-agnostic portion of a result record. Drivers add their
# own identity keys (scenario_id, dataset_id, estimator, ...) around this.
#
#   fit          : the (possibly post-refit) nlmixr2 fit object
#   true_mod     : the true-model function (carries structure + bounds attrs)
#   scenario_id  : integer scenario for the rel-err join
#   true_params  : long-format true parameter table (from refit_helpers.R)
#   runtime_sec  : wall-clock seconds for the primary fit
#   status       : "ok" | "error"
#   estimator    : grid label ("saem","focei","foceif","irlsfoceif"); passed to
#                  diagnose_fit_table3 so SAEM's FOCEi-only flags read N/A. When
#                  NULL it is auto-detected from fit$est.
#
# Returns a named list; NULL fit yields an all-NA skeleton so the aggregators
# see a consistent shape even for failed cells.
assemble_common <- function(fit, true_mod, scenario_id, true_params,
                            runtime_sec = NA_real_, status = "ok",
                            estimator = NULL) {
  bounds_spec <- attr(true_mod, "bounds_spec")
  model_type  <- fit_model_type(true_mod)

  if (is.null(fit)) {
    diag_base <- diagnose_fit(NULL)
    return(list(
      schema_version = SCHEMA_VERSION,
      model_type     = model_type,
      objf           = diag_base$objf,
      converged      = diag_base$converged,
      cond_num_cor   = diag_base$cond_num_cor,
      rel_err        = NULL,
      diag           = diag_base,
      diag_t3        = diagnose_fit_table3(NULL, bounds_spec = bounds_spec,
                                           estimator = estimator),
      parFixed       = NULL,
      cov            = NULL,
      bounds_spec    = bounds_spec,
      cont_params    = attr(true_mod, "estimated_cont_params"),
      cat_params     = attr(true_mod, "categorical_params"),
      fn_text        = attr(true_mod, "fn_text"),
      fit_runtime_sec = runtime_sec,
      status         = status
    ))
  }

  diag_t3   <- diagnose_fit_table3(fit, bounds_spec = bounds_spec,
                                   estimator = estimator)
  diag_base <- diagnose_fit(fit)
  final_est <- extract_params_long(fit)
  rel_err   <- rel_err_one(final_est, true_params, scenario_id)
  parFixed  <- .fill_parFixed_se(fit)
  cov       <- tryCatch(fit$cov, error = function(e) NULL)

  list(
    schema_version = SCHEMA_VERSION,
    model_type     = model_type,
    objf           = diag_base$objf,          # headline OBJF (FOCEi-approx)
    converged      = diag_base$converged,     # PMx-strict convergence flag
    cond_num_cor   = diag_base$cond_num_cor,  # correlation-matrix cond number
    rel_err        = rel_err,                 # carries estimate + true + rel_err
    diag           = diag_base,
    diag_t3        = diag_t3,
    parFixed       = parFixed,                # SE back-filled from cov if needed
    cov            = cov,                      # raw covariance matrix (dimnamed)
    bounds_spec    = bounds_spec,
    cont_params    = attr(true_mod, "estimated_cont_params"),
    cat_params     = attr(true_mod, "categorical_params"),
    fn_text        = attr(true_mod, "fn_text"),
    fit_runtime_sec = runtime_sec,
    status         = status
  )
}

# ---- Sidecar writer --------------------------------------------------------
# Writes THREE artefacts next to each other, all derived from out_rds:
#   <name>.rds       -- the full record list (heavy tibbles + cov matrix)
#   <name>.fit.rds   -- the raw nlmixr2 fit object (only if `fit` supplied);
#                       the authoritative source for cov/SE/residuals/tables so
#                       downstream code can recompute anything from the fit.
#   <name>.meta.json -- a tiny greppable manifest of scalar identity + status
#                       + headline numbers (objf, converged, cond_num_cor) so a
#                       directory of results can be triaged without loading RDS.
#
#   result  : the full named list to persist (must be a list)
#   out_rds : target .rds path (the sidecar paths are derived from it)
#   fit     : optional nlmixr2 fit object; when non-NULL it is saved to the
#             parallel *.fit.rds. Drivers pass it only on a successful fit.
#
# The heavy tibbles (rel_err, parFixed, cov) stay in the record RDS; the fit
# object stays in the .fit.rds; the meta holds only flat scalars.
write_fit_sidecar <- function(result, out_rds, fit = NULL) {
  stopifnot(is.list(result), is.character(out_rds), length(out_rds) == 1L)
  saveRDS(result, out_rds)

  # Full fit object -> parallel *.fit.rds (authoritative cov/SE/tables).
  fit_rds <- NA_character_
  if (!is.null(fit)) {
    fit_rds <- sub("\\.rds$", ".fit.rds", out_rds)
    if (identical(fit_rds, out_rds)) fit_rds <- paste0(out_rds, ".fit.rds")
    fit_rds <- tryCatch({ saveRDS(fit, fit_rds); fit_rds },
                        error = function(e) NA_character_)
  }

  meta_keys <- c("schema_version", "model_type", "sample_N", "scenario_id",
                 "dataset_id", "estimator", "outer_opt", "refit_estimator",
                 "boundary", "fit_runtime_sec", "refit_runtime_sec",
                 "objf", "converged", "cond_num_cor", "status")
  meta <- result[intersect(meta_keys, names(result))]

  meta_json <- sub("\\.rds$", ".meta.json", out_rds)
  if (identical(meta_json, out_rds)) meta_json <- paste0(out_rds, ".meta.json")

  ok <- tryCatch({
    con <- file(meta_json, open = "wt", encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    writeLines(.meta_to_json(meta), con)
    TRUE
  }, error = function(e) FALSE)

  invisible(list(rds  = out_rds,
                 fit  = fit_rds,
                 meta = if (ok) meta_json else NA_character_))
}

# ---- Minimal JSON encoder (no jsonlite dependency on HPCE) -----------------
# Handles the flat scalar list produced above: character, numeric, integer,
# logical, NA. Not a general encoder -- deliberately tiny and dependency-free.
.meta_to_json <- function(x) {
  esc <- function(s) gsub('"', '\\\\"', s, fixed = TRUE)
  enc_val <- function(v) {
    if (length(v) != 1L) v <- v[1L]
    if (is.null(v) || is.na(v))       return("null")
    if (is.logical(v))                return(if (v) "true" else "false")
    if (is.numeric(v))                return(format(v, scientific = FALSE,
                                                    trim = TRUE))
    paste0('"', esc(as.character(v)), '"')
  }
  if (length(x) == 0L) return("{}")
  body <- vapply(names(x), function(k)
    sprintf('  "%s": %s', k, enc_val(x[[k]])), character(1))
  paste0("{\n", paste(body, collapse = ",\n"), "\n}")
}
