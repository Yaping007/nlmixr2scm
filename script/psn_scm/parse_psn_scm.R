# ==============================================================================
# parse_psn_scm.R  --  PsN scm output -> schema-2.1 record
# ------------------------------------------------------------------------------
# Parses one PsN `scm` cell directory (see script/psn_scm/OUTPUT_SCHEMA.md) into
# the SAME schema-2.1 fields produced by package_scm_schema21() on the nlmixr2
# side, so PsN and focei_bobyqa records feed the identical compute_scm_*
# aggregators.
#
# Sources parsed per cell (cell = output/psn_scm/runs/N*/scn*/ds*):
#   logs/scmlog.txt (preferred) | scm_dir/scmlog.txt | logs/scm_console.log
#                                                                 -> selection path
#   scm_dir/**/<final model>.{mod,lst,ext}                        -> theta, OBJ, status
#   logs/refit.{lst,ext} (optional, from STAGE 2 refit)           -> SE, cond#
#   logs/timing.json  (scope=base+scm[+refit], wall/cpu/hog)      -> runtime
#   logs/base.{lst,ext} (optional)                                -> base OFV
#
# The light record is written to output/psn_scm/records/N*/scn*/ds*/ by default.
#
# Covariate/shape handling (see OUTPUT_SCHEMA.md sec 3):
#   power  state5  (cov/med)**theta          -> nlmixr2 "power", no transform
#   exp    state4  exp(theta*(cov-med))       -> nlmixr2 "lin",   no transform
#   cat    state2  (1 + theta*I)              -> nlmixr2 "cat" exp form:
#                    theta_nlmixr = log(1+theta_psn); SE = SE_psn/(1+theta_psn)
#
# PsN covar/param name -> truth parameter name in true_params_long:
#   CLBW->CLBW, CLCRCL->CLcrCL, V2BW->VcBW, V2SEX->VcSEX
#   TVCL->TVCL, TVV2->TVVc, TVQ->TVQ, TVV3->TVVp, TVKA->TVKA
#
# Usage:
#   Rscript script/psn_scm/parse_psn_scm.R \
#       --cell output/psn_scm/runs/N300/scn16/ds001 \
#       --N 300 --scenario 16 --dataset 1 \
#       --true_params Inputdataset/true_params_long.rds \
#       [--out_rds <path>]  # default: cell path with runs/ -> records/
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

CN_STRICT_CUTOFF <- 1000

# ---- expected centering anchors for continuous covariates -------------------
# PsN's scm is configured (export_one_dataset.R [code] section) to pin BW at 70
# and CrCL at 95, so TVCL/TVV2 are estimated DIRECTLY at these fixed references
# and no post-hoc re-centering is applied.  These values are kept only to VALIDATE
# that the winner .mod actually carries the expected center (a mismatch means the
# [code] section didn't take effect); see the check in build_record().  BMI and
# the categoricals are not anchored and are not checked here.
REFERENCE_VALUES <- c(BW = 70, CRCL = 95)

# ---- name maps -------------------------------------------------------------
# PsN model param name -> (nlmixr2 var, truth beta parameter)
.PSN_COVAR_TO_TRUTH <- c(
  CLBW   = "CLBW",   CLCRCL = "CLcrCL",
  V2BW   = "VcBW",   V2SEX  = "VcSEX"
)
# PsN param stem -> nlmixr2 var label
.PSN_PARAM_VAR <- c(CL = "cl", V2 = "vc")
# structural theta name (from $THETA comment) -> truth structural parameter
.PSN_STRUCT_TO_TRUTH <- c(TVCL = "TVCL", TVV2 = "TVVc", TVQ = "TVQ",
                          TVV3 = "TVVp", TVKA = "TVKA")

# ============================================================================
# 0. small helpers
# ============================================================================

## ---- condition number from a NONMEM "EIGENVALUES OF COR MATRIX" block -------
# Eigenvalues are printed in scientific notation (e.g. 2.15E-01); the index
# rows (1 2 3 ...) are plain integers, so grab only sci-notation tokens.
cond_from_eigenvalues <- function(ln) {
  ei <- grep("EIGENVALUES OF (COR|CORRELATION) MATRIX", ln, ignore.case = TRUE)
  if (!length(ei)) return(NA_real_)
  # Eigenvalues print within a short window after the header; the index row
  # (1 2 3 ...) uses plain integers, so keeping only sci-notation tokens
  # (e.g. 2.15E-01) naturally excludes it.
  blk <- ln[(ei[1] + 1L):min(length(ln), ei[1] + 12L)]
  ev  <- suppressWarnings(as.numeric(unlist(
    regmatches(blk, gregexpr("[0-9]\\.[0-9]+[eE][+-][0-9]+", blk)))))
  ev  <- ev[is.finite(ev) & ev > 0]
  if (length(ev) < 2L) return(NA_real_)
  max(ev) / min(ev)
}

## ---- GNU /usr/bin/time -v raw file -> wall/cpu/hog -------------------------
# Fallback when no timing.json is present; reads the STAGE time_*_raw.txt files.
parse_gnu_time <- function(path) {
  if (!file.exists(path)) return(NULL)
  ln <- readLines(path, warn = FALSE)
  g  <- function(pat) {
    x <- grep(pat, ln, value = TRUE)
    if (length(x)) trimws(sub(".*:\\s*", "", x[1])) else NA_character_
  }
  usr <- suppressWarnings(as.numeric(g("User time \\(seconds\\)")))
  sys <- suppressWarnings(as.numeric(g("System time \\(seconds\\)")))
  pct <- suppressWarnings(as.numeric(gsub("%", "", g("Percent of CPU"))))
  wl  <- grep("Elapsed \\(wall clock\\)", ln, value = TRUE)
  wr  <- if (length(wl)) trimws(sub("^.*\\):\\s*", "", wl[1])) else NA_character_
  wall <- NA_real_
  if (!is.na(wr)) {
    p <- suppressWarnings(as.numeric(strsplit(wr, ":")[[1]]))
    wall <- switch(as.character(length(p)),
                   "3" = p[1] * 3600 + p[2] * 60 + p[3],
                   "2" = p[1] * 60 + p[2],
                   "1" = p[1], NA_real_)
  }
  list(wall_sec = wall,
       cpu_sec  = sum(c(usr, sys), na.rm = TRUE),
       hog_factor = if (!is.na(pct)) pct / 100 else NA_real_)
}

# ============================================================================
# 1. Low-level file parsers
# ============================================================================

## ---- .ext: final estimates + optional SE -----------------------------------
# NM7 .ext special ITERATION codes:
#   -1000000000 final estimates ; -1000000001 standard errors
# Returns list(est = named numeric, se = named numeric or NULL, obj = numeric).
parse_ext <- function(ext_path) {
  stopifnot(file.exists(ext_path))
  ln <- readLines(ext_path, warn = FALSE)
  # header row starts with " ITERATION"
  hdr_i <- grep("^\\s*ITERATION", ln)[1]
  hdr   <- strsplit(trimws(ln[hdr_i]), "\\s+")[[1]]   # ITERATION THETA1 ... OBJ
  rows  <- ln[(hdr_i + 1L):length(ln)]
  mat   <- do.call(rbind, lapply(rows, function(r) {
    as.numeric(strsplit(trimws(r), "\\s+")[[1]])
  }))
  colnames(mat) <- hdr
  it <- mat[, "ITERATION"]
  grab <- function(code) {
    i <- which(it == code)
    if (length(i) == 0L) return(NULL)
    v <- mat[i[1], ]
    v[!(names(v) %in% c("ITERATION", "OBJ"))]
  }
  est <- grab(-1000000000)
  se  <- grab(-1000000001)          # absent when no $COV
  obj_i <- which(it == -1000000000)
  obj <- if (length(obj_i)) mat[obj_i[1], "OBJ"] else NA_real_
  list(est = est, se = se, obj = obj)
}

## ---- .lst: minimization status, OFV, condition number ----------------------
parse_lst <- function(lst_path) {
  if (!file.exists(lst_path)) {
    return(list(min_success = NA, obj = NA_real_, cond_num = NA_real_,
                rounding = NA, cov_step_ok = NA))
  }
  ln <- readLines(lst_path, warn = FALSE)
  has <- function(pat) any(grepl(pat, ln, ignore.case = TRUE))
  obj_ln <- grep("OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT", ln, value = TRUE)
  obj <- if (length(obj_ln)) {
    as.numeric(sub(".*:\\s*", "", obj_ln[length(obj_ln)]))
  } else NA_real_
  # condition number appears only when $COV ran
  cn_ln <- grep("COND(ITION)?\\.?\\s*(NR|NUMBER)", ln, value = TRUE,
                ignore.case = TRUE)
  cond_num <- if (length(cn_ln)) {
    suppressWarnings(as.numeric(gsub("[^0-9.eE+-]", "",
                                     sub(".*[: ]", "", cn_ln[length(cn_ln)]))))
  } else NA_real_
  # NONMEM with $COV PRINT=E prints "EIGENVALUES OF COR MATRIX" instead of a
  # single condition-number line; condition number = max(eig)/min(eig).
  if (!is.finite(cond_num)) {
    cond_num <- cond_from_eigenvalues(ln)
  }
  list(
    min_success = has("MINIMIZATION SUCCESSFUL"),
    obj         = obj,
    cond_num    = cond_num,
    rounding    = has("ROUNDING ERRORS"),
    cov_step_ok = has("STANDARD ERROR OF ESTIMATE") ||
                  has("COVARIANCE STEP.*(OK|SUCCESSFUL)")
  )
}

## ---- .mod: covariate definitions (theta index, shape, centering) -----------
# Reads ";;; <PARAMCOVAR>-DEFINITION START ... END" blocks and the $THETA
# comments to map THETA(k) -> (param, covar, shape, center).
# Returns tibble(theta_idx, psn_param, var, covar, shape, center).
parse_mod_covariates <- function(mod_path) {
  ln <- readLines(mod_path, warn = FALSE)
  # --- theta index from $THETA comments (order = theta number) ---
  theta_lines <- grep("\\$THETA|;", ln)               # rough; refine below
  # Collect every line that assigns a theta with a trailing "; NAME" comment.
  # PsN writes covariate thetas as:  $THETA  (-100,0.001,100000) ; CLBW1
  theta_names <- c()
  tk <- 0L
  in_theta <- FALSE
  for (l in ln) {
    if (grepl("^\\s*\\$THETA", l)) in_theta <- TRUE
    if (grepl("^\\s*\\$(OMEGA|SIGMA|PK|ERROR|EST|COV|PRED|DES|SUB)", l)) in_theta <- FALSE
    if (in_theta && grepl("\\d", l)) {
      # count theta entries on this line (handles one-per-line PsN output)
      nm <- sub(".*;\\s*", "", l)
      nm <- trimws(nm)
      # a $THETA line may hold the keyword; count numeric pared groups
      ngroup <- length(gregexpr("\\(|^\\s*[-0-9]", l)[[1]])
      tk <- tk + 1L
      theta_names[tk] <- nm
    }
  }
  # --- covariate definition blocks ---
  starts <- grep("-DEFINITION START", ln)
  defs <- lapply(starts, function(i) {
    tag <- sub(";;;\\s*", "", sub("-DEFINITION START.*", "", ln[i]))
    tag <- trimws(tag)                                # e.g. "CLBW", "V2SEX"
    end <- grep("-DEFINITION END", ln)
    end <- end[end > i][1]
    body <- paste(ln[(i + 1L):(end - 1L)], collapse = "\n")
    # shape + centering
    shape <- NA_character_; center <- NA_real_; thidx <- NA_integer_
    if (grepl("\\*\\*THETA", body)) {                 # (cov/med)**THETA(k)  power
      shape <- "power"
      m <- regmatches(body, regexpr("/[0-9.]+", body))
      center <- suppressWarnings(as.numeric(sub("/", "", m)))
      thidx <- as.integer(sub(".*THETA\\((\\d+)\\).*", "\\1", body))
    } else if (grepl("EXP\\(THETA", body) || grepl("\\+\\s*THETA.*\\*\\s*\\(", body)) {
      shape <- "exp"                                  # exp(theta*(cov-med))
      m <- regmatches(body, regexpr("-\\s*[0-9.]+", body))
      center <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", m)))
      thidx <- as.integer(sub(".*THETA\\((\\d+)\\).*", "\\1", body))
    } else if (grepl("1\\s*\\+\\s*THETA", body)) {    # (1 + THETA(k))  categorical
      shape <- "cat"
      thidx <- as.integer(sub(".*THETA\\((\\d+)\\).*", "\\1", body))
    }
    # split tag into param stem + covar: CLBW -> CL + BW ; V2SEX -> V2 + SEX
    stem <- if (startsWith(tag, "CL")) "CL" else if (startsWith(tag, "V2")) "V2" else NA
    covar <- sub("^(CL|V2)", "", tag)
    tibble(theta_idx = thidx, psn_param = stem,
           var = unname(.PSN_PARAM_VAR[stem] %||% NA),
           covar = covar, shape = shape, center = center,
           tag = tag)
  })
  out <- bind_rows(defs)
  attr(out, "theta_names") <- theta_names
  out
}

## ---- selection path from scmlog.txt or scm_console.log ----------------------
# Returns tibble(step_dir = forward|backward, action = add|remove,
#                covar, param, state, order).
#
# PsN 5.5 scmlog.txt records each accepted step as e.g.
#   "Parameter-covariate relation chosen in this forward step: CL-CRCL-4"
#   "Parameter-covariate relation chosen in this backward step: CL-SEX-1"
# where the trailing integer is the STATE (5=power, 4=exp, 2=cat, 1=removed).
# A backward step to state 1 == the relation is dropped.
parse_selection <- function(cell_dir) {
  cand <- c(file.path(cell_dir, "logs", "scmlog.txt"),
            file.path(cell_dir, "scm_dir", "scmlog.txt"),
            file.path(cell_dir, "logs", "scm_console.log"),
            file.path(cell_dir, "scm_console.log"))
  src <- cand[file.exists(cand)][1]
  if (is.na(src)) stop("no scm log found in ", cell_dir)
  ln  <- readLines(src, warn = FALSE)

  # ---- preferred: "relation chosen in this <dir> step: PARAM-COVAR-STATE" ----
  chosen <- grep("relation chosen in this (forward|backward) step", ln,
                 value = TRUE)
  if (length(chosen)) {
    m <- regmatches(chosen, regexec(
      "chosen in this (forward|backward) step:\\s*([A-Za-z0-9]+)-([A-Za-z0-9]+)-(\\d+)",
      chosen))
    df <- do.call(rbind, lapply(m, function(z) {
      if (length(z) < 5) return(NULL)
      data.frame(step_dir = z[2], param = z[3], covar = z[4],
                 state = as.integer(z[5]), stringsAsFactors = FALSE)
    }))
    return(as_tibble(df) |>
             mutate(action = if_else(state == 1L, "remove", "add"),
                    order  = row_number()))
  }

  # ---- legacy fallback: "Adding/Removing COVAR on PARAM state N" ----
  add <- grep("^Adding\\s+", ln, value = TRUE)
  rem <- grep("^Removing\\s+", ln, value = TRUE)
  parse_line <- function(x, action) {
    m <- regmatches(x, regexec(
      "(Adding|Removing)\\s+(\\S+)\\s+on\\s+(\\S+)\\s+state\\s+(\\d+)", x))
    do.call(rbind, lapply(m, function(z) {
      if (length(z) < 5) return(NULL)
      data.frame(action = action, covar = z[3], param = z[4],
                 state = as.integer(z[5]), stringsAsFactors = FALSE)
    }))
  }
  fwd <- parse_line(add, "add")
  bck <- parse_line(rem, "remove")
  bind_rows(
    if (!is.null(fwd)) mutate(as_tibble(fwd), step_dir = "forward") else NULL,
    if (!is.null(bck)) mutate(as_tibble(bck), step_dir = "backward") else NULL
  ) |>
    mutate(order = row_number())
}

# Count the number of candidate NONMEM model fits PsN performed during the SCM
# search -- the effort axis that pairs with cpu_sec.  In scmlog.txt every tested
# candidate is one "MODEL ... PVAL ..." row (both forward and backward tables);
# the base model (node-0) is fit once and is not a MODEL row, so we add 1.
# Returns a list(n_models_fit, n_forward, n_backward, n_base) for the record.
count_model_fits <- function(cell_dir) {
  cand <- c(file.path(cell_dir, "logs", "scmlog.txt"),
            file.path(cell_dir, "scm_dir", "scmlog.txt"),
            file.path(cell_dir, "logs", "scm_console.log"),
            file.path(cell_dir, "scm_console.log"))
  src <- cand[file.exists(cand)][1]
  if (is.na(src)) return(list(n_models_fit = NA_integer_, n_forward = NA_integer_,
                              n_backward = NA_integer_, n_base = 1L))
  ln <- readLines(src, warn = FALSE)

  # split into forward vs backward regions at the "Starting backward search" marker
  bk <- grep("Starting backward search", ln)[1]
  fwd_ln <- if (is.na(bk)) ln else ln[seq_len(bk - 1L)]
  bck_ln <- if (is.na(bk)) character(0) else ln[bk:length(ln)]

  # a candidate row = "TAG-STATE  PVAL  <base ofv> <new ofv> ..." (10 cols),
  # excluding the "MODEL ... TEST ..." header line.
  is_cand <- function(x) grepl("^\\s*[A-Za-z].*\\bPVAL\\b.*-?[0-9]", x) &
                          !grepl("\\bTEST\\b", x)
  n_fwd <- sum(is_cand(fwd_ln))
  n_bck <- sum(is_cand(bck_ln))
  list(n_models_fit = 1L + n_fwd + n_bck,   # +1 = base (node-0)
       n_forward = n_fwd, n_backward = n_bck, n_base = 1L)
}

# ============================================================================
# 2. Resolve the FINAL selected model + relation set
# ============================================================================
# Final relation set = forward-added relations minus backward-removed ones.
# (state carried = the last state each param-covar reached in forward.)
final_relations <- function(sel) {
  if (nrow(sel) == 0L) return(tibble())
  fwd <- filter(sel, step_dir == "forward")
  bck <- filter(sel, step_dir == "backward")
  # last state per (param,covar) in forward
  kept <- fwd |>
    group_by(param, covar) |>
    slice_max(order, n = 1L, with_ties = FALSE) |>
    ungroup()
  # remove those eliminated backward
  if (nrow(bck)) {
    kept <- anti_join(kept, distinct(bck, param, covar), by = c("param", "covar"))
  }
  kept |>
    transmute(param, covar, state,
              shape = dplyr::case_when(state == 5 ~ "power",
                                       state == 4 ~ "exp",
                                       state == 2 ~ "cat",
                                       TRUE ~ NA_character_))
}

# Find the .mod/.lst/.ext whose covariate-definition set == final relation set.
# Models are matched on tag + STATE (shape), not tag name alone: a forward
# candidate can share the same covariate NAMES as the backward winner but carry
# one of them at a different state (e.g. CLBW at exp/4 vs power/5), so name-only
# matching is ambiguous.  state encoding: power=5, exp=4, cat=2.
find_final_model <- function(cell_dir, fr) {
  if (nrow(fr) == 0L) return(NULL)                    # nothing selected -> no cov winner
  mods <- list.files(file.path(cell_dir, "scm_dir"), pattern = "\\.mod$",
                     recursive = TRUE, full.names = TRUE)
  mods <- mods[!grepl("NM_run|temp_dir|final_models", mods)]
  shape_to_state <- c(power = 5L, exp = 4L, cat = 2L)
  want_sig <- sort(paste0(fr$param, fr$covar, "-",
                          shape_to_state[fr$shape]))
  want <- paste(want_sig, collapse = "|")
  best <- NULL; best_score <- -Inf
  for (m in mods) {
    dv <- tryCatch(parse_mod_covariates(m), error = function(e) NULL)
    if (is.null(dv) || nrow(dv) == 0L) next
    have_sig <- sort(paste0(dv$tag, "-", shape_to_state[dv$shape]))
    have <- paste(have_sig, collapse = "|")
    # exact tag+state set preferred; otherwise reward matches but PENALISE
    # extra/missing signatures so neither a superset nor a wrong-state model
    # can outscore the true selected set.
    score <- if (identical(have, want)) 1e6 else {
      hit   <- sum(have_sig %in% want_sig)
      miss  <- length(setdiff(want_sig, have_sig))
      extra <- length(setdiff(have_sig, want_sig))
      hit - miss - extra
    }
    if (score > best_score) { best_score <- score; best <- m }
  }
  if (is.null(best)) return(NULL)
  # Reject a partial match that is MISSING any selected covariate (a subset
  # model, e.g. a backward candidate that dropped one relation): its estimates
  # would misrepresent the winner.  A superset is tolerated; a subset is not.
  best_sig <- sort(paste0(
    parse_mod_covariates(best)$tag, "-",
    shape_to_state[parse_mod_covariates(best)$shape]))
  if (length(setdiff(want_sig, best_sig)) > 0L) return(NULL)
  stem <- sub("\\.mod$", "", best)
  list(mod = best,
       lst = paste0(stem, ".lst"),
       ext = paste0(stem, ".ext"))
}

# Find the BASE model fit (node-0): the fitted .mod carrying NO covariate
# definitions.  Used when SCM selects nothing (null scenario) -- the base model
# IS the winner, so its OFV / condition number / structural typical values are
# the correct record, not NA.  Prefers a name containing "base", else the
# shallowest path; requires a matching .lst + .ext (an actually-fitted model).
find_base_model <- function(cell_dir) {
  mods <- list.files(file.path(cell_dir, "scm_dir"), pattern = "\\.mod$",
                     recursive = TRUE, full.names = TRUE)
  mods <- mods[!grepl("NM_run|temp_dir|final_models", mods)]
  if (!length(mods)) return(NULL)
  cand <- character(0)
  for (m in mods) {
    dv <- tryCatch(parse_mod_covariates(m), error = function(e) NULL)
    if (is.null(dv) || nrow(dv) == 0L) cand <- c(cand, m)
  }
  if (!length(cand)) return(NULL)
  pick <- cand[grepl("base", basename(cand), ignore.case = TRUE)]
  if (!length(pick)) pick <- cand[order(lengths(strsplit(cand, "[\\\\/]")))]
  for (m in pick) {
    stem <- sub("\\.mod$", "", m)
    if (file.exists(paste0(stem, ".lst")) && file.exists(paste0(stem, ".ext")))
      return(list(mod = m, lst = paste0(stem, ".lst"), ext = paste0(stem, ".ext")))
  }
  NULL
}

# ============================================================================
# 3. Build schema-2.1 record
# ============================================================================

build_record <- function(cell_dir, N, scenario, dataset, true_params,
                          estimator = "nonmem_scm", outer_opt = "focei",
                          structure = "advan4") {
  sel <- parse_selection(cell_dir)
  fr  <- final_relations(sel)
  fm  <- find_final_model(cell_dir, fr)
  # null-selection (e.g. scn01): SCM kept the BASE model -> use it as the winner
  # so OFV / condition number / structural typical values are recovered.
  base_winner <- FALSE
  if (is.null(fm) && nrow(fr) == 0L) {
    fm <- find_base_model(cell_dir)
    if (!is.null(fm)) base_winner <- TRUE
  }
  nfit <- count_model_fits(cell_dir)

  have_winner <- !is.null(fm)
  if (have_winner) {
    covdef <- parse_mod_covariates(fm$mod)
    # base-model winner (null selection) has no DEFINITION blocks -> bind_rows()
    # yields a 0-column tibble; replace with a correctly-typed empty covariate set.
    if (base_winner || ncol(covdef) == 0L) {
      covdef <- tibble(theta_idx = integer(), psn_param = character(),
                       var = character(), covar = character(),
                       shape = character(), center = numeric(), tag = character())
    }
    ext    <- parse_ext(fm$ext)
    lst    <- parse_lst(fm$lst)
  } else {
    # scm_dir winner .mod not available (e.g. partial download). Derive the
    # selected set from the log's final relations; covariate point estimates
    # are unavailable, but selection pattern / shapes / OFV / cond / timing
    # are still recovered below from the refit + log.
    warning("winner .mod not found under scm_dir; building reduced record ",
            "from scmlog final relations for ", cell_dir)
    covdef <- if (nrow(fr) == 0L) {
      # no covariate selected (e.g. null scenario): empty, correctly-typed set
      tibble(theta_idx = integer(), psn_param = character(),
             var = character(), covar = character(),
             shape = character(), center = numeric(), tag = character())
    } else {
      fr |>
        transmute(theta_idx = NA_integer_, psn_param = param,
                  var = unname(.PSN_PARAM_VAR[param]), covar = covar,
                  shape = shape, center = NA_real_, tag = paste0(param, covar))
    }
    ext <- list(est = NULL, se = NULL, obj = NA_real_)
    lst <- parse_lst(NA_character_)
  }

  # optional final refit (SE + cond#) from STAGE 2.  Prefer the light copy in
  # logs/, else search the refit/ or scm_dir trees.
  refit_lst_path <- file.path(cell_dir, "logs", "refit.lst")
  if (!file.exists(refit_lst_path)) {
    cand <- list.files(cell_dir, pattern = "refit.*\\.lst$", recursive = TRUE,
                       full.names = TRUE)
    cand <- cand[!grepl("NM_run|temp_dir", cand)]
    refit_lst_path <- if (length(cand)) cand[1] else NA_character_
  }
  refit <- if (!is.na(refit_lst_path) && file.exists(refit_lst_path)) {
    parse_lst(refit_lst_path)
  } else NULL
  refit_ext <- if (!is.na(refit_lst_path) && file.exists(refit_lst_path)) {
    e <- sub("\\.lst$", ".ext", refit_lst_path)
    if (file.exists(e)) parse_ext(e) else NULL
  } else NULL

  # ---- REFIT VALIDITY GUARD ---------------------------------------------------
  # The refit MUST be the covariance evaluation of the WINNER (MAXEVAL=0 seeded
  # at the winner .ext), so its OFV equals the winner OFV.  A stale/leftover
  # refit_run (e.g. a base-model .lst from an earlier attempt) would otherwise
  # silently supply SE + condition number for the WRONG model.  If the refit OFV
  # disagrees with the winner OFV by > 0.01, reject the refit (cov_done=FALSE)
  # rather than record a mismatched condition number.
  if (have_winner && !is.null(refit) &&
      is.finite(ext$obj %||% NA_real_) &&
      is.finite(refit$obj %||% NA_real_) &&
      abs((refit$obj %||% NA_real_) - ext$obj) > 0.01) {
    warning(sprintf(
      "refit OFV (%.4f) != winner OFV (%.4f) in %s; rejecting stale refit ",
      refit$obj, ext$obj, cell_dir),
      "(cond#/SE dropped, cov_done=FALSE). Re-run seed_refit + parser.")
    refit <- NULL; refit_ext <- NULL
  }

  # ---- selected covariates with estimates (+ categorical back-transform) ----
  th <- ext$est
  th_se <- (refit_ext %||% ext)$se
  if (have_winner) {
    selected <- covdef |>
      filter(!is.na(theta_idx)) |>
      rowwise() |>
      mutate(
        theta_name = paste0("THETA", theta_idx),
        est_psn    = unname(th[paste0("THETA", theta_idx)]),
        se_psn     = if (!is.null(th_se)) unname(th_se[paste0("THETA", theta_idx)]) else NA_real_,
        # back-transform categorical (1+theta) -> exp form log(1+theta)
        estimate   = dplyr::if_else(shape == "cat", log(1 + est_psn), est_psn),
        se         = dplyr::if_else(shape == "cat", se_psn / (1 + est_psn), se_psn),
        truth_param = unname(.PSN_COVAR_TO_TRUTH[tag])
      ) |>
      ungroup() |>
      select(var, covar, shape, theta_name, estimate, se, est_psn, center,
             truth_param)
  } else {
    # reduced record: selection pattern only, estimates unavailable
    selected <- covdef |>
      transmute(var, covar, shape,
                theta_name = NA_character_,
                estimate   = NA_real_, se = NA_real_, est_psn = NA_real_,
                center,
                truth_param = unname(.PSN_COVAR_TO_TRUTH[tag]))
  }

  # ---- structural estimates (TVCL, TVV2->TVVc, ...) ----
  # structural thetas are THETA1..THETA5 in fixed order
  struct_est <- if (!is.null(th)) unname(th[paste0("THETA", 1:5)]) else rep(NA_real_, 5)
  struct <- tibble(
    psn_name    = c("TVCL", "TVV2", "TVQ", "TVV3", "TVKA"),
    theta_name  = paste0("THETA", 1:5),
    estimate    = struct_est,
    truth_param = unname(.PSN_STRUCT_TO_TRUTH[c("TVCL","TVV2","TVQ","TVV3","TVKA")])
  )

  # ---- structural typical values are already at the fixed references ----------
  # PsN's [code] section pins BW at 70 and CrCL at 95 (see export_one_dataset.R),
  # so TVCL/TVV2 are estimated DIRECTLY at those references -- no post-hoc
  # re-centering is applied here.  (theta is center-invariant; only the intercept
  # anchor changed, and that change now happens inside NONMEM.)

  # ---- random-effect estimates (BSV variances/covariance + residual error) ----
  # Preferentially read the OMEGA/SIGMA point estimates from the REFIT .ext (the
  # MAXEVAL=0 covariance evaluation seeded at the winner): its estimates equal
  # the winner's by construction but come from the model that actually ran the
  # $COV step.  Fall back to the winner search .ext when no refit is available.
  # Model wiring (see export_one_dataset.R $PK / $OMEGA / $SIGMA):
  #   ETA(1) -> CL, ETA(2) -> V2/Vc ; $OMEGA BLOCK(2):
  #     OMEGA(1,1) = var(eta.CL)  -> var_CL    (truth 0.10, variance scale)
  #     OMEGA(2,2) = var(eta.Vc)  -> var_Vc    (truth 0.10, variance scale)
  #     OMEGA(2,1) = cov          -> cov_VcCL  (truth 0.02, covariance scale)
  #   $SIGMA(1,1) = proportional error VARIANCE (truth 0.01); the DGP defines the
  #     residual on the SD scale (prop.err = 0.1 in true_model_factory.R), so
  #     ResErr = sqrt(SIGMA(1,1)) to compare like-for-like against truth 0.1.
  th_re <- (refit_ext %||% ext)$est
  getre <- function(nm) if (!is.null(th_re) && nm %in% names(th_re))
    unname(th_re[[nm]]) else NA_real_
  sigma11 <- getre("SIGMA(1,1)")
  reff <- tibble(
    estimate    = c(getre("OMEGA(1,1)"), getre("OMEGA(2,2)"),
                    getre("OMEGA(2,1)"),
                    if (is.finite(sigma11)) sqrt(sigma11) else NA_real_),
    truth_param = c("var_CL", "var_Vc", "cov_VcCL", "ResErr")
  )

  # ---- relative error vs truth ----
  truth <- true_params |>
    filter(scenario == !!scenario) |>
    select(parameter, true_value)
  rel_err <- bind_rows(
    struct |> transmute(parameter = truth_param, estimate),
    selected |> transmute(parameter = truth_param, estimate),
    reff |> transmute(parameter = truth_param, estimate)
  ) |>
    inner_join(truth, by = "parameter") |>
    mutate(abs_err = estimate - true_value,
           rel_err = ifelse(true_value != 0, abs_err / true_value, NA_real_),
           rel_err_pct = 100 * rel_err)

  # ---- convergence / condition number ----
  # The final refit is now a MAXEVAL=0 covariance evaluation SEEDED at the winner
  # (see seed_refit.R): its OFV equals the SCM winner's by construction, and it
  # supplies SE + condition number.  But MAXEVAL=0 is NOT a minimisation, so the
  # refit .lst carries no "MINIMIZATION SUCCESSFUL" -- convergence status must
  # come from the winner's SEARCH .lst (which actually minimised).
  cond_num <- refit$cond_num %||% lst$cond_num
  if (have_winner) {
    objf      <- ext$obj                       # winner search OFV (== refit OFV)
    if (!is.finite(objf %||% NA_real_)) objf <- refit$obj %||% NA_real_
    converged <- isTRUE(lst$min_success) && is.finite(ext$obj)
    min_succ  <- isTRUE(lst$min_success)
    rounding  <- isTRUE(lst$rounding)
  } else {
    # partial download (no winner .lst): fall back to whatever the refit gives.
    objf      <- refit$obj %||% ext$obj
    converged <- isTRUE(refit$min_success) && is.finite(objf %||% NA_real_)
    min_succ  <- isTRUE(refit$min_success)
    rounding  <- isTRUE(refit$rounding)
  }

  # ---- timing (single-block total, scope=base+scm[+refit]) ----
  tj_path <- file.path(cell_dir, "logs", "timing.json")
  if (!file.exists(tj_path)) tj_path <- file.path(cell_dir, "timing.json")  # legacy
  tj <- if (file.exists(tj_path)) fromJSON(tj_path) else list()
  # fallback: GNU /usr/bin/time -v raw files (time_scm_raw.txt + time_refit_raw.txt)
  if (length(tj) == 0L) {
    t_scm   <- parse_gnu_time(file.path(cell_dir, "logs", "time_scm_raw.txt"))
    t_refit <- parse_gnu_time(file.path(cell_dir, "logs", "time_refit_raw.txt"))
    wall <- sum(c(t_scm$wall_sec, t_refit$wall_sec), na.rm = TRUE)
    cpu  <- sum(c(t_scm$cpu_sec,  t_refit$cpu_sec),  na.rm = TRUE)
    tj <- list(
      wall_sec   = if (isTRUE(wall > 0)) wall else NA_real_,
      cpu_sec    = if (isTRUE(cpu  > 0)) cpu  else NA_real_,
      hog_factor = t_scm$hog_factor %||% NA_real_,
      scope      = if (!is.null(t_refit)) "scm+refit" else if (!is.null(t_scm)) "scm" else NA_character_
    )
  }

  list(
    schema_version = "2.1",
    model_type     = structure,
    sample_N       = N,
    scenario_id    = scenario,
    dataset_id     = dataset,
    estimator      = estimator,
    outer_opt      = outer_opt,
    objf           = objf,
    converged      = converged,
    min_success    = min_succ,
    cond_num_cor   = cond_num,
    cn_below_cutoff = isTRUE(is.finite(cond_num) && cond_num <= CN_STRICT_CUTOFF),
    rounding_error = rounding,
    cov_done       = !is.null(refit) && isTRUE(refit$cov_step_ok),
    rel_err        = rel_err,
    scm = list(
      selected  = selected,
      step_hist = sel,
      final_relations = fr,
      final_model = if (have_winner) basename(fm$mod) else NA_character_,
      winner_found = have_winner && nrow(fr) > 0L,
      n_models_fit = nfit$n_models_fit,
      n_forward_fit = nfit$n_forward,
      n_backward_fit = nfit$n_backward,
      n_base_fit = nfit$n_base
    ),
    runtime = list(
      total_sec = tj$wall_sec %||% NA_real_,
      cpu_sec   = tj$cpu_sec  %||% NA_real_,
      scope     = tj$scope    %||% NA_character_
    ),
    cpu = list(
      cpu_sec    = tj$cpu_sec    %||% NA_real_,
      hog_factor = tj$hog_factor %||% NA_real_
    ),
    identity = list(sample_N = N, scenario_id = scenario, dataset_id = dataset,
                    estimator = estimator, outer_opt = outer_opt,
                    structure = structure)
  )
}

# ============================================================================
# 4. CLI
# ============================================================================
# Only run when this file is executed DIRECTLY (Rscript parse_psn_scm.R ...),
# not when it is source()d by another script (e.g. parse_pilot_scn16_N300.R).
# sys.nframe() == 0 at top level of a direct run; > 0 inside source().
if (!interactive() && sys.nframe() == 0L) {
  argv <- commandArgs(trailingOnly = TRUE)
  getopt <- function(flag, default = NULL) {
    i <- which(argv == flag)
    if (length(i) && i < length(argv)) argv[i + 1L] else default
  }
  cell   <- getopt("--cell")
  N      <- as.integer(getopt("--N", "300"))
  scn    <- as.integer(getopt("--scenario", "16"))
  ds     <- as.integer(getopt("--dataset", "1"))
  # structure tag recorded into the record (model_type / identity$structure):
  #   "advan4" (default, ADVAN4 analytic) | "ode" (ADVAN13 general-ODE).
  #   Distinguished downstream from runSCM's structure="ode" by estimator
  #   (nonmem_scm vs the runSCM label), so both ODE structures share the tag.
  structure <- tolower(getopt("--structure", "advan4"))
  tp_path <- getopt("--true_params", "Inputdataset/true_params_long.rds")
  stopifnot(!is.null(cell))
  # default record path mirrors the cell under records/ instead of runs/
  rec_dir_default <- sub("([\\\\/])runs([\\\\/])", "\\1records\\2", cell)
  out_rds <- getopt("--out_rds", file.path(rec_dir_default, "psn_scm_record.rds"))
  dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)

  true_params <- readRDS(tp_path)
  rec <- build_record(cell, N, scn, ds, true_params, structure = structure)
  saveRDS(rec, out_rds)

  # flat meta.json manifest (mirrors nlmixr2 res_ds*.meta.json headline fields)
  meta <- list(
    schema_version = rec$schema_version, model_type = rec$model_type,
    sample_N = rec$sample_N, scenario_id = rec$scenario_id,
    dataset_id = rec$dataset_id, estimator = rec$estimator,
    outer_opt = rec$outer_opt, objf = rec$objf, converged = rec$converged,
    min_success = rec$min_success, cond_num_cor = rec$cond_num_cor,
    wall_total_sec = rec$runtime$total_sec, cpu_total_sec = rec$cpu$cpu_sec,
    hog_factor = rec$cpu$hog_factor,
    n_selected = nrow(rec$scm$selected),
    selected = paste(sprintf("%s~%s(%s)", rec$scm$selected$var,
                             rec$scm$selected$covar, rec$scm$selected$shape),
                     collapse = "; ")
  )
  writeLines(toJSON(meta, auto_unbox = TRUE, pretty = TRUE, digits = 8),
             sub("\\.rds$", ".meta.json", out_rds))

  message(sprintf("[parse_psn_scm] %s -> %s", cell, out_rds))
  message(sprintf("  selected: %s", meta$selected))
  message(sprintf("  objf=%.2f converged=%s cond=%s wall=%.0fs",
                  rec$objf, rec$converged, rec$cond_num_cor %||% NA,
                  rec$runtime$total_sec %||% NA))
}
