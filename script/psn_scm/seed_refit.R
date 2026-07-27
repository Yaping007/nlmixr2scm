# ==============================================================================
# seed_refit.R  --  update_inits replacement for the PsN 5.5 pharmpy-missing bug
# ------------------------------------------------------------------------------
# PsN's candidate .mod files carry SCREENING initial values (covariate thetas at
# ~0.001), NOT the winner's converged estimates -- those live only in the
# winner's .ext.  Normally PsN's `update_inits` (via pharmpy) injects the final
# estimates back as initials before the covariance refit; that binary is MISSING
# in this build, so `execute winner.mod` re-minimises from screening initials and
# drifts to a DIFFERENT optimum (wrong OFV -> wrong SE / condition number).
#
# This script does update_inits ourselves: it seeds the winner .mod's
# $THETA/$OMEGA/$SIGMA from the winner .ext final estimates, forces MAXEVAL=0 so
# NONMEM evaluates + computes the covariance AT the winner (OFV identical to the
# SCM winner, by construction), fixes $DATA, and appends $COVARIANCE.
#
# Usage:
#   Rscript seed_refit.R --winner_mod <path.mod> --winner_ext <path.ext> \
#       --out refit/refit.mod [--data ../data.csv] \
#       [--cov_file final_refit_cov.txt] [--maxeval 0]
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || all(is.na(a))) b else a

argv <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default = NULL) {
  i <- which(argv == flag)
  if (length(i) && i < length(argv)) argv[i + 1L] else default
}
winner_mod <- getopt("--winner_mod")
winner_ext <- getopt("--winner_ext")
out_path   <- getopt("--out")
data_rel   <- getopt("--data", "../data.csv")
cov_file   <- getopt("--cov_file", "final_refit_cov.txt")
maxeval    <- getopt("--maxeval", "0")
stopifnot(!is.null(winner_mod), !is.null(winner_ext), !is.null(out_path))
stopifnot(file.exists(winner_mod), file.exists(winner_ext))

# ---- parse winner .ext final estimates (row ITERATION == -1000000000) --------
ext_ln <- readLines(winner_ext, warn = FALSE)
hdr_i  <- grep("^\\s*ITERATION", ext_ln)[1]
hdr    <- strsplit(trimws(ext_ln[hdr_i]), "\\s+")[[1]]
rows   <- ext_ln[(hdr_i + 1L):length(ext_ln)]
mat    <- do.call(rbind, lapply(rows, function(r)
  as.numeric(strsplit(trimws(r), "\\s+")[[1]])))
colnames(mat) <- hdr
fin <- mat[mat[, "ITERATION"] == -1000000000, ]
if (is.matrix(fin)) fin <- fin[1, ]
val  <- function(nm) unname(fin[nm])
have <- function(nm) all(nm %in% names(fin))

fmt <- function(x) formatC(x, format = "g", digits = 8)
num_re <- "[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?"

# split a data line into code + trailing comment
split_cmt <- function(line) {
  if (grepl(";", line)) list(code = sub(";.*$", "", line),
                             cmt  = sub("^[^;]*", "", line))
  else                  list(code = line, cmt = "")
}

# Replace the INIT term of a single-theta line, preserving (low,.,up) / FIX.
seed_theta_line <- function(line, v) {
  sc <- split_cmt(line); code <- sc$code
  if (grepl("\\bFIX\\b", code, ignore.case = TRUE)) return(line)      # fixed: leave
  if (grepl("\\([^,]*,[^,]*,[^,]*\\)", code)) {                       # (low, init, up)
    code <- sub(paste0("(\\([^,]*,)\\s*", num_re, "(\\s*,)"),
                paste0("\\1", fmt(v), "\\3"), code)
  } else if (grepl("\\([^,]*,[^,]*\\)", code)) {                      # (low, init)
    code <- sub(paste0("(\\([^,]*,)\\s*", num_re, "(\\s*\\))"),
                paste0("\\1", fmt(v), "\\3"), code)
  } else {                                                           # bare init
    code <- sub(num_re, fmt(v), code)
  }
  paste0(code, sc$cmt)
}

# Replace the first length(vals) numeric tokens in a line's code with new values.
# Uses digit-free private-use sentinels so an inserted value can't be re-matched
# as the next numeric token.
seed_numbers <- function(line, vals) {
  sc <- split_cmt(line); code <- sc$code
  sent <- vapply(seq_along(vals), function(i) intToUtf8(0xE000 + i), "")
  for (i in seq_along(vals)) code <- sub(num_re, sent[i], code)
  for (i in seq_along(vals)) code <- sub(sent[i], fmt(vals[i]), code, fixed = TRUE)
  paste0(code, sc$cmt)
}

# ---- read winner model from the first $PROB/$SIZES ---------------------------
mln <- readLines(winner_mod, warn = FALSE)
p0  <- grep("^\\s*\\$(PROB|SIZES)", mln)[1]
mln <- mln[p0:length(mln)]

# ---- walk the model, seeding $THETA / $OMEGA / $SIGMA in parameter order ------
out <- character(0)
sec <- ""; tk <- 0L; ok <- 0L
for (line in mln) {
  if      (grepl("^\\s*\\$THETA", line)) sec <- "THETA"
  else if (grepl("^\\s*\\$OMEGA", line)) sec <- "OMEGA"
  else if (grepl("^\\s*\\$SIGMA", line)) sec <- "SIGMA"
  else if (grepl("^\\s*\\$[A-Z]", line)) sec <- "other"

  is_theta_data <- sec == "THETA" && grepl("[0-9]", line)
  is_omega_data <- sec == "OMEGA" && grepl("[0-9]", line) && !grepl("BLOCK", line)
  is_sigma_data <- sec == "SIGMA" && grepl("[0-9]", line)

  if (is_theta_data) {
    tk <- tk + 1L
    nm <- paste0("THETA", tk)
    if (nm %in% names(fin)) line <- seed_theta_line(line, val(nm))
  } else if (is_omega_data) {
    ok <- ok + 1L
    if (ok == 1L && have("OMEGA(1,1)")) {
      line <- seed_numbers(line, val("OMEGA(1,1)"))
    } else if (ok == 2L && have(c("OMEGA(2,1)", "OMEGA(2,2)"))) {
      line <- seed_numbers(line, c(val("OMEGA(2,1)"), val("OMEGA(2,2)")))
    }
  } else if (is_sigma_data && have("SIGMA(1,1)")) {
    line <- seed_numbers(line, val("SIGMA(1,1)"))
  }

  # force MAXEVAL=0: evaluate + covariance AT the seeded winner (OFV == winner)
  if (grepl("^\\s*\\$EST", line)) {
    if (grepl("MAXEVAL\\s*=", line, ignore.case = TRUE)) {
      line <- sub("MAXEVAL\\s*=\\s*[0-9]+", paste0("MAXEVAL=", maxeval), line, ignore.case = TRUE)
    } else {
      line <- paste0(line, " MAXEVAL=", maxeval)
    }
  }

  # fix $DATA path to the cell's data.csv (winner lived deep in scm_dir)
  if (grepl("^\\s*\\$DATA", line)) {
    line <- sub("(\\$DATA\\s+)(\\S+)", paste0("\\1", data_rel), line)
  }

  out <- c(out, line)
}

# ---- append $COVARIANCE (from template if present, else R,S sandwich) ---------
if (!any(grepl("^\\s*\\$COV", out))) {
  cov_line <- if (file.exists(cov_file)) paste(readLines(cov_file, warn = FALSE), collapse = "\n")
              else "$COVARIANCE UNCONDITIONAL MATRIX=RSR PRINT=E"
  out <- c(out, "", cov_line)
}

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
writeLines(out, out_path)
cat(sprintf("[seed_refit] %s (winner OFV=%s) -> %s  MAXEVAL=%s\n",
            basename(winner_mod), fmt(val("OBJ") %||% NA_real_), out_path, maxeval))
