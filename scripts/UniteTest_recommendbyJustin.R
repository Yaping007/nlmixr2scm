library(nlmixr2)
library(nlmixr2data)
library(nlmixr2scm)
devtools::load_all("C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm")

pkdata <- nlmixr2data::warfarin[nlmixr2data::warfarin$dvid == "cp", ]
warf_pk <- function() {
  ini({
    tka <- log(1.15) # log absorption rate constant (h^-1)
    tcl <- log(0.135) # log clearance (L/h)
    tv <- log(7.0) # log volume of distribution (L)
    eta.ka ~ 0.40
    eta.cl ~ 0.25
    eta.v ~ 0.10
    prop.err <- 0.10
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ prop(prop.err)
  })
}

fit_base <- nlmixr2(
  warf_pk, pkdata,
  est = "focei",
  control = nlmixr2est::foceiControl(print = 0),
  table = nlmixr2est::tableControl(cwres = TRUE)
)

#-----prespecified explicit route----------
##1. Call .buildPairs() — creates the "long" pairs table with all combinations of var-covar-shape
pairsVec_step1 <- list(
  list(var = "cl", covar = "wt",  shapes = c("power", "lin")),
  list(var = "v",  covar = "wt",  shapes = c("power", "lin")),
  list(var = "cl", covar = "sex", shapes = "cat"),
  list(var = "v",  covar = "sex", shapes = "cat")
)
pairs_v2 <- nlmixr2scm:::buildPairs(pairsVec = pairsVec_step1)
pairs_v2

#2 Call .enrichPairs() — adds data-derived metadata + categorical expansion
pairs_step2 <- nlmixr2scm:::.enrichPairs(
  pairs_v2,
  data         = pkdata,
  catvarsVec   = catvarsVec,
  missingToken = NA
)
pairs_step2

#3. Call .expandShapes() — the row-multiplier (where covar gets renamed)
pairs_step3 <- nlmixr2scm:::.expandShapes(
  pairs_step2,
  shapes       = c("power", "lin"),
  customShapes = NULL,
  inits        = NULL
)
pairs_step3[, c("var", "covar", "shape", "covExpr", "init", "lower", "upper")]
pairs_step3




#---------AUTO ROUTE with categorical covariate [passed]----------
# Step 1: buildPairs with varsVec + covarsVec (no categorical yet)
varsVec    <- c("cl", "v")
covarsVec  <- "wt"
catvarsVec <- "sex"
covarsVec <- c(covarsVec, catvarsVec) # ← THIS IS THE KEY: include catvarsVec in covarsVec for the auto route to know

pairs_step1_auto <- nlmixr2scm:::buildPairs(
  varsVec   = varsVec,
  covarsVec = covarsVec
)
cat("Step 1 — buildPairs (auto):\n")
print(pairs_step1_auto)

# ============================================================================
# Step2 makeSCMData() trace — prepares categorical variables for SCM and determines reference levels
# ============================================================================
args(nlmixr2scm:::.makeSCMData) #function (data, catvarsVec, fit, catCutoff = 0.05, covarsVec = NULL) 

# Call makeSCMData on the base fit
scm_data <- nlmixr2scm:::.makeSCMData(pkdata, catvarsVec = "sex", fit_base)
cat("\nOutput columns:", paste(colnames(scm_data$data), collapse = ", "), "\n")
cat("catLevels :", deparse(scm_data$catLevels), "\n")
cat("catRef    :", deparse(scm_data$catRef), "\n")
cat("catDropped:", deparse(scm_data$catDropped), "\n")

# Prepare SCM data first so catLevels (reference levels) are consistent
# between .enrichPairs() and the indicator columns in scm_data$data.
# Without catLevels, .enrichPairs() uses alphabetical reference (sex_male),
# but .makeSCMData() uses most-frequent reference (sex_female) → mismatch.
scm_data_auto <- nlmixr2scm:::.makeSCMData(pkdata, catvarsVec = catvarsVec, fit_base)
cat("catRef from .makeSCMData():", deparse(scm_data_auto$catRef), "\n")

# Step 3: enrichPairs — pass catLevels so reference level matches scm_data
pairs_step2_auto <- nlmixr2scm:::.enrichPairs(
  pairs_step1_auto,
  data         = pkdata,
  catvarsVec   = catvarsVec,
  missingToken = NA,
  catLevels    = scm_data_auto$catLevels  # ensures same non-ref levels as scm_data$data
)
cat("\nStep 2 — enrichPairs (now includes sex):\n")
print(pairs_step2_auto)

# Step 4: expandShapes
pairs_step3_auto <- nlmixr2scm:::.expandShapes(
  pairs_step2_auto,
  shapes       = c("power", "lin", "log", "identity"),
  customShapes = NULL,
  inits        = NULL
)
pairs_step3_auto

# ============================================================================
# Step 4 — .rebuildUiFromPairs(): inspect the generated UI
#   Pick exactly ONE continuous and ONE categorical row from pairs_step3_auto
#   so we can inspect the generated model body in isolation.
# ============================================================================

# Select one power-shape row (continuous) and one cat row (categorical)
rebuild_pairs <- pairs_step3_auto[
  pairs_step3_auto$var == "cl" &
    pairs_step3_auto$covar %in% c("wt_power", "sex_male"),
  c("var", "covar", "shape", "covExpr", "init", "lower", "upper")
]
cat("Pairs fed into .rebuildUiFromPairs():\n")
print(rebuild_pairs)

# .rebuildUiFromPairs() needs the clean base UI, not a fit object
base_ui <- fit_base$ui
rebuilt_ui <- nlmixr2scm:::.rebuildUiFromPairs(base_ui, rebuild_pairs)

#All passed as expected


# ============================================================================
# .fitCandidatePairs() trace — fits one candidate and returns stats
#   Using a single forward step: cl ~ wt_power (continuous, no context)
# ============================================================================
cat("\nFunction signature:\n")
args(nlmixr2scm:::.fitCandidatePairs)

# Pick ONE candidate row: cl ~ wt (power shape)
candidate_pair <- pairs_step3_auto[
  pairs_step3_auto$var == "cl" & pairs_step3_auto$covar == "wt_power",
  c("var", "covar", "shape", "covExpr", "init", "lower", "upper")
]
cat("\nCandidate pair fed to .fitCandidatePairs():\n")
print(candidate_pair)

# Use scm_data_auto$data — indicator columns match catLevels used in enrichPairs
# (catRef = male → indicator column is sex_female, not sex_male)
cat("\nFitting candidate — this calls nlmixr2() once...\n")
fit_results <- nlmixr2scm:::.fitCandidatePairs(
  pairs         = candidate_pair,
  base_ui       = base_ui,
  context_pairs = NULL,    # first step: no covariates included yet
  fit           = fit_base,
  data          = scm_data_auto$data,
  pVal          = 0.05,
  stepIdx       = 1L,
  add           = TRUE     # forward step
)

# Result is a list with one element per candidate
res1 <- fit_results[[1]]

cat("\n--- Result structure ---\n")
cat("Names:", paste(names(res1), collapse = ", "), "\n")
cat("var  :", res1$var, "\n")
cat("covar:", res1$covar, "\n")

cat("\n--- stats data frame ---\n")
print(t(res1$stats))  # transpose for readability

cat("\n--- Key metrics ---\n")
cat(sprintf("  dOFV        : %.3f\n", res1$stats$deltObjf))
cat(sprintf("  p-value     : %.4f\n", res1$stats$pchisqr))
cat(sprintf("  AIC (cand)  : %.3f  vs base: %.3f\n", res1$stats$AIC, fit_base$AIC))
cat(sprintf("  covarEffect : %.4f\n", res1$stats$covarEffect))
cat(sprintf("  bsvReduction: %.4f\n", res1$stats$bsvReduction))
cat(sprintf("  significant : %s (pVal = 0.05)\n",
            if (res1$stats$pchisqr < 0.05) "YES" else "NO"))


# ---- Candidate 2: sex_female ~ cl (categorical, forward, no context) ----
# catRef = "male" (most frequent), so non-reference indicator is sex_female
cat("\n\n--- Candidate 2: sex_female ~ cl (categorical) ---\n")

candidate_sex <- pairs_step3_auto[
  pairs_step3_auto$var == "cl" & pairs_step3_auto$covar == "sex_female",
  c("var", "covar", "shape", "covExpr", "init", "lower", "upper")
]
cat("Candidate pair:\n")
print(candidate_sex)

fit_results_sex <- nlmixr2scm:::.fitCandidatePairs(
  pairs         = candidate_sex,
  base_ui       = base_ui,
  context_pairs = NULL,
  fit           = fit_base,
  data          = scm_data_auto$data,
  pVal          = 0.05,
  stepIdx       = 1L,
  add           = TRUE
)

res2 <- fit_results_sex[[1]]

cat("\n--- stats data frame ---\n")
print(t(res2$stats))

cat("\n--- Key metrics ---\n")
cat(sprintf("  dOFV        : %.3f\n", res2$stats$deltObjf))
cat(sprintf("  p-value     : %.4f\n", res2$stats$pchisqr))
cat(sprintf("  AIC (cand)  : %.3f  vs base: %.3f\n", res2$stats$AIC, fit_base$AIC))
cat(sprintf("  covarEffect : %.4f\n", res2$stats$covarEffect))
cat(sprintf("  bsvReduction: %.4f\n", res2$stats$bsvReduction))
cat(sprintf("  significant : %s (pVal = 0.05)\n",
            if (res2$stats$pchisqr < 0.05) "YES" else "NO"))


# ---- Side-by-side comparison ----
cat("\n\n=== COMPARISON: wt_power ~ cl vs sex_female ~ cl ===\n")
comp <- rbind(res1$stats, res2$stats)[, c("covar", "var", "shape", "deltObjf", "pchisqr", "AIC", "covarEffect", "bsvReduction")]
print(comp)