library(nlmixr2)
library(nlmixr2data)
library(nlmixr2scm)
devtools::load_all("C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm")

pkdata <- nlmixr2data::warfarin[nlmixr2data::warfarin$dvid == "cp", ]

## The example varsVec / covarsVec / catvarsVec / shapes (from README) - 
## these are exactly what the user-facing call passes to the pipeline.
varsVec    <- c("cl", "v")
covarsVec  <- "wt"
catvarsVec <- "sex"
shapes     <- c("power", "lin")

pairs_step1 <- nlmixr2scm:::buildPairs(
  varsVec   = varsVec,
  covarsVec = covarsVec
)
pairs_step1

#pre-build a pairsVec:
pairsVec_step1 <- list(
  list(var = "cl", covar = "wt",  shapes = c("power", "lin")),
  list(var = "v",  covar = "wt",  shapes = c("power", "lin")),
  list(var = "cl", covar = "sex", shapes = "cat"),
  list(var = "v",  covar = "sex", shapes = "cat")
)
pairs_v2 <- nlmixr2scm:::buildPairs(pairsVec = pairsVec_step1)
pairs_v2

# Call .enrichPairs() — adds data-derived metadata + categorical expansion
pairs_step2 <- nlmixr2scm:::.enrichPairs(
  pairs_v2,
  data         = pkdata,
  catvarsVec   = catvarsVec,
  missingToken = NA,
  catLevels    = NULL
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

#inspect function 
nlmixr2scm:::.fitCandidatePairs   # paste this and read top 30 lines