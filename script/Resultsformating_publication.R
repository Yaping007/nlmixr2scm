##### 
# This file is to present and visualize the aggreated results for publication
# Output files cohorts(N)*scenarios*boundaries*estimation method*outoptimizer


#1. Aggregate the refit results for N300
setwd("C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm")
source("script/aggregate_refit_results.R")

res <- aggregate_refit_run(
  root    = "output",
  out_dir = "output/refit_aggregated_N300",
  cohorts = "N300"
)

# Quick sanity: cell coverage
dplyr::count(res$diag_long, scenario, boundary) |> print(n = Inf)

# Table 3 view
table3_diag_wide(res$diag_rates, "N300",
                 c("MinSuc_pct","CovStep_pct","PhysBnd_pct","MedCN","StCN"))

#2. bash commmand
#cohorts(N)*scenarios*boundaries*
#Rscript script/aggregate_refit_results.R --root output

