# HPCE VAE covariate-selection pilot

One-run VAE covariate selection (`nlmixr2(..., est = "vae")`, `covariateSelection = TRUE`) across the full pilot grid on the LSF cluster. **`runSCM` is not used** — the stepwise search is internal to the VAE training loop (BICc-ELBO penalty).

------------------------------------------------------------------------

# Theory & implementation primer (for review / presentation)

> Self-contained explainer of **what a VAE is**, **why it applies to population PK/PD**, and **exactly how nlmixr2est turns a raw covariate column into a tested covariate–parameter relationship**. Each block pairs a **technical** statement with a **plain-language** gloss for Q&A.

## 1. The population PK estimation problem

**Technical.** A nonlinear mixed-effects (NLME) model describes observation $y_{ij}$ for subject $i$ at time $t_{ij}$ as

$$y_{ij} = f\!\left(t_{ij},\, \theta,\, \eta_i,\, x_i\right) + \varepsilon_{ij},
\qquad \eta_i \sim \mathcal{N}(0,\Omega),\quad
\varepsilon_{ij}\sim\mathcal{N}(0,\sigma^2).$$

$f$ is the structural (ODE / `linCmt`) model, $\theta$ the population ("typical") parameters, $\eta_i$ the subject-level random effects (BSV), $x_i$ the covariates (BW, CrCL, …). Individual parameters are $P_i = g(\theta,\eta_i,x_i)$, e.g. $\mathrm{CL}_i = \exp(\mathrm{lTVCL} + \beta\log(\mathrm{BW}_i/70) + \eta_{i,\mathrm{CL}})$. Fitting maximizes the marginal likelihood $\int p(y_i\mid\eta_i)\,p(\eta_i)\,d\eta_i$, which has **no closed form** because $f$ is nonlinear in $\eta_i$.

**Plain language.** Every patient has their own clearance and volume, drawn around a population average. We only observe drug concentrations, not those personal values. To learn the averages and the spread we must "integrate out" each patient's hidden parameters — an integral algebra can't solve. FOCEi approximates it by linearization; VAE approximates it with a neural-network encoder.

## 2. What a VAE is (variational autoencoder)

**Technical.** A VAE is a latent-variable model trained by **variational inference**. It replaces the intractable posterior $p(\eta_i\mid y_i)$ with a tractable **approximate posterior** $q_\phi(\eta_i\mid y_i)=\mathcal{N}(\mu_\phi(y_i),\Sigma_\phi(y_i))$ produced by an **encoder** network with parameters $\phi$, and maximizes the **Evidence Lower BOund (ELBO)**:

$$\mathrm{ELBO}(\theta,\phi)=
\underbrace{\mathbb{E}_{q_\phi}\!\big[\log p_\theta(y_i\mid\eta_i)\big]}_{\text{reconstruction / fit}}
-\underbrace{\mathrm{KL}\!\big(q_\phi(\eta_i\mid y_i)\,\|\,p(\eta_i)\big)}_{\text{regularizer toward }\mathcal N(0,\Omega)}.$$

Since $\log p(y_i)\ge \mathrm{ELBO}$, raising the ELBO raises the true likelihood. Gradients flow through sampling via the **reparameterization trick** $\eta_i=\mu_\phi+\Sigma_\phi^{1/2}z,\; z\sim\mathcal N(0,I)$.

-   **Encoder** ("recognition model"): in nlmixr2est an **LSTM** reading each subject's dose/observation sequence **plus covariates**, emitting $(\mu_\phi,\Sigma_\phi)$.
-   **Decoder** ("generative model"): **the ordinary `rxode2` structural model** $f$. There is **no learned decoder network and no torch** — the "decoder" is your PK ODE solved with the sampled $\eta_i$.

**Plain language.** An autoencoder squeezes data through a bottleneck and rebuilds it; here the bottleneck **is** the patient's random effects $\eta_i$. A neural net (encoder) reads a patient's curve + covariates and *guesses* their CL/V offsets; the PK model (decoder) turns that guess back into a predicted curve. Training rewards two things: predictions matching data, and guessed offsets staying believable (close to the population distribution — the KL term). The reparameterization trick lets us train through random sampling with ordinary gradient descent.

## 3. Why "one run is all you need" for covariates

**Technical.** Classic covariate building (SCM/PsN) is a **discrete outer loop**: fit → likelihood-ratio test each candidate → add/remove → re-fit, dozens of times. VAE folds covariate **selection** into the **same gradient-based training** that estimates $\theta,\Omega$. Once per latent dimension per iteration it solves an **L0-penalized model selection** of the latent means on the candidate covariate columns with a **BIC-type objective**

$$\mathcal{J}(S)=\frac{\mathrm{RSS}(S)}{\omega} + \lambda\,|S|,$$

where $S$ is the chosen subset, $\mathrm{RSS}(S)$ the residual of regressing the SAEM sufficient statistic (an EMA of posterior means $\mu_\phi$) on $[\,1\mid X_S\,]$, $\omega$ the current random-effect variance, $\lambda|S|$ the sparsity penalty. Solved **exactly by branch-and-bound** when small (\< `covSelectMaxExact` bits), else via `L0Learn` proposals re-scored exactly. Mutual exclusion (one shape per covariate per parameter) makes it stepwise-equivalent — inside one fit.

**Plain language.** SCM/PsN fit hundreds of models, one covariate at a time. VAE instead asks every step: "given my current guess of each patient's CL and V, which covariates best explain the person-to-person differences — and is each one worth its complexity cost?" It answers with a penalized regression preferring fewer covariates unless the data demand more. Estimation and selection happen together.

## 4. From a data column to a tested relationship (the covariate formula)

The part fixed in **7.0.2** and configured by these scripts.

### 4.1 Candidate discovery

**Technical.** Every non-reserved column **constant within each subject** becomes a candidate (reserved = `ID,TIME,DV,EVID,AMT,CMT,MDV,SS,II,ADDL,RATE,DUR,…`). Numeric with \> 2 unique values → **continuous**; else **categorical** (one 0/1 indicator per level, rare levels below `catCutoff` merged into reference). `vaeCovariates(data, …)` previews this without fitting.

**Plain language.** VAE scans for anything describing the *patient* (same value at every time) — weight, kidney function, sex — and lists each as testable.

### 4.2 Shapes = functional form of the effect

**Technical.** With center $c$, the coefficient $\beta$ multiplies:

| shape      | transform                                | family             |
|------------|------------------------------------------|--------------------|
| `power`    | $\beta\log(\mathrm{cov}/c)$              | log                |
| `log`      | $\beta\log(\mathrm{cov})$                | log                |
| `lin`      | $\beta(\mathrm{cov}-c)$                  | linear             |
| `identity` | $\beta\,\mathrm{cov}$                    | linear             |
| `center`   | $\beta(\mathrm{cov}/c)$                  | linear             |
| `hockey`   | two-armed piecewise-linear, knee at $c$  | hockey             |
| `cat`      | indicator $\exp(\beta I_{\text{level}})$ | categorical (auto) |

Shapes spanning the same space collapse to a **family** (`power`≡`log`; `lin`≡`identity`≡`center`), so selection sees one column per family per covariate — mirroring PsN `[valid_states]` (power = 5, exponential/`lin` = 4, categorical = 2).

**Plain language.** "Shape" is *how* a covariate bends a parameter: allometric power law $(\mathrm{BW}/70)^{\beta}$, linear, or two-slope "hockey stick". Categoricals like sex switch a parameter by a fixed factor. VAE can test several and keep the best.

### 4.3 Center / reference value

**Technical.** `power`, `lin`, `center` subtract/divide by reference $c$. `covCenterType` picks the statistic (`"median"` default / `"mean"`); `covCenter=c(BW=70, CrCL=95)` **overrides** with fixed physiological anchors (matching PsN's `[code]` section). The **coefficient** $\beta$ is centering-invariant — only the structural intercept absorbs a change of $c$ — so selection and $\beta$ stay comparable across engines.

**Plain language.** The reference is the "typical patient" the effect is measured around (a 70 kg adult). Changing it only shifts what baseline CL means; the *strength* of the effect ($\beta$) is unchanged. We pin BW = 70, CrCL = 95 so numbers line up one-to-one with NONMEM/PsN.

### 4.4 Selection objective and write-back

**Technical.** Per iteration VAE regresses each $\eta$-dimension's running mean on the candidates and scores subsets with $\mathcal J(S)=\mathrm{RSS}/\omega+\lambda|S|$ (branch-and-bound / `L0Learn`). The winner's **family** is fixed; the `shapes=` choice decides **write-back** as a new coefficient theta named (7.0.2) **`beta.<param>.<cov>.<shape>`** (e.g. `beta.lTVCL.CRCL.power`), categorical `beta.<param>.<cov>`; the structural theta is adjusted jointly so predictions are unchanged. Dropped candidates → `0`. Selected coefficients surface in **`fit$parFixedDf` / `fit$parFixed`** (rownames = coef, column `Estimate`) and `fit$ui$iniDf`.

**Plain language.** Each step, VAE fits a quick penalized regression of "patient CL offsets" on the covariates and keeps the best fit-vs-simplicity subset. Winners get injected back as new model terms; losers get coefficient zero. The final model is read from the fit's parameter table.

## 5. How this maps to the SCM / PsN benchmark (scenario 16)

**Technical.** To search the **identical candidate set** as `runSCM`/PsN, `vaeControl()` is configured (see `script/_smoke_vae_covsel_scn16.R`):

``` r
vaeControl(
  covariateSelection = TRUE,
  shapes = list(BW = c("power","lin"), CrCL = c("power","lin"),
                BMI = c("power","lin"), SEX = TRUE, RACE = TRUE,
                fixCov = TRUE),               # restrict to exactly these 5 covariates
  covCenter = c(BW = 70, CrCL = 95),          # PsN physiological anchors
  covCenterType = "median")                   # BMI -> median (SCM/PsN default)
```

Yields **16 candidates**: `cl,vc × {BW,CrCL,BMI}` in `power`+`lin` (12) plus `cl,vc × {SEX,RACE}` categorical (4) — matching PsN `continuous_covariates=BW,CRCL,BMI` / `categorical_covariates=SEX,RACE`, `[valid_states] continuous=1,4,5 / categorical=1,2`. True scn16 effects: `cl~BW (power)`, `cl~CrCL (power)`, `vc~BW (power)`, `vc~SEX (cat)`; `BMI`, `RACE`, any `lin` pick are distractors → false positives.

**Plain language.** VAE gets the exact same covariate menu and reference weights the NONMEM/PsN search used, so any difference is due to *method*, not setup. We then score how many of the four true effects it recovers and how many spurious ones it adds.

## 6. Glossary

| term | meaning |
|-----------------------------|-------------------------------------------|
| **BSV /** $\eta$ | between-subject variability; each patient's offset from the typical parameter |
| $\Omega$ | variance–covariance matrix of the $\eta$'s |
| **ELBO** | evidence lower bound; VAE training objective (fit − KL) |
| **encoder /** $q_\phi$ | LSTM guessing a patient's $\eta$ from data + covariates |
| **decoder** | the `rxode2` PK model (no neural net, no torch) |
| **KL term** | penalty keeping guessed $\eta$'s close to $\mathcal N(0,\Omega)$ |
| **reparameterization trick** | $\eta=\mu+\Sigma^{1/2}z$ so gradients flow through sampling |
| **shape** | functional form of a covariate effect (power, lin, hockey, cat, …) |
| **center / reference** | value a covariate effect is measured around (BW = 70) |
| **L0 / BIC penalty** | sparsity term $\lambda|S|$ discouraging extra covariates |
| **branch-and-bound** | exact search over covariate subsets in the selection step |
| **`beta.<param>.<cov>.<shape>`** | 7.0.2 name of a selected covariate coefficient |

**Primary references.** Rohleff et al. (2025), *VAEs for population PK modeling*; Kingma & Welling (2014), *Auto-Encoding Variational Bayes*; nlmixr2est `NEWS.md` 7.0.1 → 7.0.2 (shapes / `covCenter` / `fixCov`).

------------------------------------------------------------------------

## Grid

| axis                    | values                    | count |
|-------------------------|---------------------------|-------|
| cohort `N`              | 40, 80, 300               | 3     |
| scenario                | 1 … 16                    | 16    |
| structure               | linCmt, ode               | 2     |
| datasets / cell (pilot) | 1 … `NDS` (default **5**) | 5     |

-   One **LSF array per `(N, scenario, structure)` cell**; the array index (`$LSB_JOBINDEX`) **is the dataset id**.
-   Full pilot = 3 × 16 × 2 = **96 arrays** × 5 datasets = **480 tasks**.

## Files

### Submission (this directory)

| file | role |
|------------------------------------|------------------------------------|
| `vae_covsel_array.lsf` | LSF array template; reads `SAMPLE_N/SCN/STRUCTURE` from env, calls `$SCRIPTS_DIR/vae_covsel_driver.R --dataset $LSB_JOBINDEX` (`SCRIPTS_DIR` auto-detected by the submitter, e.g. `script/hpce_vae_covsel`). |
| `submit_one_array.sh` | Submit ONE cell: `bsub -J "vaecov_N<N>_scn<SS>_<struct>[1-NDS]%MAXPAR"`. |
| `submit_all_arrays.sh` | Triple loop over `NS × SCENARIOS × STRUCTURES`. |

### R scripts (under `script/`, required at runtime)

| file | role |
|------------------------------------|------------------------------------|
| `vae_covsel_driver.R` | **Entry point** run by each array task. Loads the cell's dataset, fits VAE with `covariateSelection = TRUE`, derives the per-scenario `true_set`, and writes the schema-2.1 sidecars. |
| `refit_helpers.R` | `to_nm_dataset()`, `PsN_scenarios` (scenario→covariate indicators). |
| `true_model_factory.R` | Loaded for parity with the smoke file (model-attr helpers). |
| `scm_bench_helpers.R` | Base models `base_2cmt_oral_linCmt` / `base_2cmt_oral_ode`. |
| `output_schema.R` | `assemble_common()` (schema 2.1) and `write_fit_sidecar()`. |
| `estimator_factory.R` | `seed_for_dataset()`, `seed_all()` for reproducible seeding. |

The driver `source()`s the five helpers automatically via a location-robust resolver (`.source_helper` searches `.script_dir`, its parent, and `script/`), so it works whether it lives in `script/` or `script/hpce_vae_covsel/`, and only `vae_covsel_driver.R` is named on the command line.

### R package prerequisites

`nlmixr2 >= 3.0`, `nlmixr2est >= 7.0.2` (provides `est = "vae"` / `vaeControl`
with `shapes=` / `covCenter=` / `fixCov`, and dot-separated coefficient names),
`rxode2` (dev build exporting `getIndCmt`), `nlmixr2scm`, `dplyr`, `tibble`,
`tidyr`.

**No `torch` / libtorch.** `est = "vae"` is a native C++/Armadillo engine. The only optional extra is **`L0Learn`** (CRAN), used solely for large covariate searches (`covSelectMethod = "l0learn"`/`"auto"`). The four-covariate grid here uses the exact branch-and-bound, so `L0Learn` is not required — but it is cheap to install and future-proofs wider grids.

## Output layout (schema 2.1)

```         
output/vae_covsel_pilot/N<N>/scn<SS>_<structure>/covsel/
    res_ds<DDD>.rds        # full record (rel_err, diag, cov, $covsel block)
    res_ds<DDD>.fit.rds    # raw nlmixr2 fit
    res_ds<DDD>.meta.json  # greppable scalar manifest
    res_ds<DDD>_ERROR.txt  # only if the task failed (missing input, fit error)
```

The `$covsel` block carries `selected` (VAE's promoted `beta.<param>.<cov>.<shape>`
terms mapped to `var/covar/shape`) and `true_set` (the scenario's true covariates
**with their true shape** — continuous → `power`, categorical → `cat` — derived
from the `PsN_scenarios` indicators) — the basis for per-cell TP/FN/FP once
aggregated. **Scoring is shape-aware** (see *Selection scoring* below): a true
pair recovered on the wrong shape is *not* a free hit.

## Engine prerequisite (no torch)

`est = "vae"` is a **native C++/Armadillo** LSTM encoder with an analytic backward pass; the decoder is your ordinary `rxode2` model. There is **no `torch` / libtorch dependency** — nothing to download on the login node.

What each node does need is a working **`nlmixr2est >= 7.0.2`** install (which links against `rxode2`, `RcppParallel`/`tbb`, `stringfish`). If a node's build is broken you will see a `LoadLibrary`/`getIndCmt` error at `library()` time, not a torch error. Optionally install `L0Learn` for large covariate searches:

``` r
install.packages("L0Learn")   # optional; only for covSelectMethod="l0learn"/"auto"
```

## Recommended workflow: probe first, then full sweep

VAE resource use on the cluster is unknown, so submit **one array** and read its actual usage before launching all 96.

### 1. Single probe (known-good cell: N80 × scn16 × linCmt × 5 ds)

``` bash
cd ~/nlmixr2scm && git pull
# normalize line endings + clear any stale ERROR sidecars first
sed -i 's/\r$//' script/hpce_vae_covsel/*.sh script/hpce_vae_covsel/*.lsf script/*.R
rm -f output/vae_covsel_pilot/N80/scn16_linCmt/covsel/*_ERROR.txt #Remove them from the probe cell:
find output/vae_covsel_pilot -name '*_ERROR.txt' -delete # clear stale ERROR files across the whole pilot tree

NS=80 SCENARIOS=16 STRUCTURES=linCmt \
  bash script/hpce_vae_covsel/submit_all_arrays.sh 5 20
# equivalently:
bash script/hpce_vae_covsel/submit_one_array.sh 80 16 linCmt 5 20



FULL sweep: all cohorts x all 16 scenarios x both structures,
# datasets 6..250 (245 per cell), 40 concurrent tasks per array
bash script/hpce_vae_covsel/submit_all_arrays.sh 245 40 6
#submit_all_arrays.sh takes [n_ds] [maxpar] [ds_start], and DS_END = DS_START + NDS - 1, so ds_start=6, n_ds=245 gives indices 6–250.

NS="40" bash script/hpce_vae_covsel/submit_all_arrays.sh 245 40 6   
#Submitted 32 array(s), 7840 task(s) total. 2 model parameterization *16 scenarios- each scenario have 245 tasks. 
#The testing queue at novartis has 25 hosts/nodes; each node has 64-core. ~1000 jobs = 20-host limits 
#Occupy at most 20 nodes at once. With mostly 64-core nodes → ~1,280 cores → matches your ~1,066 running tasks. 
#So you're already near that ceiling; the rest pend until tasks finish.

bgadd -L 300 /liuya8j/vaecov 2>/dev/null   # # never run more than 150 tasks at once
NS="80" bash script/hpce_vae_covsel/submit_all_arrays.sh 245 10 6  
NS="300" bash script/hpce_vae_covsel/submit_all_arrays.sh 245 10 6


NS="40 80 300" SCENARIOS=16 STRUCTURES="linCmt ode" \
OUT_ROOT=output/vae_covsel_full0723_est710 \
  bash script/hpce_vae_covsel/submit_all_arrays.sh 245 40 1

  OUT_ROOT=output/vae_covsel_full0807_est703 bash script/hpce_vae_covsel/submit_all_arrays.sh 100 1 1
```

Monitor and inspect resource usage:

``` bash
bjobs -J 'vaecov_*'                 # running/pending
bjobs -l <jobid>                    # live per-task detail
bacct -l <jobid>                    # peak MEM / CPU / walltime AFTER completion
ls output/vae_covsel_pilot/N80/scn16_linCmt/covsel/   # res_ds00{1..5}.*
```

Check for task failures:

``` bash
find output/vae_covsel_pilot -name '*_ERROR.txt' -exec sed -n '1,6p' {} +
```

### 2. Tune resources

From `bacct -l`, read **MAX MEM** and **CPU TIME**; adjust `#BSUB -M` and `#BSUB -W` in `vae_covsel_array.lsf` if the 6000 MB / 3 h probe values are off.

### 3. Full pilot

``` bash
bash script/hpce_vae_covsel/submit_all_arrays.sh 5 20
bjobs -J 'vaecov_*'
```

## Aggregation

Once (some or all of) the sweep has landed on disk, roll the per-dataset `res_ds*.rds` records up into operating characteristics with `script/hpce_vae_covsel/aggregate_vae_covsel.R`. It scans `output/<sub>/N<N>/scn<SS>_<structure>/covsel/res_ds*.rds` and writes 11 CSVs + 1 bundled `.rds` to `output/vae_covsel_aggregated/`.

**Key argument — `--sub` / `sub`** selects which run tree to scan. The full sweep lives in `output/vae_covsel_full/`, so pass `--sub vae_covsel_full` (the default is `vae_covsel_pilot`). `root` stays `output` — it is the common parent of both the input tree (`output/<sub>`) and the output tree (`output/vae_covsel_aggregated`).

### Approach A — command line (HPCE-friendly)

``` bash
Rscript script/hpce_vae_covsel/aggregate_vae_covsel.R --root output --sub vae_covsel_full
# optional custom output dir:
# Rscript script/hpce_vae_covsel/aggregate_vae_covsel.R --root output --sub vae_covsel_full --out_dir output/vae_covsel_full_aggregated
```

### Approach B — interactive R

``` r
source("script/hpce_vae_covsel/aggregate_vae_covsel.R")
res <- aggregate_vae_covsel_run(root = "output", sub = "vae_covsel_full")

res$power                                                 # Power / PowerCN / PowerMinSuc per cell
res$diag_rates                                            # %converged, CN, runtime
subset(res$estim_all, param_class == "covariate_beta")    # clean covariate-beta rel-err
res$covsel_by_cov                                         # per-covariate detection rate
res$relpower                                              # per-k fraction_at_least_k
```

Both approaches produce identical files in `output/vae_covsel_aggregated/`:

| file | contents |
|---------------------------|---------------------------------------------|
| `vae_file_index.csv` | one row per discovered RDS (+ `has_error` flag) |
| `vae_diag_long.csv` | per fit: convergence, CN, timing, selection tally |
| `vae_rse_long.csv` | per (fit, parameter): `rel_err` + `param_class` |
| `vae_covsel_long.csv` | per (fit, var, covar, shape): `in_true`/`in_vae`/`verdict` (TP/FN/FP) |
| `vae_diag_rates.csv` | per cell: %Converged(Strict), CN, MedObjF, runtime |
| `vae_estim_all.csv` | per (cell, parameter): MedRE / MARE / RMRSE (all fits) |
| `vae_estim_success.csv` | same, strict-converged fits only |
| `vae_estim_cond.csv` | same, exact-match fits only |
| `vae_power.csv` | per cell: Power / PowerCN / PowerMinSuc |
| `vae_relpower.csv` | per (cell, k): fraction recovering ≥ k true covariates |
| `vae_covsel_by_covar.csv` | per (cell, var, covar, shape): detection rate + TP/FP/FN |
It aggregates **whatever is on disk at call time** — partially-complete cells just contribute fewer datasets (smaller `N` in `vae_power`/`vae_diag_rates`), so re-run it after the sweep finishes for final numbers.

**Reference caveat** (baked into `param_class`): the DGP centres continuous covariates at fixed references (`BW/70`, `CrCL/95`), while both VAE and runSCM median-centre. Covariate-β coefficients are centring-invariant, so `param_class == "covariate_beta"` rel-err is the clean, comparable estimation metric. `structural_intercept` (TVCL/TVVc) rel-err is reference-dependent — compare with care.

## Notes

-   **Selection scoring is shape-aware (parity with `runSCM` / PsN).** Because 7.0.2 searches a **multi-shape** candidate space (`power` + `lin` per continuous covariate), functional shape is a real degree of freedom, so TP/FN/FP are scored on **`(var, covar, shape)`**, not just `(var, covar)`. A true `power` effect that VAE recovers as `lin` counts as **FN** (true term missed) **+ FP** (spurious `lin` term) — exactly as PsN's strict-power criterion. Shape is normalised (`hockey*` → `hockey`; blank → `NA`), and an `NA` shape on either side acts as a **wildcard**, so legacy records that never stored a shape fall back to `(var, covar)` matching instead of scoring all-FN. `Power = proportion of datasets with exact match` therefore requires the **right pair on the right shape**, matching the SCM operating-characteristic definition. The per-covariate rollup (`vae_covsel_by_covar.csv`) and the covsel heatmap now group on `(var, covar, shape)`, structurally identical to the SCM outputs.
-   **`assemble_common` is untouched** — the VAE-specific `beta.<param>.<cov>.<shape>` → true-label mapping and `rel_err` backfill live only in the driver, so the FOCEi HPCE pipeline is unaffected.
-   **Correct `detection_rate` denominator (fixed at source).** In `vae_covsel_by_covar.csv`, `n_datasets` is the cell's **true dataset (fit) count**, taken from `vae_diag_long.csv` (one row per fit) — *not* the number of datasets in which a pair happened to appear. A distractor only appears when it is falsely selected, so an appearance count would make `n_FP / n_datasets ≈ 1` for every distractor; using the fit count makes `detection_rate = n_detected / n_datasets` the correct per-relationship **sensitivity** (true terms) and **false-positive rate** (distractors). The fit count also includes fits that selected nothing or failed, matching `vae_diag_rates.csv$n_total`. The covsel heatmap therefore no longer needs a `vae_diag_rates.csv` work-around — it divides by `n_datasets` straight from the CSV. (The SCM pipeline carries the identical fix.)
-   The **CRCL uppercase alias** (`CRCL = CrCL`) is added in the driver because VAE uppercases covariate names when building `beta.lTVCL.CRCL.power`.
-   **Pinned search space (7.0.2).** The driver's `vaeControl()` fixes the search to the SCM/PsN candidate set — `cl,vc × BW,CrCL,BMI` in `{power, lin}` + `SEX,RACE` categorical (16 candidates), anchored at `covCenter = c(BW = 70, CrCL = 95)` with BMI at the data median. This is **scenario-independent** (the same candidate menu for all 16 scenarios); only `true_set` varies per scenario. It mirrors `script/_smoke_vae_covsel_scn16.R`.

## ⚠️ Known regression: `nlmixr2est ≥ 7.0.0` and `fit$theta`

**Symptom.** After a manager-side `nlmixr2est` update (≥ 7.0.0), a fresh full re-run scored **`Power = 0` / all-FN** for every cell — even scenarios where the raw fit clearly promoted the correct covariates (verified by hand for `N300 × scn16 × ode`).

**Root cause.** The updated `nlmixr2est` makes **`fit$theta` return `NULL`**. The old driver derived `covsel$selected` by grepping `names(fit$theta)` for `^beta_<PARAM>_<COV>`, so it silently produced `selected = NULL` — no covariate was ever counted as selected. A `tryCatch` masked it as "nothing selected" rather than an error, so **no `_ERROR.txt` was written** (a loud failure became a silent one). The promoted coefficients were never lost — they still live in the record's **`parFixed`** table (7.0.2 rownames `beta.lTVCL.BW.power`, `beta.lTVCL.CRCL.power`, `beta.lTVVc.BW.power`, `beta.lTVVc.SEX`; 7.0.1 used underscores `beta_lTVCL_BW` …; column `Estimate`).

**Fixes (both applied).**

1.  **Driver (`vae_covsel_driver.R`) — prevents recurrence.** `selected` is now read via a version-robust accessor `.named_theta_estimates()` that prefers `parFixed` (→ `fit$parFixedDf` → `fit$parFixed` → legacy `fit$theta`), and parsed by a **dot-aware `.parse_beta_theta()`** that matches both the 7.0.2 dot form (`beta.<param>.<cov>.<shape>`) and the 7.0.1 underscore form, capturing a `shape` column. It also emits a **loud `warning()`** when 0 promoted coefficients are extracted for a scenario that *has* true covariates, so a future API change cannot silently zero the results again.

2.  **Aggregation (`aggregate_vae_covsel.R`) — salvages existing runs.** `.ensure_selected()` rebuilds `covsel$selected` from `parFixed` whenever the stored value is `NULL`/empty (via a **dot-aware `.recover_selected_from_parfixed()`** matching both `beta.<param>.<cov>.<shape>` and legacy `beta_<param>_<cov>`), and `.ensure_relerr_backfill()` then fills the NA covariate-β rows of `rel_err` (`CLBW`, `CLcrCL`, `VcBW`, `VcSEX`) from the same recovered estimates — so **power, FP/FN, detection, AND per-β relative error / RMRSE / MARE all come back**. **No re-fit needed** — just re-aggregate the affected tree, e.g.:

    ``` bash
    Rscript script/hpce_vae_covsel/aggregate_vae_covsel.R --root output --sub vae_covsel_full0722
    ```

    Validation on the hand-checked cell: `N300 × scn16 × ode` went from `Power = 0` to **`Power = 0.90`** (27/30 exact), all four true effects at `detection_rate = 1.0`, and covariate-β `rel_err` fully populated (e.g. `CLBW` −17 %, `CLcrCL` +10 %, `VcBW` −27 %, `VcSEX` −26 %).

**Scope.** The `diag` / `diag_t3` blocks and all scalars (convergence, CN, `cov_ok`, `min_suc`, objf, runtime) were **unaffected** by the update, so convergence / CN gating and the estimation-metric denominators never broke. The only two casualties — `selected` and the covariate-β `rel_err` — are both recovered from `parFixed` at aggregation time. **The VAE aggregation is now fully compatible with `nlmixr2est ≥ 7.0.0`; every output (Power/PowerCN/ PowerMinSuc, FP/FN/detection, MedRE/MARE/RMRSE, relpower) populates correctly.** Only the distractor-β rows (`CLBMI`, `VcBMI`, `VcCrCL`, `VcRACE`) remain `NA` — by design, since they have no `true_value` to compare against.

> **Note.** This regression and its recovery are specific to the **VAE** covariate-selection workflow. The FOCEi / `runSCM` benchmark is a separate pipeline and is not addressed here.