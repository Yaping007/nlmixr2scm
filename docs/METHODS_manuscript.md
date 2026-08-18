# Methods

## 1. Overview and study aim

We conducted a simulation study to evaluate automated covariate-model-building for population pharmacokinetic (PopPK) analysis in the open-source **nlmixr2** ecosystem, benchmarked against the **NONMEM/PsN** reference standard. The study is organised around two comparison axes:

1.  **Open-source FOCE-I vs. NONMEM.** Stepwise covariate modelling (SCM) implemented natively in `nlmixr2` (the `nlmixr2scm::runSCM` tool, FOCE-I with the BOBYQA outer optimiser) is compared against **PsN SCM** on NONMEM (Perl-speaks-NONMEM), the established reference implementation. NONMEM/PsN serves as the reference benchmark against which the open-source FOCE-I workflow is calibrated.
2.  **FOCE-I vs. VAE within nlmixr2.** The stepwise FOCE-I search is compared against a **variational-autoencoder (VAE)** estimator (`nlmixr2est`, `est = "vae"`) that performs simultaneous parameter estimation and covariate selection in a single fit, without a forward/backward search.

The design replicates and extends the simulation framework of Khandelwal et al. (2019), in which a known covariate structure is imposed on a virtual population, data are simulated across a factorial grid of covariate-effect scenarios and sample sizes, and each workflow is evaluated on its ability to recover the true covariate model. Workflows were compared on five families of operating characteristics: (i) numerical convergence, (ii) statistical power to detect true covariates, (iii) covariate-selection pattern (per-covariate true-positive and false-positive behaviour), (iv) parameter-estimation accuracy, and (v) run-time efficiency.

## 2. Virtual population and covariate distributions

A virtual population of 100 replicate datasets of up to 300 subjects each was generated (`virtual_population_250x300`). Covariate distributions and their joint correlation structure were derived from the **NHANES** adult survey so that the simulated subjects reflect a realistic adult patient population. For every subject, five candidate covariates were sampled:

| Covariate                   | Type        | Role                            |
|-------------------------|------------------|-----------------------------|
| Body weight (BW)            | continuous  | true covariate on CL and Vc     |
| Creatinine clearance (CrCL) | continuous  | true covariate on CL            |
| Body mass index (BMI)       | continuous  | distractor (never truly active) |
| Sex (SEX)                   | categorical | true covariate on Vc            |
| Race (RACE)                 | categorical | distractor (never truly active) |

BMI and RACE were included as **non-informative distractor covariates**: they are correlated with the true covariates (notably BMI with BW) but never carry a true effect. Their inclusion lets us measure collinearity-driven false-positive selection (e.g. BMI→CL capturing the BW→CL slot). Covariate correlation structure was preserved across replicates and is reported in Table 1 (`output/table1_covariate_correlation`).

Three cohort sizes were studied — **N = 40, 80 and 300 subjects** — to characterise small-, moderate-, and large-sample behaviour.

## 3. Structural pharmacokinetic model

Concentration–time data were simulated from a **two-compartment model with first-order oral absorption** following a single 100 mg dose. The model was parameterised in terms of clearance (CL), central volume (Vc), inter-compartmental clearance (Q), peripheral volume (Vp) and first-order absorption rate (KA). True population values were:

| Parameter | True value           |
|-----------|----------------------|
| CL        | 0.6 (L/h)            |
| Q         | 1.8 (L/h)            |
| Vc        | 20 (L)               |
| Vp        | 80 (L)               |
| KA        | 0.7 (1/h), **fixed** |

Between-subject variability (BSV) was placed on CL and Vc as a correlated bivariate log-normal random effect (variances 0.1 and 0.1, covariance 0.02). Residual variability was proportional (0.1). KA was fixed at its true value because the sparse sampling design does not identify absorption. Observations followed a **sparse 6-sample-per-subject** schedule, consistent with the information-poor conditions under which covariate selection is most challenging.

The structural model was implemented in two mathematically equivalent forms: an **analytic solution** (`linCmt()`) used for the primary analysis, and an **explicit ODE** form used to enable analytic-gradient estimators (see §6). Results are reported separately for the `linCmt` and `ode` structural implementations.

## 4. Covariate model and the 16-scenario factorial

The true covariate model followed the log-additive form of Khandelwal et al. (2019). Continuous covariates entered CL and Vc as power (allometric-style) terms centred at population reference values, and sex entered Vc as a log-linear categorical shift:

$$
\mathrm{CL} = \exp\!\big(\theta_{CL} + \theta_{BW,CL}\log(BW/70) + \theta_{CrCL,CL}\log(CrCL/95) + \eta_{CL}\big)
$$

$$
\mathrm{Vc} = \exp\!\big(\theta_{Vc} + \theta_{BW,Vc}\log(BW/70) + \theta_{SEX,Vc}\,\mathrm{SEX} + \eta_{Vc}\big)
$$

The four true covariate–parameter relationships and their coefficients were:

| Relationship          | True coefficient |
|-----------------------|------------------|
| BW → CL (power)       | 0.75             |
| CrCL → CL (power)     | 0.50             |
| BW → Vc (power)       | 1.00             |
| SEX → Vc (log-linear) | log(1.5) ≈ 0.405 |

Each of the four effects was independently switched **on or off**, giving a **2⁴ = 16-scenario full-factorial grid**. Scenario 1 is the **null model** (no true covariate effects; used to measure the family-wise false-positive rate), and scenario 16 contains **all four true effects** simultaneously. The full indicator grid is stored in `PsN_scenarios` and the corresponding true parameter table in `true_params_long`.

The complete simulation grid was therefore **3 sample sizes × 16 scenarios × 250 replicate datasets** per workflow and structural form.

## 5. Reference "true-model" refit

To establish an accuracy ceiling independent of model selection, the **true data-generating model** (with the correct covariates already in place) was refitted to every simulated dataset. This "given the truth, can the estimator recover it?" analysis was performed with **FOCE-I + BOBYQA**. It provides the reference bias and relative root-mean-square error (RMRSE) against which the covariate-search workflows are compared.

## 6. Covariate-selection workflows

### 6.1 nlmixr2 SCM (`runSCM`)

Automated SCM was performed with `nlmixr2scm::runSCM`, a from-scratch `nlmixr2` re-implementation of the PsN SCM algorithm (Jonsson & Karlsson, 1998; Lindbom et al., 2004). Starting from a covariate-free two-compartment base model, the search proceeded in two phases:

-   **Forward inclusion** — at each step, every remaining candidate covariate–parameter–shape relationship was added in turn; the relationship producing the largest likelihood-ratio-test (LRT) improvement significant at **α = 0.05** (ΔOFV \> 3.84, χ²₁) was retained. Steps continued until no candidate reached significance.
-   **Backward elimination** — starting from the forward-final model, each included relationship was tested for removal at the stricter threshold **α = 0.01** (ΔOFV \> 6.63). The backward-eliminated model is the **final selected model**.

The candidate space comprised the two disposition parameters **{CL, Vc}** crossed with three continuous covariates **{BW, CrCL, BMI}** under **power** and **linear** shapes, plus two categorical covariates **{SEX, RACE}** — the scenario-16 superset applied to all scenarios so that both true effects and distractors were always eligible. `runSCM` builds the covariate terms internally with user-specified shape and reference for continuous covariates and generating categorical indicators, so no manual model editing was required.

For efficiency, each step used a **two-tier** fitting scheme. Every candidate was first evaluated with a fast **screening fit** that skips the covariance step and output-table construction (`covMethod = ""`, `calcTables = FALSE`), since only the objective-function value is needed for the likelihood-ratio test. The two tiers share identical estimation settings (the same sigdig, derivative step, and ODE tolerances); once the winning candidate is chosen, it is refitted in a **final tier** that adds the covariance step and diagnostic tables (`covMethod = "r,s"`, `calcTables = TRUE`) for reporting. The screening tier is therefore purely a run-time optimisation, not a change in estimation accuracy. Candidate covariate coefficients were started at the **nlmixr2 default** initial estimate. On sparse, information-poor data a gradient outer optimiser can stall at this flat zero-effect starting point and reject a true covariate. To guard against this, a candidate whose nested-model objective-function value fails to improve on its parent (ΔOFV ≤ 0, which cannot happen at a genuine optimum and therefore signals the optimiser never left the init) triggers a **1-D frozen-base profile** for that single coefficient: all structural θ are fixed at the parent estimates and the between-subject variability is held fixed but present, so the covariate slope is profiled under the correct mixed-effects likelihood. The profiled value — obtained by a cheap bounded scalar line search in the spirit of Brent's method (a boundary solution is rejected and the default init retained) — then warm-starts a rescue refit, supplying gradient optimisers with a nonzero, gradient-informative starting value.

### 6.2 PsN SCM (NONMEM reference benchmark)

The identical search specification (same base model, candidate space, and forward/backward significance levels) was executed with **PsN SCM** on NONMEM. NONMEM/PsN is treated as the **reference benchmark** against which the open-source `nlmixr2` FOCE-I workflow is compared (`output/psn_scm`, `output/psn_scm_full0727`). The PsN reference is currently available for the analytic (`linCmt`) structural form; the corresponding NONMEM ODE benchmark is pending.

### 6.3 VAE simultaneous estimation and selection

The VAE workflow used the native `nlmixr2est` variational-autoencoder estimator (`est = "vae"`), a C++/Armadillo LSTM-encoder implementation with an rxode2 decoder. Unlike SCM, VAE performs **parameter estimation and covariate selection jointly in a single fit** (`covariateSelection = TRUE`); no forward/backward search is run. For the four-covariate problem the exact branch-and-bound selection branch was used (L0-regularised selection is available for larger searches). VAE therefore constitutes a fundamentally different, search-free alternative to stepwise selection.

## 7. Diagnostics and operating characteristics

### 7.1 Convergence (PMx-strict)

A single, harmonised convergence rule was applied across all workflows. Because gradient optimisers report code-8 "false convergence" on good fits, convergence was defined by **numerical evidence** rather than the raw optimiser code:

```         
converged        = is.finite(OFV) & covariance-step succeeded & is.finite(condition_number)
converged_strict = converged & (condition_number ≤ 1000) & (no θ estimated at a boundary)
```

The `condition_number ≤ 1000` and no-boundary criteria follow the pharmacometrics-strict definition of Khandelwal et al. (2019, Table 3). `NA` diagnostics were treated as failures (`%in% TRUE` idiom) to make aggregation NA-safe.

### 7.2 Power and covariate-selection pattern

For each (workflow, scenario, N) cell, **power** was defined as the proportion of replicates recovering the exact true covariate set. Ancillary definitions included power conditional on acceptable conditioning (`PowerCN`) and on minimal successful selection (`PowerMinSuc`). In the null scenario, power reduces to 1 − family-wise false-positive rate. The **selection pattern** was summarised per candidate covariate as an *error* rate against the truth: **false-negative** rate for the four true relationships (a true effect not recovered, or recovered in the wrong shape) and **false-positive** rate for the BMI and RACE distractors (a null covariate wrongly selected). These per-covariate error rates were displayed as covariate × scenario heatmaps.

### 7.3 Estimation accuracy

Parameter-estimation accuracy was quantified against the known true values using **relative estimation error**, its median (bias), **relative root-mean-square error (RMRSE)** and **mean absolute relative error (MARE)**. Following Khandelwal et al. (2019), the **primary** metrics are **unconditional**: they are computed over **all fits** with a finite estimate, without restricting to converged runs, so that the reported accuracy reflects the workflow as actually deployed. For structural parameters this uses every dataset; for a covariate coefficient it uses every dataset in which that covariate was selected. Two secondary variants are reported as sensitivity analyses: one restricted to **strictly converged** fits, and a **selection-conditional** variant restricted to datasets in which the exact true covariate model was recovered (the δ = 1 indicator of Khandelwal et al. 2019).

### 7.4 Run-time efficiency

Wall-clock time per fit was recorded for every task and summarised per (workflow, scenario, N, structure) as median, mean and 90th-percentile run-time and total compute-hours, on a log scale to span the estimator range.

## 8. Result synthesis and statistical comparison of methods

The five operating-characteristic families are summarised into a single **scorecard** so that the head-to-head message can be read at a glance while the underlying per-replicate evidence still supports every claim. The scorecard is deliberately built in two complementary layers — a **descriptive** layer (what the numbers are) and an **inferential** layer (whether the differences are real) — because these answer different questions and are computed on different scales.

### 8.1 Reference frame and the six scalars

All synthesis is reported at the moderate cohort size **N = 80**, separately for the **Analytic** (`linCmt`/ADVAN4) and **ODE** structural forms, with **PsN-SCM taken as the reference method** so that every comparison reads as "method vs the NONMEM/PsN standard". Each workflow × structure cell is reduced to six scalars, one per operating-characteristic family:

| Dimension | Scalar (per method × structure) | Better |
|-------------|----------------------------------------------|-------------|
| Convergence | median % of replicates minimised | higher |
| Cov. step | median % with a successful covariance step | higher |
| Power | median exact-match power (%) | higher |
| False positive | median distractor false-positive rate (%) | lower |
| MARE | median mean-absolute-relative-error of the 4 covariate βs | lower |
| Runtime | median wall-clock minutes per fit | lower |

The scenario dimension is collapsed by the **median over the 16 scenarios**, chosen for robustness to the wide scenario-to-scenario spread (power alone ranges from \~90% to \~18% across scenarios). **Runtime is compared head-to-head for every method** on a common endpoint — *wall-clock time until covariate selection is complete* — which is the full forward/backward stepwise pipeline for SCM and the single joint fit for VAE; VAE finishing selection in one pass is the property under test, not a reason to exclude it. The only remaining like-for-like caveat is that the **Analytic** stratum pairs `linCmt` (nlmixr2) against ADVAN4 (NONMEM); this is noted in the caption rather than hidden by omission.

### 8.2 The three scorecard views

-   **Layer A — raw table** (`scorecard_table.csv`/`.html`): the six scalars in natural units, six rows (3 methods × 2 structures).
-   **Layer B — heatmap** (`fig_scorecard_heat`): the inferential object. The **number** in each cell is the descriptive median (natural units); the **colour** is that median expressed as a **relative change vs PsN**, direction-corrected and clipped to ±100% so that blue always means "better than PsN", red "worse", and white "equal to PsN"; a **significance marker** (see §8.3) reports whether the method-vs-PsN difference is statistically real *and* practically non-negligible.
-   **Layer C — radar** (`fig_scorecard_radar`): a glanceable **direction** summary only. Each spoke is the relative median vs PsN; the two structures are overlaid, with a dashed mid-ring at the PsN reference. The radar deliberately carries **no significance markers and no single weighted composite** — it answers "what is the shape of the trade-off", not "how much" or "is it significant". Because radar radius cannot faithfully co-display a 5% power gap and a 45× runtime gap on one scale, runtime magnitude is read from the dedicated log-scale run-time figure, not from the radar.

A subtlety that the caption states explicitly: on the heatmap the **cell number (a median) and the significance marker (a model effect) are different estimators of the same "vs PsN" question**. They agree in direction and in the obvious cases but need not agree in value, because the median lives in natural units whereas the model effect lives on a link scale (see §8.4). The number is chosen for readability; the marker for defensibility. For the accuracy (MARE) row this split is sharper still: the displayed number collapses the four covariate coefficients by their **median** (robust to the occasional exploded relative error at a small true $\beta$), whereas the tested response averages them by their **mean** within each replicate. This is a deliberate, documented choice — the same accuracy construct summarised two ways — not an inconsistency; keeping the display as a median-of-medians preserves a directly interpretable "the typical coefficient is recovered to within $x\%$" reading that a per-replicate mean would sacrifice.

### 8.3 Paired, structure-stratified significance testing

The comparison is **paired by construction**: the three workflows analyse the *same* simulated datasets, so each metric is tested on the raw per-replicate records (keyed by `dataset_id`) rather than on the aggregated medians, pairing on the composite key **(scenario, dataset_id)** and restricting to the `dataset_id` intersection shared by the two methods being compared. Testing is **stratified by structure** (fit separately within Analytic and within ODE), because the two structural forms are different experiments — most starkly, PsN collapses on ODE — so a per-cell result maps one-to-one onto a heatmap row. **Multiplicity** across the tested comparisons is controlled with the **Holm** step-down correction. The marker tiers are `*` (p\<0.05), `**` (p\<0.01), `***` (p\<0.001), `ns` (not significant), and `≈` (statistically significant but within the practical-equivalence band of §8.6 — reliable but too small to matter).

### 8.4 Primary model: paired mixed-effects, and the meaning of its coefficient

Although three workflows are compared (PsN-SCM, nlmixr2-SCM, VAE), the `method` factor in each model has **exactly two levels, not three**. PsN-SCM is fixed as the reference, and every non-reference workflow is tested against it in its **own separate two-level model** — one fit for *nlmixr2-SCM vs PsN-SCM* and a second, independent fit for *VAE vs PsN-SCM*. There is no single three-level model carrying two contrasts at once. This pairwise design is deliberate: each comparison is refitted on **only the `dataset_id` intersection shared by that particular pair**, so a method is never penalised for datasets the other method happened not to run, and each contrast keeps its own maximal paired sample rather than being restricted to the datasets common to all three workflows. The reference row (PsN-SCM against itself) is not modelled. Multiplicity across this family of pairwise comparisons is then controlled by the Holm adjustment (§8.3).

The primary test for each metric × structure is a **mixed-effects model** with PsN-SCM as the reference level of this two-level `method` factor (levels `ref` = PsN-SCM, `cmp` = the single workflow under test):

$$
g\big(\mathbb{E}[y_{sdm}]\big) \;=\; \beta_0 \;+\; \beta_{\text{method}}\,\mathbb{1}(m = \text{method}) \;+\; u_s \;+\; v_{sd},
\qquad u_s \sim \mathcal{N}(0,\sigma^2_{\text{scen}}), \; v_{sd} \sim \mathcal{N}(0,\sigma^2_{\text{data}})
$$

for outcome $y_{sdm}$ on dataset $d$ within scenario $s$ under method $m$, where the indicator $\mathbb{1}(m=\text{method})$ selects the single non-reference level (`cmp`), $u_s$ is a **scenario** random intercept — `(1 | scenario)` — and $v_{sd}$ is a **dataset-within-scenario** random intercept — `(1 | scenario:dataset_id)`. The dataset term is what makes the model the multi-scenario generalisation of a paired test: it links the two fits of one dataset, exactly as McNemar/paired-Wilcoxon do for a single stratum. `dataset_id` is **nested**, not crossed (dataset 5 of scenario 1 is a different unit from dataset 5 of scenario 2), and structure is a **fixed stratifier**, not a random grouping (a dataset is never re-used across structures). The link $g$ and the interpretation of $\beta_{\text{method}}$ depend on the metric family:

-   **Binary metrics** (convergence, cov. step, power, false-positive) — `glmer(..., family = binomial)`; $g=\text{logit}$, so $\beta_{\text{method}}$ is a **log-odds difference** vs PsN and $e^{\beta}$ is an **odds ratio**.
-   **Runtime** — `lmer` on $\log(\text{minutes})$; $\beta_{\text{method}}$ is a **log-fold difference**, so $e^{\beta}$ is a **fold-change** (e.g. VAE on ODE $e^{\beta}\approx0.03$, i.e. \~33× faster than PsN).
-   **MARE** — `lmer` on the percentage; $\beta_{\text{method}}$ is a **mean difference** in error (percentage points).

In words, $\beta_{\text{method}}$ is the **shrunk, dataset-paired mean effect of the method versus PsN on the link scale, after adjusting for how hard each scenario is**. It is *not* the same quantity as the radar/heatmap median: the displayed median is the median across scenarios of the per-scenario (method − PsN) difference in natural units, whereas $\beta$ is a mean-like, partially-pooled effect on a log-odds / log-fold scale. They routinely agree in sign; they are not expected to agree in magnitude.

Collapsing scenarios to one number *before* testing is deliberately avoided: it would discard the per-dataset pairing and within-scenario variability that give the test its power. The mixed model performs that collapse **internally and correctly**, yielding one pooled effect while honouring both the pairing and the scenario heterogeneity visible across the x-axis of the per-scenario figures. The same fitted model therefore serves as a single source of truth for both the heatmap markers and any per-panel annotation of the corresponding per-metric line plot (per structure).

### 8.5 Fallback and ceiling-effect handling

Two situations make the Wald inference from the mixed model unreliable, and both fall back to the **exact paired tests that the mixed model generalises**:

1.  When `lme4` is unavailable (e.g. an offline session), every comparison uses **exact McNemar** (binomial test on the discordant pairs) for binary metrics and the **paired Wilcoxon signed-rank** test for continuous metrics. These are also reported as a sensitivity column when the mixed model is used.
2.  **Ceiling / separation.** Several binary metrics sit at 100% for the nlmixr2 workflows (e.g. convergence, cov. step), producing complete or quasi-separation that inflates the `glmer` Wald standard error and drives the p-value spuriously toward 1 despite an obvious raw difference. A guard detects the degenerate fit — a non-finite or implausibly large $|\beta|$, or a non-significant p accompanying a clear raw gap — and substitutes the **exact McNemar** result (labelled *separation fallback*) for that cell, which remains valid on the discordant pairs.

The complete inferential table (effect estimate, 95% CI, effect scale, raw and Holm-adjusted p, test used, number of paired datasets, and marker) is written to `output/figures/scorecard/scorecard_stats.csv`; the implementation is `script/viz/scorecard_stats.R`, and the descriptive scorecard and its figures are produced by `script/viz/scorecard.R`.

### 8.6 Exact definitions of the three cell channels

Each heatmap cell encodes three deliberately-separated quantities. For method $m$, structure $s$ and metric $k$, let $x_{k}(m,s,c,d)$ be the raw outcome on scenario $c\in\{1,\dots,16\}$ and dataset $d\in\{1,\dots,100\}$, let $r=\text{PsN-SCM}$ be the reference, and let $\text{dir}_k=+1$ for higher-is-better metrics (convergence, cov. step, power) and $-1$ for lower-is-better metrics (false-positive, MARE, runtime).

**(A) The number — descriptive median.** Datasets are first collapsed to a per-scenario value $v_k$, then scenarios are collapsed by the median:

$$
N_k(m,s) \;=\; \operatorname{median}_{c=1..16} v_k(m,s,c),
$$

with the per-scenario value being a percentage for the rate metrics, $v = 100\cdot\tfrac{1}{100}\sum_d \mathbb{1}[\text{success}_{cd}]$; a median wall-time for runtime, $v = \operatorname{median}_d(\text{wall}_{cd}/60)$ minutes; and, for accuracy, the **displayed** value is a fully median-based collapse of the four covariate coefficients: each coefficient's absolute relative error is first collapsed across datasets by its median, $\text{MARE}_{cj} = 100\cdot\operatorname{median}_d\lvert(\hat\beta_{jcd}-\beta_j)/\beta_j\rvert$, and the cell number is the median of these over the pooled $16\times 4$ (scenario $\times$ coefficient) values. The median (not the mean) is used throughout for robustness to the wide scenario spread and to the heavy right tail of relative error near small true $\beta$. **Note the display and the significance test collapse the four coefficients differently** (see §8.4): the number takes a *median* over the four $\beta$s for robustness, whereas the paired mixed model of §8.4 uses the per-replicate *mean* of the four $\beta$s, $\text{MARE}_{cd} = \tfrac{100}{4}\sum_{j=1}^{4}\lvert(\hat\beta_{jcd}-\beta_j)/\beta_j\rvert$, as its response. Both summarise the same accuracy construct; the median is chosen for display readability and robustness, the mean for a well-behaved model response.

**(B) The colour — relative change vs the reference, clipped.** The direction-corrected gap and its relative form are

$$
\Delta_k(m,s) = \text{dir}_k\big(N_k(m,s) - N_k(r,s)\big),
\qquad
\rho_k =
\begin{cases}
0 & N_k(r,s)=0,\ \Delta=0\\
\operatorname{sign}(\Delta) & N_k(r,s)=0,\ \Delta\neq0\\
\Delta_k / \lvert N_k(r,s)\rvert & \text{otherwise,}
\end{cases}
$$

and the fill value is clipped to $\pm1$ and mapped through a diverging red–white–blue scale, $\text{fill}_k=\max(\min(\rho_k,1),-1)$, with $-1$ red (worse), $0$ white ($=$ reference), $+1$ blue (better). Because $\text{dir}_k$ has already flipped the sign, **blue is "better than PsN" for every metric**. Relative (rather than absolute or per-column) scaling keeps a fraction-of-a-point MARE gap pale while a tens-of-points power or runtime gap saturates.

**(C) The significance marker — paired test, then a practical gate.** The primary per-cell test is the mixed-effects model of §8.4; its coefficient gives $z=\hat\beta/\operatorname{SE}(\hat\beta)$ and $p=2\,\Phi(-\lvert z\rvert)$. For binary metrics at a ceiling, complete separation is detected by

$$
\lvert\hat\beta\rvert>10 \;\lor\; \neg\,\text{finite}(\hat\beta) \;\lor\; \big(p>0.05 \,\land\, \lvert\bar y_m-\bar y_r\rvert>0.05\big),
$$

in which case the cell falls back to the **exact McNemar** test on the discordant pair counts $b$ (method wins) and $c$ (reference wins), $p = 2\sum_{i\ge\max(b,c)}\binom{b+c}{i}0.5^{\,b+c}$. All tested p-values are Holm-adjusted to $\tilde p$. A **practical-equivalence band** then gates the glyph so that a bare star can never signal a trivial effect: a difference is negligible if it is small on *either* the relative *or* an absolute scale,

$$
\text{negl} \;=\; \big(\lvert\rho_k\rvert<\delta_{\text{rel}}\big)\;\lor\;\big(\text{unit}=\%\ \land\ \lvert\Delta_k\rvert<\delta_{\text{abs}}\big),
\qquad \delta_{\text{rel}}=0.10,\ \delta_{\text{abs}}=2\text{ pp},
$$

and the marker is

$$
\text{marker}=
\begin{cases}
\text{“”} & m=r\ (\text{reference row})\\
\approx & \tilde p<0.05 \,\land\, \text{negl}\\
\text{***} & \tilde p<0.001 \,\land\, \neg\text{negl}\\
\text{**} & \tilde p<0.01 \,\land\, \neg\text{negl}\\
\text{*} & \tilde p<0.05 \,\land\, \neg\text{negl}\\
\text{ns} & \tilde p\ge0.05.
\end{cases}
$$

The dual band handles both regimes of baseline: the **relative** threshold governs large-baseline metrics (e.g. MARE near 12%), while the **absolute floor** protects near-zero baselines — a false-positive rate of 3% vs the reference 2% is a $+50\%$ relative change but only a $1$-pp absolute gap, so it is correctly rendered `≈` (reliable but negligible) rather than `***`. Colour therefore answers *how big* while the marker answers *how certain*, and a star is printed only when a difference clears **both** a statistical and a practical bar. The thresholds $\delta_{\text{rel}}$ and $\delta_{\text{abs}}$ are exposed as arguments (`practical_delta`, `absolute_floor`) so the minimally-important effect is an explicit, auditable choice rather than an implicit consequence of sample size.

## 9. Computation and reproducibility

All fits were executed on an **IBM Spectrum LSF** high-performance compute cluster as job arrays over the (N × scenario × dataset × workflow) grid, with each task pinned to a single node (4 CPU, 6 GB RAM, 3 h wall-time). To avoid fork-over-OpenMP corruption of the ODE solver, `runSCM` candidate parallelism and rxode2 threading were mutually constrained (rxThreads forced to 1 whenever candidate workers \> 1). Every fit wrote a versioned schema-2.1 record (fit object, diagnostics, timing, and a JSON manifest) for reproducible aggregation. Results were aggregated with dedicated R pipelines (`aggregate_*`), applying the PMx-strict convergence rule uniformly, and rendered into the convergence, power, selection-pattern, accuracy and run-time figures reported here.

Analyses were performed in R 4.5.3 with `nlmixr2`, `nlmixr2est`, `rxode2`, and `nlmixr2scm`.

## References

-   Khandelwal A. et al. (2019). *\[covariate model building simulation study\]*.
-   Jonsson E.N. & Karlsson M.O. (1998). Automated covariate model building within NONMEM. *Pharmaceutical Research*, 15(9), 1463–1468.
-   Lindbom L., Ribbing J., & Jonsson E.N. (2004). Perl-speaks-NONMEM (PsN). *Computer Methods and Programs in Biomedicine*, 75, 85–94.
-   Ribbing J. & Jonsson E.N. (2004). Power, selection bias and predictive performance of the population pharmacokinetic covariate model. *JPKPD*, 31(2), 109–134.