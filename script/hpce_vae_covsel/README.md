# HPCE VAE covariate-selection pilot

One-run VAE covariate selection (`nlmixr2(..., est = "vae")`,
`covariateSelection = TRUE`) across the full pilot grid on the LSF cluster.
**`runSCM` is not used** — the stepwise search is internal to the VAE training
loop (BICc-ELBO penalty).

## Grid

| axis | values | count |
|------|--------|-------|
| cohort `N` | 40, 80, 300 | 3 |
| scenario | 1 … 16 | 16 |
| structure | linCmt, ode | 2 |
| datasets / cell (pilot) | 1 … `NDS` (default **5**) | 5 |

- One **LSF array per `(N, scenario, structure)` cell**; the array index
  (`$LSB_JOBINDEX`) **is the dataset id**.
- Full pilot = 3 × 16 × 2 = **96 arrays** × 5 datasets = **480 tasks**.

## Files

| file | role |
|------|------|
| `vae_covsel_array.lsf`   | LSF array template; reads `SAMPLE_N/SCN/STRUCTURE` from env, calls `script/vae_covsel_driver.R --dataset $LSB_JOBINDEX`. |
| `submit_one_array.sh`    | Submit ONE cell: `bsub -J "vaecov_N<N>_scn<SS>_<struct>[1-NDS]%MAXPAR"`. |
| `submit_all_arrays.sh`   | Triple loop over `NS × SCENARIOS × STRUCTURES`. |
| `../vae_covsel_driver.R` | The R driver (generalized from `_smoke_vae_covsel_scn16.R`). |

## Output layout (schema 2.1)

```
output/vae_covsel_pilot/N<N>/scn<SS>_<structure>/covsel/
    res_ds<DDD>.rds        # full record (rel_err, diag, cov, $covsel block)
    res_ds<DDD>.fit.rds    # raw nlmixr2 fit
    res_ds<DDD>.meta.json  # greppable scalar manifest
    res_ds<DDD>_ERROR.txt  # only if the task failed (e.g. torch missing)
```

The `$covsel` block carries `selected` (VAE's promoted `beta_*` terms mapped to
`var/covar`) and `true_set` (the scenario's true covariates, derived from the
`PsN_scenarios` indicators) — the basis for per-cell TP/FN/FP once aggregated.

## ⚠️ Torch prerequisite

VAE needs the **`torch`** R package **and** a working **libtorch** backend.
The driver runs a preflight (`requireNamespace("torch")` +
`torch::torch_is_installed()`); on failure it writes `res_ds<DDD>_ERROR.txt`
and exits non-zero (loud fail, not a silent hang).

If the probe returns torch errors, install once on an **online login node**
(libtorch is downloaded — it will fail on an offline compute node):

```r
# on a login node with internet
install.packages("torch")          # or the repo's pinned source
torch::install_torch()             # downloads libtorch backend
torch::torch_is_installed()        # must be TRUE
```

## Recommended workflow: probe first, then full sweep

VAE resource use on the cluster is unknown, so submit **one array** and read
its actual usage before launching all 96.

### 1. Single probe (known-good cell: N80 × scn16 × linCmt × 5 ds)

```bash
cd ~/nlmixr2scm && git pull
NS=80 SCENARIOS=16 STRUCTURES=linCmt \
  bash script/hpce_vae_covsel/submit_all_arrays.sh 5 20
# equivalently:
# bash script/hpce_vae_covsel/submit_one_array.sh 80 16 linCmt 5 20
```

Monitor and inspect resource usage:

```bash
bjobs -J 'vaecov_*'                 # running/pending
bjobs -l <jobid>                    # live per-task detail
bacct -l <jobid>                    # peak MEM / CPU / walltime AFTER completion
ls output/vae_covsel_pilot/N80/scn16_linCmt/covsel/   # res_ds00{1..5}.*
```

Check for torch failures:

```bash
find output/vae_covsel_pilot -name '*_ERROR.txt' -exec sed -n '1,6p' {} +
```

### 2. Tune resources

From `bacct -l`, read **MAX MEM** and **CPU TIME**; adjust `#BSUB -M` and
`#BSUB -W` in `vae_covsel_array.lsf` if the 6000 MB / 3 h probe values are off.

### 3. Full pilot

```bash
bash script/hpce_vae_covsel/submit_all_arrays.sh 5 20
bjobs -J 'vaecov_*'
```

## Notes

- **`assemble_common` is untouched** — the VAE-specific `beta_*` → true-label
  mapping and `rel_err` backfill live only in the driver, so the FOCEi HPCE
  pipeline is unaffected.
- The **CRCL uppercase alias** (`CRCL = CrCL`) is added in the driver because
  VAE uppercases covariate names when building `beta_lTVCL_CRCL`.
- Aggregation (`aggregate_vae_covsel.R`) is intentionally deferred until the
  first pilot batch returns and the record shape is confirmed against real
  cluster output.
