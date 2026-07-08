Table 1. Preservation of pairwise Pearson correlations of covariates in
bootstrap virtual populations, referenced to NHANES 2021-2023.

|Covariate pair                | NHANES (n = 3576) | N = 40 (250 datasets) | N = 80 (250 datasets) | N = 300 (250 datasets) |
|:-----------------------------|:-----------------:|:---------------------:|:---------------------:|:----------------------:|
|BMI–BW                        |       0.88        |      0.88 ± 0.03      |      0.88 ± 0.02      |      0.88 ± 0.01       |
|BMI–CrCL                      |       0.50        |      0.50 ± 0.07      |      0.50 ± 0.05      |      0.49 ± 0.03       |
|BW–CrCL                       |       0.58        |      0.59 ± 0.06      |      0.59 ± 0.04      |      0.58 ± 0.02       |
|BMI–RACE                      |       0.18        |      0.19 ± 0.07      |      0.18 ± 0.05      |      0.18 ± 0.02       |
|BW–RACE                       |       0.23        |      0.23 ± 0.06      |      0.24 ± 0.04      |      0.23 ± 0.02       |
|CrCL–RACE                     |       0.04        |      0.04 ± 0.07      |      0.04 ± 0.06      |      0.04 ± 0.03       |
|BMI–SEX                       |       -0.04       |     -0.05 ± 0.07      |     -0.04 ± 0.05      |      -0.05 ± 0.03      |
|BW–SEX                        |       0.29        |      0.28 ± 0.07      |      0.29 ± 0.05      |      0.29 ± 0.02       |
|CrCL–SEX                      |       0.08        |      0.08 ± 0.07      |      0.08 ± 0.05      |      0.08 ± 0.03       |
|RACE–SEX                      |       0.00        |     -0.03 ± 0.00      |      0.01 ± 0.00      |      0.00 ± 0.03       |
|Pairs within ±0.05            |         —         |         20/20         |         20/20         |         20/20          |
|Max &#124;Δr&#124; vs NHANES  |         —         |         0.036         |         0.011         |         0.005          |
|Mean &#124;Δr&#124; vs NHANES |         —         |         0.010         |         0.004         |         0.002          |

Notes: NHANES column shows single-sample Pearson r. Cohort columns show
pooled r across all subjects +/- SD across 250 bootstrap datasets.
Delta r = |pooled r - NHANES r|. Stratified bootstrap by SEX x RACE with
rejection gate |cor - NHANES| <= 0.05 on the upper triangle; tolerance
scaled by sqrt(300 / N).
