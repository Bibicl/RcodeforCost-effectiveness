# R code for cost-effectiveness analysis

`upload.R` contains R code defining a three-state cohort model, costs and health utilities, deterministic sensitivity analysis, and probabilistic sensitivity analysis distributions.

## Execution requirements

The script expects an externally prepared R environment and CSV inputs in the working directory. It does not include package-loading or full environment setup instructions. Required package versions have not been established.

Before execution, the following inputs and definitions must be supplied:

- Survival parameter objects: `PFS_Pemb`, `PFS_Pacl`, `OS_Pemb`, and `OS_Pacl`; the script accesses columns by numeric position.
- Spline knot objects: `knots1` and `knots2`.
- The cumulative hazard functions `Hllogis` and `Hsurvspline2`, and model APIs used in the script.
- The death-state definition `state_D`.
- Helper functions `Low25`, `High25`, `Low10`, `High10`, `get_param_value`, `get_normal_sd`, `get_gamma_sd`, `get_alpha`, and `get_beta`.

These external inputs and setup definitions are not included in this upload. The analysis horizon must also be confirmed: the `run_model()` call does not explicitly specify `cycles`.

## Script behavior

The script clears the global R environment, reads CSV files from the working directory, runs the base-case model and deterministic sensitivity analysis, and writes `res_dsa.RData`. The final section defines PSA distributions; it does not call `run_psa()`.

## Validation status

End-to-end execution and reproduction of manuscript results have not been verified. Supplying this script alone does not establish reproducibility of the reported results.
