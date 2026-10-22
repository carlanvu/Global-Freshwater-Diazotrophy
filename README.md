# Temperature-threshold model selection for NCD drivers

Exhaustive multiple linear regression (all predictor subsets) below and above a series of temperature thresholds.
Models are ranked by AIC. The best model per regime and threshold is then refit, summarized and plotted.

## Repository layout
```
NCDdrivers_TempSelection.R   main script (run from the repo root)
R/calculate_Reg_v2.R         model-fitting helper (RMSE, AdjR2, AIC, n)
data/source_data.csv         input table: download from the paper's supplementary files and save here with this name (not included in this repo)
results/                     created by the script (tables)
figures/                     created by the script (PNG, 600 dpi)
```

## Input data requirements
- One CSV, one row per sample, **no missing values** (complete cases only).
- All columns are numeric. Every column other than the target is used as a candidate predictor.
- Required columns: the target (`CFG$target`, e.g. `non_cyanoMM`) and the temperature column (`CFG$temp_var`, default `Temp`, in °C).
- Do not include ID, date or other non-predictor columns. Remove them first.
- Temperature is also a predictor if it is in the table, so remove it if you do not want it as one.

## How the temperature split works
For each threshold *T*, samples with `Temp < T` form the "Below" regime and `Temp >= T` the "Above" regime.
Each regime is z-scored (mean 0, SD 1, all columns including the target) and fitted separately.
Same-size subsets are compared within a regime.

## Settings (top of the script, `CFG`)
`target`, `temp_var`, `thresholds`, `min_preds`, `chunk_size`, `workers`, `overfit_cut`, `winning_thr`, `zoom_ids`, `run_search`.

## Outputs
- `TempThreshold_<T>.csv`: all models (raw), `_fixed.csv` (with AICc), `_clean.csv` (valid models only)
- `NCD_AIC_BestModels_AllThresholds(_withPvals).csv`: best model per regime and threshold
- `TempSelectionNCDSummary.csv`: long table used for the AIC vs adjusted R² figures
  (`ID, Group, Temp, AdjR2, n, k, RSS, AIC, p`)
- `NCD_RelImp_Threshold<T>.csv/png`, regression coefficient plots, AIC vs R² figures

## Definitions
- **measured**: input columns
- **computed**: RMSE, adjusted R², AIC, n, k, overall F-test p (`lm` on z-scored data)
- **derived**: RSS = RMSE² × n; AICc = AIC + 2K(K+1)/(n−K−1), K = k + 1
- AIC is `stats::AIC()` (Gaussian likelihood; counts k predictors, the intercept and the residual variance, i.e. k + 2 parameters).
- The AICc correction uses K = k + 1 (as in the published analysis). AICc is only used to exclude models with too few degrees of freedom; best models are selected by AIC.
- Models with n ≤ k + 1 get AIC = NA; AICc is NA if n ≤ k + 2.
- Models with adjusted R² > `overfit_cut` (0.99) in either regime are excluded from selection.

## Caveats
- Subsets of all predictors grow as 2^p: runtime and memory scale accordingly. Use `workers` and `chunk_size`.
- With small n relative to k (e.g. n = 17, k = 14) adjusted R² is near 1 and not reliable. Interpret with AICc/validity filters in mind.
- Selection among many models inflates apparent fit. P-values are from the selected model and are not corrected for the search.

## Requirements
R ≥ 4.2; packages: future.apply, tictoc, ggplot2, ggbreak, relaimpo, gridExtra, ggstatsplot.

## Run
```r
source("NCDdrivers_TempSelection.R")
```
