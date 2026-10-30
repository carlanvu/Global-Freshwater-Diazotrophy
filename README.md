# Global freshwater diazotrophy: analysis and figure code

Code for the model selection (temperature-threshold multiple linear regressions) and for the figures of the manuscript.

## Repository layout
```
scripts/
  NCDdrivers_TempSelection.R      exhaustive threshold model selection (NCD) -> results/
  fig1_bubbles_and_insets.R       Fig 1B, 1D + insets 1B, 1D, 1E
  fig1F_nifH_cluster_violin.R     Fig 1F
  fig2_taxonomic_composition.R    Fig 2B, 2C-i, 2C-ii, 2C-iii
  fig3_fig4_driver_panels.R       Fig 3 and 4, panels A-C (set FIG at the top)
R/
  calculate_Reg_v2.R              model-fitting helper (RMSE, AdjR2, AIC, n)
  supplement_io.R                 readers for the supplementary tables + variable-name map
  figure_helpers.R                shared colours, palettes, helpers
data/                             supplementary files (not stored in this repo)
results/, figures/                created by the scripts
```
Run every script from the repo root, e.g. `source("scripts/fig1_bubbles_and_insets.R")`.

## Input files (supplementary data; save in `data/`)
| File | Used by |
|---|---|
| `20261019_GFREDIZ23_metadata.xlsx` (full metadata, one row per site) | all scripts |
| `FullTaxCounts_GF23_20250814OutMM.txt` (ASV counts + taxonomy) | Fig 2 |
| `GF23CountsTaxaClusters_clean.csv` (ASV counts + nifH cluster) | Fig 1F |
| `Figure3_metadata.xlsx`, `Figure4_metadata.xlsx` (optional) | pre-computed importance values (`importance_source = "sheet"`) and reproducibility check |

Column names in the supplement sheets are descriptive (e.g. "Phototrophs abundance"); `R/supplement_io.R`
maps them to the column names of the full metadata (`PhotoBA`, ...).

## Figures: script map
| Figure | Script | Notes |
|---|---|---|
| 1A, 1C, 1E (maps), 2A map | – | ArcGIS Pro / Excel, not scripted |
| 1B, 1D, insets 1B/1D/1E | `fig1_bubbles_and_insets.R` | |
| 1F | `fig1F_nifH_cluster_violin.R` | needs ggplot2 >= 3.5.0 |
| 2B, 2C-i, 2C-ii, 2C-iii | `fig2_taxonomic_composition.R` | |
| 3 (N2 fixation) A-C | `NCDdrivers_TempSelection.R` (RUN = "N2fix"), then `fig3_fig4_driver_panels.R`, FIG = "fig3" | |
| 4 (NCD) A-C | `NCDdrivers_TempSelection.R` (RUN = "NCD"), then `fig3_fig4_driver_panels.R`, FIG = "fig4" | |

Final layouts (inset overlays, labels, significance marks) were assembled in Canva from the layers saved in `figures/`.

## Model selection
`scripts/NCDdrivers_TempSelection.R` runs the search for either response: `RUN = "NCD"` (target `non_cyanoMM`, thresholds 8-12 -> `results/fig4/`) or `RUN = "N2fix"` (target `N2FixationRates`, thresholds 15-25 -> `results/fig3/`). Predictor sets and thresholds are in the `PRESET` list at the top.
For each threshold *T*, samples with `Temp < T` form the "Below" regime and `Temp >= T` the "Above" regime.
Each regime is z-scored (all numeric columns) and every predictor subset is fitted by `lm`; models are ranked by AIC.
The best model per regime and threshold is refit (overall F-test p) and summarized in `results/<run>/TempSelectionSummary.csv`
(`ID, Group, Temp, AdjR2, n, k, RSS, AIC, p, Combination`). Missing values are handled per model (complete cases).

### Definitions
- **measured**: input columns
- **computed**: RMSE, adjusted R², AIC, n, k, overall F-test p (`lm` on z-scored data)
- **derived**: RSS = RMSE² × n; AICc = AIC + 2K(K+1)/(n−K−1), K = k + 1
- AIC is `stats::AIC()` (Gaussian likelihood; k predictors + intercept + residual variance = k + 2 parameters).
- The AICc correction uses K = k + 1 (as in the published analysis). AICc is only used to exclude models with too few degrees of freedom; best models are selected by AIC.
- Models with n ≤ k + 1 get AIC = NA; AICc is NA if n ≤ k + 2.
- Models with adjusted R² > `overfit_cut` (0.99) in either regime are excluded from selection.

## Notes on Figures 3 and 4 (`fig3_fig4_driver_panels.R`)
- Panel A: AIC and adjusted R2 come from the model-selection table. `model_source = "pipeline"` uses the full run (`results/<run>/TempSelectionSummary.csv`); `"published"` plots directly from the supplementary model table (`Figure3/4_metadata.xlsx`) without re-running the search. The table, n and p are then taken as supplied (p is refitted if the table has no p column). Axis limits, breaks and axis breaks are derived from the range of the plotted values (whole table for the main plot, selected models for the inset); adapt them if another layout is wanted. Inset layers encode p (fill), k and n (outer rings). n is the number of observations of each refitted model.
- Panel B: LMG importances are calculated from the refitted models of the chosen threshold (`importance_source = "calculate"`), or read from a sheet computed before (`"sheet"`). Squares ordered by importance; colour = share (%) of the model R2 explained by each predictor.
- Panel C: PCA variables are those explaining at least 2% (Fig 3) or 5% (Fig 4) of the model R2 in either model; `DA` is left out of the Fig 3 PCA. Arrows: red = important only above the threshold, blue = only below, grey = both. The response is added with `vegan::envfit` (999 permutations, seed 1).
- The script saves a check of the tabulated model statistics against a refit (`results/<run>/*_refit_check.csv`).
- Both responses use identical definitions (`stats::AIC()`, Below = `Temp < T`, Above = `Temp >= T`, overfit filter 0.99, n = observations of the model); only target, predictors and thresholds differ. Running the N2 search therefore regenerates the Fig 3 table and may differ from earlier runs that used other versions of the code or of the taxonomy variables.
- No sites are excluded by ID: a site is left out of a panel only when the value needed for it is missing.

## Caveats
- Subsets of all predictors grow as 2^p: runtime and memory scale accordingly.
- With small n relative to k (e.g. n = 15, k = 12) adjusted R² is near 1 and not reliable.
- Selection among many models inflates apparent fit. P-values are from the selected model and are not corrected for the search.

## Requirements
R >= 4.2; ggplot2 >= 3.5.0 (Fig 1F); readxl, dplyr, tidyr, scales, reshape2, vegan, relaimpo, ggbreak, future.apply, tictoc.
