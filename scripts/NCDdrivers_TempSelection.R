###############################################################################
# Temperature-threshold model selection (NCD prevalence, Fig 4; or N2 fixation, Fig 3)
# Set RUN below: "NCD" or "N2fix". Same code, different target / predictors / thresholds.
#
# Pipeline
#   PART 1  Exhaustive subset regressions below/above each temperature threshold
#   PART 2  Post-processing: AICc, validity filters, best model per threshold,
#           p-values (refit)
#   PART 3  Bridge: best-model table -> `model_table` (long format)
#   (figures: scripts/fig3_fig4_driver_panels.R)
#
# Requirements
#   - readxl (reads the supplementary workbook)
#   - Run from the project root (working directory = repo root).
#   - R/calculate_Reg_v2.R must be in the repo. It must return a named vector
#     with RMSE, AdjR2, AIC, n.
#   - Input: data/20261019_GFREDIZ23_metadata.xlsx (supplementary full metadata)
#
# Both runs use exactly the same definitions: AIC = stats::AIC() (R/calculate_Reg_v2.R),
# Below = Temp < T, Above = Temp >= T, overfit filter AdjR2 <= 0.99, n = observations of the model.
# Only target, predictor set and thresholds differ (PRESET).
#
# Provenance of numbers (see README):
#   measured   : data columns
#   computed   : RMSE, AdjR2, AIC, n, k, p (from lm() on z-scored data)
#   derived    : RSS = RMSE^2 * n ; AICc (formula below)
###############################################################################

# ---- 0. Setup ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(future.apply)
  library(tictoc)
})

## Config -- edit here only
RUN <- "NCD"      # "NCD" (Fig 4: target non_cyanoMM) or "N2fix" (Fig 3: target N2FixationRates)
source("R/supplement_io.R")
PRESET <- list(
  NCD   = list(run_name = "fig4", target = "non_cyanoMM", predictors = NCD_PREDICTORS,
               thresholds = 8:12, winning_thr = 9, overfit_cut = 0.99,
               published_sheet = "NCD model selection params", published_file = "Figure4_metadata.xlsx"),
  N2fix = list(run_name = "fig3", target = "N2FixationRates", predictors = N2_PREDICTORS,
               thresholds = 15:25, winning_thr = 21, overfit_cut = 0.99,
               published_sheet = "Model selection params", published_file = "Figure3_metadata.xlsx")
)[[RUN]]

CFG <- list(
  data_file    = file.path("data", "20261019_GFREDIZ23_metadata.xlsx"),
  published_table = file.path("data", PRESET$published_file),   # optional: reproducibility check
  published_sheet = PRESET$published_sheet,
  reg_function = file.path("R", "calculate_Reg_v2.R"),
  res_dir      = file.path("results", PRESET$run_name),
  target       = PRESET$target,
  predictors   = PRESET$predictors,  # Temp is both a predictor and the split variable
  temp_var     = "Temp",
  thresholds   = PRESET$thresholds,  # degrees C; same set is used in every step
  min_preds    = 2,                  # smallest predictor-set size
  chunk_size   = 50000,
  workers      = 12,                 # adjust to available cores
  overfit_cut  = PRESET$overfit_cut, # drop models with AdjR2 > cutoff in either regime
  run_search   = TRUE                # FALSE = reuse TempThreshold_<thr>.csv in res_dir
)

dir.create(CFG$res_dir, showWarnings = FALSE, recursive = TRUE)
res_path <- function(f) file.path(CFG$res_dir, f)


# ---- Data -------------------------------------------------------------------
full <- read_full_metadata(CFG$data_file)
stopifnot(CFG$target %in% names(full), all(CFG$predictors %in% names(full)),
          CFG$temp_var %in% CFG$predictors)
data <- full[, c(CFG$target, CFG$predictors)]       # missing values are handled per model by lm()
source(CFG$reg_function)  # defines calculate_Reg(formula, data, target)

target     <- CFG$target
predictors <- setdiff(names(data), target)

###############################################################################
# PART 1. Exhaustive subset search
###############################################################################
if (CFG$run_search) {

  # All predictor combinations of size min_preds .. n
  all_combos <- unlist(
    lapply(CFG$min_preds:length(predictors),
           function(k) combn(predictors, k, simplify = FALSE)),
    recursive = FALSE
  )
  chunks <- split(all_combos, ceiling(seq_along(all_combos) / CFG$chunk_size))

  plan(multisession, workers = CFG$workers)

  for (thr in CFG$thresholds) {
    message("Processing threshold = ", thr)
    tic(paste("Threshold", thr))

    # Split by temperature and z-score within each regime
    data_below_scaled <- as.data.frame(scale(subset(data, data[[CFG$temp_var]] <  thr)))
    data_above_scaled <- as.data.frame(scale(subset(data, data[[CFG$temp_var]] >= thr)))

    chunk_results <- future_lapply(chunks, function(chunk) {
      out <- lapply(chunk, function(vars) {
        f  <- as.formula(paste(target, "~", paste(vars, collapse = "+")))
        sb <- calculate_Reg(f, data_below_scaled, target)
        sa <- calculate_Reg(f, data_above_scaled, target)
        data.frame(
          Threshold   = thr,
          Combination = paste(vars, collapse = "+"),
          k           = length(vars),
          RMSE_Below  = sb["RMSE"], AdjR2_Below = sb["AdjR2"],
          AIC_Below   = sb["AIC"],  n_Below     = sb["n"],
          RMSE_Above  = sa["RMSE"], AdjR2_Above = sa["AdjR2"],
          AIC_Above   = sa["AIC"],  n_Above     = sa["n"],
          stringsAsFactors = FALSE, row.names = NULL
        )
      })
      do.call(rbind, out)
    }, future.seed = NULL)

    threshold_df <- do.call(rbind, chunk_results)
    write.csv(threshold_df, res_path(paste0("TempThreshold_", thr, ".csv")),
              row.names = FALSE)
    toc()
    rm(chunk_results, threshold_df); gc()
  }
  plan(sequential)
}

###############################################################################
# PART 2. Post-processing and best-model selection
###############################################################################

## 2.1 Validity filter and AICc -----------------------------------------------
# AIC is set to NA when n <= k_AIC (k_AIC = k + 1, i.e. predictors + intercept);
# the raw table contains -Inf AIC for such models.
# AICc = AIC + 2*k_AIC*(k_AIC + 1) / (n - k_AIC - 1).
# It is finite only if n - k_AIC - 1 > 0; otherwise AICc is set to NA.
# Rationale: some models have n = k + 2 with AdjR2 ~ 1; AICc penalizes these.
# NOTE: AIC comes from stats::AIC(lm), which counts k predictors + intercept + residual
# variance (k + 2 parameters). The AICc correction below uses K = k + 1 (as in the
# original analysis). Best models are chosen by AIC; AICc is used only as a validity filter.
add_aicc <- function(df) {
  df$k_AIC <- df$k + 1
  for (g in c("Below", "Above")) {
    n_col <- paste0("n_", g); aic_col <- paste0("AIC_", g)
    aicc_col <- paste0("AICc_", g)
    df[[aic_col]][df[[n_col]] <= df$k_AIC] <- NA
    df[[aicc_col]] <- df[[aic_col]] +
      (2 * df$k_AIC * (df$k_AIC + 1)) / (df[[n_col]] - df$k_AIC - 1)
    df[[aicc_col]][df[[n_col]] <= df$k_AIC + 1] <- NA
  }
  df
}

valid_cols <- c("AdjR2_Below", "AdjR2_Above", "AIC_Below", "AIC_Above",
                "AICc_Below", "AICc_Above")

threshold_tables <- list()
for (thr in CFG$thresholds) {
  df       <- add_aicc(read.csv(res_path(paste0("TempThreshold_", thr, ".csv"))))
  df_clean <- df[complete.cases(df[, valid_cols]), ]
  write.csv(df,       res_path(paste0("TempThreshold_", thr, "_fixed.csv")), row.names = FALSE)
  write.csv(df_clean, res_path(paste0("TempThreshold_", thr, "_clean.csv")), row.names = FALSE)
  threshold_tables[[paste0("T", thr)]] <- df_clean
  message("Threshold ", thr, ": ", nrow(df_clean), " valid models")
}

## 2.2 Best model per regime and threshold (min AIC) ---------------------------
# Overfit filter: models with AdjR2 > overfit_cut in EITHER regime are removed.
get_best <- function(d, g) {
  aic <- paste0("AIC_", g)
  if (nrow(d) == 0 || all(is.na(d[[aic]]))) {
    message("  WARNING: no valid models for regime ", g)
    return(data.frame(Combination = NA, k = NA, AIC = NA, AdjR2 = NA,
                      RMSE = NA, n = NA, stringsAsFactors = FALSE))
  }
  r <- d[which.min(d[[aic]]), ]
  data.frame(Combination = r$Combination, k = r$k, AIC = r[[aic]],
             AdjR2 = r[[paste0("AdjR2_", g)]], RMSE = r[[paste0("RMSE_", g)]],
             n = r[[paste0("n_", g)]], stringsAsFactors = FALSE)
}

best_df <- do.call(rbind, lapply(names(threshold_tables), function(nm) {
  df <- threshold_tables[[nm]]
  df <- df[df$AdjR2_Below <= CFG$overfit_cut & df$AdjR2_Above <= CFG$overfit_cut, ]
  message(nm, ": ", nrow(df), " models after overfit filter")
  a <- get_best(df, "Above"); b <- get_best(df, "Below")
  data.frame(
    Threshold   = as.numeric(sub("T", "", nm)),
    Combo_Above = a$Combination, k_Above = a$k, AIC_Above = a$AIC,
    AdjR2_Above = a$AdjR2, RMSE_Above = a$RMSE, n_Above = a$n,
    Combo_Below = b$Combination, k_Below = b$k, AIC_Below = b$AIC,
    AdjR2_Below = b$AdjR2, RMSE_Below = b$RMSE, n_Below = b$n,
    Combined_AIC = a$AIC + b$AIC, stringsAsFactors = FALSE
  )
}))
best_df <- best_df[order(best_df$Threshold), ]
rownames(best_df) <- NULL
write.csv(best_df, res_path("NCD_AIC_BestModels_AllThresholds.csv"), row.names = FALSE)

cat("\n== Ranked by ABOVE AIC ==\n")
print(best_df[order(best_df$AIC_Above),
              c("Threshold", "AIC_Above", "AdjR2_Above", "k_Above", "n_Above",
                "AIC_Below", "AdjR2_Below", "Combined_AIC")])
cat("\n== Ranked by BELOW AIC ==\n")
print(best_df[order(best_df$AIC_Below),
              c("Threshold", "AIC_Below", "AdjR2_Below", "k_Below", "n_Below",
                "AIC_Above", "AdjR2_Above", "Combined_AIC")])
cat("\n== Ranked by COMBINED AIC ==\n")
print(best_df[order(best_df$Combined_AIC),
              c("Threshold", "Combined_AIC", "AIC_Above", "AdjR2_Above",
                "AIC_Below", "AdjR2_Below")])

## 2.3 p-values: refit best models, overall F-test -----------------------------
get_pvalue <- function(combo, scaled_data) {
  if (is.na(combo)) return(NA_real_)
  tryCatch({
    s <- summary(lm(as.formula(paste(target, "~", combo)), data = scaled_data))
    unname(pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE))
  }, error = function(e) NA_real_)
}

pvals <- do.call(rbind, lapply(seq_len(nrow(best_df)), function(i) {
  thr <- best_df$Threshold[i]
  below <- as.data.frame(scale(subset(data, data[[CFG$temp_var]] <  thr)))
  above <- as.data.frame(scale(subset(data, data[[CFG$temp_var]] >= thr)))
  data.frame(pval_Above = get_pvalue(best_df$Combo_Above[i], above),
             pval_Below = get_pvalue(best_df$Combo_Below[i], below))
}))
best_df <- cbind(best_df, pvals)
write.csv(best_df, res_path("NCD_AIC_BestModels_AllThresholds_withPvals.csv"),
          row.names = FALSE)

###############################################################################
# PART 3. BRIDGE: best_df  ->  model_table (input to the figures in PART 4)
#
# One row per regime x threshold:
#   ID, Group, Temp, AdjR2, n, k, RSS, AIC, p
# RSS is derived: RSS = RMSE^2 * n (checked on TempSelectionSummary_AddedStatistics.csv).
###############################################################################
build_model_table <- function(best_df) {
  one <- function(g) {
    n <- best_df[[paste0("n_", g)]]
    data.frame(
      ID    = paste0(g, best_df$Threshold),
      Group = g,
      Temp  = best_df$Threshold,
      AdjR2 = best_df[[paste0("AdjR2_", g)]],
      n     = n,
      k     = best_df[[paste0("k_", g)]],
      RSS   = best_df[[paste0("RMSE_", g)]]^2 * n,
      AIC   = best_df[[paste0("AIC_", g)]],
      p     = best_df[[paste0("pval_", g)]],
      Combination = best_df[[paste0("Combo_", g)]],
      stringsAsFactors = FALSE
    )
  }
  mt <- rbind(one("Above"), one("Below"))
  mt$Group <- factor(mt$Group, levels = c("Above", "Below"))
  rownames(mt) <- NULL
  mt
}

model_table <- build_model_table(best_df)
write.csv(model_table, res_path("TempSelectionSummary.csv"), row.names = FALSE)

# Reproducibility check against the published (supplementary) model table
if (file.exists(CFG$published_table)) {
  pub <- read_model_table(CFG$published_table, CFG$published_sheet)
  cmp <- merge(model_table, pub[, c("ID", "AdjR2", "n", "k", "RSS", "AIC")], by = "ID",
               suffixes = c(".new", ".published"))
  for (v in c("AdjR2", "n", "k", "RSS", "AIC")) {
    dd <- max(abs(cmp[[paste0(v, ".new")]] - cmp[[paste0(v, ".published")]]), na.rm = TRUE)
    message(sprintf("Check vs published table | %-5s max abs diff = %.3g", v, dd))
  }
}

# Figures for these results: scripts/fig3_fig4_driver_panels.R
sessionInfo()
