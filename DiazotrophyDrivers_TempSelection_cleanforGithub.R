###############################################################################
# NCD drivers: temperature-threshold model selection
#
# Pipeline
#   PART 1  Exhaustive subset regressions below/above each temperature threshold
#   PART 2  Post-processing: AICc, validity filters, best model per threshold,
#           p-values (refit)
#   PART 3  Bridge: best-model table -> `model_table` (long format)
#   PART 4  Figures: AIC vs adjusted R2, insets (p, n, k encodings)
#   PART 5  Relative importance + coefficient plots for the chosen threshold
#
# Requirements
#   - Run from the project root (working directory = repo root).
#   - R/calculate_Reg_v2.R must be in the repo. It must return a named vector
#     with RMSE, AdjR2, AIC, n.
#   - Input: data/source_data.csv
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
  library(ggplot2)
  library(ggbreak)
  library(relaimpo)
  library(gridExtra)
  library(ggstatsplot)
})



## Config -- edit here only
CFG <- list(
  data_file    = file.path("data", "source_data.csv"),
  reg_function = file.path("R", "calculate_Reg_v2.R"),
  res_dir      = "results",
  fig_dir      = "figures",
  target       = "non_cyanoMM",   # or "N2fixRates"
  temp_var     = "Temp",
  thresholds   = 8:12,            # degrees C; same set is used in every step
  min_preds    = 2,               # smallest predictor-set size
  chunk_size   = 50000,
  workers      = 12,              # adjust to available cores
  overfit_cut  = 0.99,            # drop models with AdjR2 > cutoff in either regime
  winning_thr  = 9,               # threshold carried into PART 5
  zoom_ids     = paste0("Below", 8:12),  # models shown in the inset figures
  run_search   = TRUE             # FALSE = reuse TempThreshold_<thr>.csv in res_dir
)

dir.create(CFG$res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(CFG$fig_dir, showWarnings = FALSE, recursive = TRUE)
res_path <- function(f) file.path(CFG$res_dir, f)
fig_path <- function(f) file.path(CFG$fig_dir, f)

col_above <- "#800026"
col_below <- "#001a66"

# ---- Data -------------------------------------------------------------------
data <- read.csv(CFG$data_file)
stopifnot(CFG$target %in% names(data), CFG$temp_var %in% names(data))
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
# NOTE: AIC from calculate_Reg uses 2k; the AICc term uses k_AIC = k + 1.
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
      stringsAsFactors = FALSE
    )
  }
  mt <- rbind(one("Above"), one("Below"))
  mt$Group <- factor(mt$Group, levels = c("Above", "Below"))
  rownames(mt) <- NULL
  mt
}

model_table <- build_model_table(best_df)
write.csv(model_table, res_path("TempSelectionNCDSummary.csv"), row.names = FALSE)

# Optional check against a previously saved table
legacy <- file.path("data", "TempSelectionNCDSummary.csv")
if (file.exists(legacy)) {
  old <- read.csv(legacy)
  cmp <- merge(model_table, old, by = "ID", suffixes = c(".new", ".old"))
  for (v in c("AdjR2", "n", "k", "RSS", "AIC", "p")) {
    d <- max(abs(cmp[[paste0(v, ".new")]] - cmp[[paste0(v, ".old")]]), na.rm = TRUE)
    message(sprintf("Check vs legacy table | %-5s max abs diff = %.3g", v, d))
  }
}

###############################################################################
# PART 4. Figures: AIC vs adjusted R2
###############################################################################
pal_col  <- c(Above = col_above, Below = col_below)
pal_fill <- c(Above = scales::alpha(col_above, 0.3), Below = scales::alpha(col_below, 0.3))

## 4.1 Scatter with regression lines
p1 <- ggplot(model_table, aes(AIC, AdjR2, color = Group)) +
  geom_point(alpha = 0.6, size = 5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  scale_color_manual(values = pal_col) +
  labs(x = "AIC", y = "Adjusted R²", color = "Group") +
  theme_classic()
ggsave(fig_path("AICvsR2.png"), p1, width = 10, height = 8, dpi = 600)

## 4.2 Threshold-labelled bubbles with broken x axis
p2 <- ggplot(model_table, aes(AIC, AdjR2)) +
  geom_point(aes(fill = Group, color = Group), size = 45, shape = 21, stroke = 1) +
  geom_text(aes(label = Temp), vjust = 0.35, size = 15, color = "black") +
  scale_fill_manual(values = pal_fill) +
  scale_color_manual(values = pal_col) +
  labs(x = "AIC", y = "Adjusted R²", color = "Group", fill = "Group") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_break(c(80, 150)) +
  scale_x_continuous(limits = c(-50, 300),
                     breaks = c(-50, 0, 50, 150, 200, 250, 300))
ggsave(fig_path("AICvsR2_NCD_TempLabels.png"), p2, width = 25, height = 15, dpi = 600)

## 4.3 Insets (models listed in CFG$zoom_ids)
model_table_sub <- subset(model_table, ID %in% CFG$zoom_ids)

inset_axes <- list(
  theme_classic(),
  scale_x_continuous(limits = c(-50, 50), breaks = seq(-50, 50, by = 10)),
  scale_y_continuous(limits = c(0.75, 1.02), breaks = seq(0.8, 1, by = 0.1)),
  scale_y_break(c(0.82, 0.95))
)

# 4.3a p-value as fill gradient
p_zoom_p <- ggplot(model_table_sub, aes(AIC, AdjR2)) +
  geom_point(aes(fill = p), size = 45, shape = 21, stroke = 0, alpha = 0.8) +
  geom_text(aes(label = Temp), vjust = 0.35, size = 15, color = "black") +
  scale_fill_gradient(low = "#B3D9FF", high = "#001a66", name = "p") +
  labs(x = "AIC", y = "Adjusted R²") +
  inset_axes
ggsave(fig_path("AICvsR2_NCD_TempLabels_zoomed_p-gradientfill.png"),
       p_zoom_p, width = 25, height = 15, dpi = 600)

# 4.3b sample size as outline gradient
p_zoom_n <- ggplot(model_table_sub, aes(AIC, AdjR2)) +
  geom_point(aes(color = n), size = 45, shape = 21, stroke = 6, fill = NA) +
  scale_color_gradient(low = "#fdefb2", high = "#edc001", name = "n") +
  labs(x = "AIC", y = "Adjusted R²") +
  inset_axes
ggsave(fig_path("AICvsR2_NCD_n-outline_zoomed.png"),
       p_zoom_n, width = 25, height = 15, dpi = 600)

# 4.3c number of predictors as outline gradient
p_zoom_k <- ggplot(model_table_sub, aes(AIC, AdjR2)) +
  geom_point(aes(color = k), size = 45, shape = 21, stroke = 6, fill = NA) +
  scale_color_gradient(low = "lightgrey", high = "black", name = "k") +
  labs(x = "AIC", y = "Adjusted R²") +
  inset_axes
ggsave(fig_path("AICvsR2_NCD_k-outline_zoomed.png"),
       p_zoom_k, width = 25, height = 15, dpi = 600)

###############################################################################
# PART 5. Winning threshold: relative importance and coefficient plots
###############################################################################
thr <- CFG$winning_thr
row <- best_df[best_df$Threshold == thr, ]
stopifnot(nrow(row) == 1)

data_below_scaled <- as.data.frame(scale(subset(data, data[[CFG$temp_var]] <  thr)))
data_above_scaled <- as.data.frame(scale(subset(data, data[[CFG$temp_var]] >= thr)))

fit_below <- lm(as.formula(paste(target, "~", row$Combo_Below)), data = data_below_scaled)
fit_above <- lm(as.formula(paste(target, "~", row$Combo_Above)), data = data_above_scaled)

cat("\n== Below model (Temp <", thr, ") ==\n"); print(summary(fit_below))
cat("\n== Above model (Temp >=", thr, ") ==\n"); print(summary(fit_above))

## 5.1 Relative importance (LMG, absolute shares of R2)
relimp_tbl <- function(fit, label) {
  ri <- calc.relimp(fit, rela = FALSE)
  data.frame(Predictor = names(ri@lmg), Rel_Importance = as.numeric(ri@lmg),
             Regime = label, stringsAsFactors = FALSE)
}
relimp_below_df <- relimp_tbl(fit_below, paste0("Below ", thr, "°C"))
relimp_above_df <- relimp_tbl(fit_above, paste0("Above ", thr, "°C"))
relimp_df <- rbind(relimp_below_df, relimp_above_df)
relimp_df <- relimp_df[order(relimp_df$Regime, -relimp_df$Rel_Importance), ]
write.csv(relimp_df, res_path(paste0("NCD_RelImp_Threshold", thr, ".csv")), row.names = FALSE)

relimp_plot <- function(d, lo, hi, title) {
  ggplot(d, aes(reorder(Predictor, Rel_Importance), Rel_Importance * 100,
                fill = Rel_Importance)) +
    geom_col() + coord_flip() +
    scale_fill_gradient(low = lo, high = hi) +
    labs(x = "Predictor", y = "Relative importance (%)", title = title) +
    theme_bw(base_size = 12) + theme(legend.position = "none")
}
pb <- relimp_plot(relimp_below_df, "aliceblue", col_below, paste0("Below ", thr, "°C"))
pa <- relimp_plot(relimp_above_df, "mistyrose", col_above, paste0("Above ", thr, "°C"))
relimp_combined <- gridExtra::arrangeGrob(pb, pa, ncol = 2)
ggsave(fig_path(paste0("NCD_RelImp_Threshold", thr, ".png")),
       relimp_combined, width = 12, height = 6, dpi = 600)

## 5.2 Coefficient plots
coef_plot <- function(fit, col) {
  ggcoefstats(fit,
              errorbar.args = list(color = "lightsteelblue4", linewidth = 1.2),
              point.args    = list(color = col, size = 6),
              stats.labels  = FALSE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}
ggsave(fig_path(paste0("Best_", thr, "_below_regplot_raw.png")),
       coef_plot(fit_below, col_below), width = 12, height = 8, dpi = 600)
ggsave(fig_path(paste0("Best_", thr, "_above_regplot_raw.png")),
       coef_plot(fit_above, col_above), width = 12, height = 8, dpi = 600)

sessionInfo()
