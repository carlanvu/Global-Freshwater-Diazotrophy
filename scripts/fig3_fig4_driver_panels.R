###############################################################################
# Figures 3 and 4: temperature-threshold models and their drivers
#   A  AIC vs adjusted R2 for every candidate threshold (+ three inset layers)
#   B  variable importance (LMG share of model R2), above- and below-threshold models
#   C  PCA biplot of the important predictors + envfit of the response
#
# Fig 3 = N2 fixation rates (highlight: above threshold)   run "fig3"
# Fig 4 = NCD prevalence    (highlight: below threshold)   run "fig4"
# Same code for both; only the settings block differs. Set FIG below.
#
# Run from the repo root:  source("scripts/fig3_fig4_driver_panels.R")
# Inputs
#   data/20261019_GFREDIZ23_metadata.xlsx   full metadata (used for refits, importance, PCA)
# Two independent choices (see CFG below):
#   model_source      "pipeline"  results/<FIG>/TempSelectionSummary.csv, written by
#                                 scripts/NCDdrivers_TempSelection.R (RUN = "N2fix" -> fig3, "NCD" -> fig4)
#                     "published" the model-selection sheet of data/Figure3_metadata.xlsx / Figure4_metadata.xlsx
#   importance_source "calculate" LMG calculated from the refitted models (relaimpo)
#                     "sheet"     importance values computed before (sheet of Figure3/4_metadata.xlsx)
# The final layouts (inset overlays, bar labels, significance marks) were assembled in
# Canva from the layers saved in figures/.
###############################################################################
FIG <- "fig4"      # "fig3" or "fig4"

suppressPackageStartupMessages({ library(ggplot2); library(vegan); library(readxl); library(relaimpo) })
has_break <- requireNamespace("ggbreak", quietly = TRUE)
if (!has_break) warning("'ggbreak' not installed: axis breaks are omitted.")
source("R/figure_helpers.R"); source("R/supplement_io.R")

SET <- list(
  fig3 = list(label = "N\u2082 fixation", target = "N2FixationRates", thr = 21,
              pct_cut = 0.02,              # important = explains >= 2% of the model R2
              pca_exclude = "DA",          # DA dominates both models; left out of the PCA
              envfit_resp = "N2FixationRates",
              zoom_ids = c("Above18", "Above19", "Above20", "Above21"),
              published_file = file.path("data", "Figure3_metadata.xlsx"),
              published_sheet = "Model selection params", importance_sheet = "Selected model var imp"),
  fig4 = list(label = "NCD prevalence", target = "non_cyanoMM", thr = 9,
              pct_cut = 0.05,              # >= 5% of the model R2
              pca_exclude = character(0),
              envfit_resp = "non_cyanoMM",
              zoom_ids = paste0("Below", 8:12),
              published_file = file.path("data", "Figure4_metadata.xlsx"),
              published_sheet = "NCD model selection params", importance_sheet = "NCD selected model var imp")
)[[FIG]]

CFG <- list(
  meta_file  = file.path("data", "20261019_GFREDIZ23_metadata.xlsx"),
  temp_var   = "Temp",
  fig_dir    = file.path("figures", FIG),
  res_dir    = file.path("results", FIG),
  model_source = "pipeline",      # "pipeline" (full run) or "published" (supplementary model table)
  # variable importance: "calculate" = LMG from refitted models (relaimpo)
  #                      "sheet"     = values computed before, read from the supplement sheet
  importance_source = "calculate",
  recompute_n = NA,               # n = observations of each refitted model; NA = TRUE for "pipeline", FALSE for "published"
  seed = 1, nperm = 999
)
dir.create(CFG$res_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(CFG$seed)
fp <- function(f) out_path(CFG$fig_dir, paste0(FIG, f))
rp <- function(f) file.path(CFG$res_dir, paste0(FIG, f))

pal_pt   <- c(Above = COL$above, Below = COL$below)
pal_fill <- c(Above = alpha(COL$above, 0.3), Below = alpha(COL$below, 0.3))
meta <- read_full_metadata(CFG$meta_file)
thr  <- SET$thr
num_cols <- names(meta)[vapply(meta, is.numeric, logical(1))]
# z-score every numeric column within the regime (as in the model-selection pipeline)
subset_regime <- function(g, t) as.data.frame(scale(
  if (g == "Below") meta[meta[[CFG$temp_var]] <  t, num_cols] else meta[meta[[CFG$temp_var]] >= t, num_cols]))
fit_model <- function(g, t, combo)
  lm(as.formula(paste(SET$target, "~", combo)), data = subset_regime(g, t))

###############################################################################
# Model table
###############################################################################
if (CFG$model_source == "pipeline") {
  mt <- read.csv(file.path("results", FIG, "TempSelectionSummary.csv"))   # from the model-selection script
  mt$Group <- factor(mt$Group, levels = c("Above", "Below"))
} else {
  mt <- read_model_table(SET$published_file, SET$published_sheet)         # supplementary table
}
fits <- Map(function(g, t, cb) fit_model(as.character(g), t, cb), mt$Group, mt$Temp, mt$Combination)
sm   <- lapply(fits, summary)
mt$p_refit <- vapply(sm, function(s) unname(pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)), numeric(1))
chk <- data.frame(ID = mt$ID, AdjR2 = mt$AdjR2, AdjR2_refit = vapply(sm, function(s) s$adj.r.squared, numeric(1)),
                  n = mt$n, n_refit = vapply(fits, nobs, numeric(1)))
write.csv(chk, rp("_refit_check.csv"), row.names = FALSE)
if (any(abs(chk$AdjR2 - chk$AdjR2_refit) > 0.01))
  warning("Some models do not refit to the tabulated adjusted R2 (see ", rp("_refit_check.csv"), ").")
if (isTRUE(CFG$recompute_n) || (is.na(CFG$recompute_n) && CFG$model_source == "pipeline")) mt$n <- chk$n_refit
pcol <- intersect(c("p", "p-value", "P"), names(mt))
mt$p <- if (length(pcol)) as.numeric(mt[[pcol[1]]]) else mt$p_refit
write.csv(mt[, c("ID", "Group", "Temp", "AdjR2", "n", "k", "RSS", "AIC", "p", "Combination")],
          rp("_model_table.csv"), row.names = FALSE)

###############################################################################
# PANEL A: AIC vs adjusted R2
# Axis limits, breaks and axis breaks are derived from the values being plotted
# (AIC, adjusted R2): the whole table for the main plot and the selected models for
# the inset. Adjust AX_* below if a different layout is wanted.
###############################################################################
nice_limits <- function(x, pad = 0.05) { r <- range(x, na.rm = TRUE); d <- diff(r); c(r[1] - pad * d, r[2] + pad * d) }
# largest gap between consecutive sorted values; returns c(lower, upper) of the gap if it
# is wider than `frac` of the total range (candidate location for an axis break)
find_gap <- function(x, frac) {
  x <- sort(unique(x)); if (length(x) < 3) return(NULL)
  g <- diff(x); i <- which.max(g)
  if (g[i] > frac * diff(range(x))) c(x[i], x[i + 1]) else NULL
}
break_inside <- function(gap, shrink = 0.15) {          # keep the break away from the points
  w <- diff(gap); c(gap[1] + shrink * w, gap[2] - shrink * w)
}
auto_axes <- function(d, gap_frac) {
  xl <- nice_limits(d$AIC); yl <- c(max(0, min(d$AdjR2) - 0.05), min(1.02, max(d$AdjR2) + 0.05))
  list(xlim = xl, ylim = yl, xticks = pretty(xl, 8), yticks = pretty(yl, 5),
       xbreak = { g <- find_gap(d$AIC, gap_frac); if (is.null(g)) NULL else break_inside(g) },
       ybreak = { g <- find_gap(d$AdjR2, gap_frac); if (is.null(g)) NULL else break_inside(g) })
}
sub <- subset(mt, ID %in% SET$zoom_ids); stopifnot(nrow(sub) == length(SET$zoom_ids))
AX_main <- auto_axes(mt, gap_frac = 0.20)
AX_zoom <- auto_axes(sub, gap_frac = 0.30)
cat("Main plot axes:\n"); str(AX_main); cat("Inset axes:\n"); str(AX_zoom)

p_main <- ggplot(mt, aes(AIC, AdjR2)) +
  geom_point(aes(fill = Group, color = Group), size = 45, shape = 21, stroke = 1) +
  geom_text(aes(label = Temp), vjust = 0.35, size = 15, color = "black") +
  scale_fill_manual(values = pal_fill) + scale_color_manual(values = pal_pt) +
  labs(x = "AIC", y = "Adjusted R\u00b2", color = "Group", fill = "Group") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = AX_main$xlim, breaks = AX_main$xticks) +
  theme_classic()
if (has_break && !is.null(AX_main$xbreak)) p_main <- p_main + ggbreak::scale_x_break(AX_main$xbreak)
ggsave(fp("A_main_AICvsR2.png"), p_main, width = 25, height = 15, dpi = 600)

zoom_axes <- function(p) {
  p <- p + labs(x = "AIC", y = "Adjusted R\u00b2") + theme_classic() +
    scale_x_continuous(limits = AX_zoom$xlim, breaks = AX_zoom$xticks) +
    scale_y_continuous(limits = AX_zoom$ylim, breaks = AX_zoom$yticks)
  if (has_break && !is.null(AX_zoom$ybreak)) p <- p + ggbreak::scale_y_break(AX_zoom$ybreak)
  p
}
z_p <- zoom_axes(ggplot(sub, aes(AIC, AdjR2)) +
  geom_point(aes(fill = p), size = 45, shape = 21, stroke = 0, alpha = 0.8) +
  geom_text(aes(label = Temp), vjust = 0.35, size = 15, color = "black") +
  scale_fill_gradient(low = "#B3D9FF", high = "#001a66", name = "p"))
z_k <- zoom_axes(ggplot(sub, aes(AIC, AdjR2)) +
  geom_point(aes(color = k), size = 45, shape = 21, stroke = 6, fill = NA) +
  scale_color_gradient(low = "lightgrey", high = "black", name = "k"))
z_n <- zoom_axes(ggplot(sub, aes(AIC, AdjR2)) +
  geom_point(aes(color = n), size = 45, shape = 21, stroke = 6, fill = NA) +
  scale_color_gradient(low = "#fdefb2", high = "#edc001", name = "n"))
ggsave(fp("A_zoom_p-fill.png"), z_p, width = 25, height = 15, dpi = 600)
ggsave(fp("A_zoom_k-ring.png"), z_k, width = 25, height = 15, dpi = 600)
ggsave(fp("A_zoom_n-ring.png"), z_n, width = 25, height = 15, dpi = 600)

###############################################################################
# PANEL B: variable importance (LMG) as coloured squares ordered by importance
# Importance = LMG share of the model R2 (absolute); colour/label = % of the model R2
###############################################################################
if (CFG$importance_source == "calculate") {
  lmg_tbl <- function(g) {
    row <- mt[mt$Group == g & mt$Temp == thr, ]; stopifnot(nrow(row) == 1)
    fit <- fit_model(g, thr, row$Combination)
    ri  <- calc.relimp(fit, type = "lmg", rela = FALSE)
    list(tab = data.frame(var = names(ri@lmg), value = as.numeric(ri@lmg)), R2 = summary(fit)$r.squared)
  }
  la <- lmg_tbl("Above"); lb <- lmg_tbl("Below")
  imp <- merge(setNames(la$tab, c("var", "Above")), setNames(lb$tab, c("var", "Below")), by = "var", all = TRUE)
  imp$Variable <- pretty_name(imp$var)
  attr(imp, "total_R2") <- c(Above = la$R2, Below = lb$R2)
} else {
  imp <- read_importance(SET$published_file, SET$importance_sheet)   # values computed before
}
tot <- attr(imp, "total_R2")
imp$pct_Above <- imp$Above / tot["Above"] * 100
imp$pct_Below <- imp$Below / tot["Below"] * 100
write.csv(imp, rp("B_variable_importance.csv"), row.names = FALSE)
write.csv(data.frame(Regime = names(tot), R2 = as.numeric(tot)), rp("B_total_R2.csv"), row.names = FALSE)

tiles <- function(regime, lo, hi) {
  d <- imp[!is.na(imp[[regime]]), ]
  d$pct <- d[[paste0("pct_", regime)]]
  d <- d[order(-d$pct), ]; d$Variable <- factor(d$Variable, levels = rev(d$Variable))
  ggplot(d, aes(x = "Importance", y = Variable, fill = pct)) +
    geom_tile(width = 0.9, height = 0.9, color = "white") +
    geom_text(aes(label = sprintf("%.1f", pct)), size = 4) +
    scale_fill_gradient(low = lo, high = hi, name = "% of model R\u00b2") +
    labs(title = sprintf("%s %s\u00b0C (R\u00b2 = %.2f)", regime, thr, tot[regime]), x = NULL, y = NULL) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())
}
ggsave(fp("B_importance_above.png"), tiles("Above", "mistyrose", COL$above), width = 5, height = 8, dpi = 300)
ggsave(fp("B_importance_below.png"), tiles("Below", "aliceblue", COL$below), width = 5, height = 8, dpi = 300)

###############################################################################
# PANEL C: PCA biplot (+ envfit)
# Variables: those explaining >= pct_cut of the model R2 in either model.
# Arrow colour: important only above (red), only below (blue), both (grey).
###############################################################################
imp$imp_above <- !is.na(imp$pct_Above) & imp$pct_Above >= SET$pct_cut * 100
imp$imp_below <- !is.na(imp$pct_Below) & imp$pct_Below >= SET$pct_cut * 100
if (!"var" %in% names(imp)) imp$var <- to_var(imp$Variable)
sel <- imp[(imp$imp_above | imp$imp_below) & !imp$var %in% SET$pca_exclude, ]
pca_vars <- sel$var
cat("PCA variables (", nrow(sel), "): ", paste(pca_vars, collapse = ", "), "\n", sep = "")
write.csv(sel[, c("Variable", "var", "pct_Above", "pct_Below", "imp_above", "imp_below")], rp("C_pca_variables.csv"), row.names = FALSE)

pc <- na.omit(meta[, unique(c(pca_vars, SET$envfit_resp, CFG$temp_var))])
pc$TempGroup <- factor(ifelse(pc[[CFG$temp_var]] >= thr, "Above", "Below"), levels = c("Above", "Below"))
cat("PCA samples:", nrow(pc), "\n")

pca <- prcomp(scale(pc[, pca_vars]), center = TRUE, scale. = TRUE)
ef  <- envfit(pca, pc[SET$envfit_resp], permutations = CFG$nperm)
print(ef); capture.output(print(ef), file = rp("C_envfit.txt"))

scores <- as.data.frame(pca$x[, 1:2]); scores$TempGroup <- pc$TempGroup
ld <- as.data.frame(pca$rotation[, 1:2]); ld$var <- rownames(ld)
sf <- min(max(abs(scores$PC1)) / max(abs(ld$PC1)), max(abs(scores$PC2)) / max(abs(ld$PC2))) * 0.7
ld$PC1s <- ld$PC1 * sf; ld$PC2s <- ld$PC2 * sf
g <- sel[match(ld$var, sel$var), ]
ld$arrow_group <- ifelse(g$imp_above & g$imp_below, "Shared", ifelse(g$imp_above, "Above", "Below"))
write.csv(ld, rp("C_pca_loadings.csv"), row.names = FALSE)

ve <- summary(pca)$importance[2, 1:2] * 100
p_pca <- ggplot(scores, aes(PC1, PC2, colour = TempGroup, fill = TempGroup)) +
  stat_ellipse(geom = "polygon", type = "norm", alpha = 0.04, linewidth = 0.36, linetype = 2) +
  geom_point(alpha = 0.2, size = 4) +
  geom_segment(data = ld, aes(x = 0, y = 0, xend = PC1s, yend = PC2s, colour = arrow_group),
               arrow = arrow(length = unit(0.25, "cm"), type = "open"), linewidth = 1.2, inherit.aes = FALSE) +
  scale_color_manual(values = c(pal_pt, Shared = "grey50"), name = "Temperature group") +
  scale_fill_manual(values = pal_pt, name = "Temperature group") +
  labs(title = paste("Environmental drivers of", SET$label),
       x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2])) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        legend.title = element_text(face = "bold"), axis.title = element_text(face = "bold"),
        panel.grid = element_blank(), plot.background = element_rect(fill = "white", color = NA))
ggsave(fp("C_PCA_biplot.png"), p_pca, width = 16, height = 8, dpi = 600)

sessionInfo()
