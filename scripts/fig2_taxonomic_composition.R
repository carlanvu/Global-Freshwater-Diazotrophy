###############################################################################
# Figure 2: cyanobacteria vs non-cyanobacterial diazotrophs (NCDs)
#   2B     pies by habitat type (lentic / lotic): cyano vs NCD, and top-3 orders of each
#   2C-i   pies by N2 fixation category (zero / low / medium / high)
#   2C-ii  boxplots, top 3 cyanobacterial orders per category
#   2C-iii boxplots, top 3 NCD orders per category
# Map 2A and its inset bar chart were made in ArcGIS Pro / Excel (not reproduced here).
#
# Run from the repo root:  source("scripts/fig2_taxonomic_composition.R")
# Inputs
#   data/20261019_GFREDIZ23_metadata.xlsx        full metadata (Type, N2FixationRates, cyanoMM, non_cyanoMM)
#   data/FullTaxCounts_GF23_20250814OutMM.txt    tab-delimited ASV table; row names = ASV IDs;
#       one column per sample (site IDs); last 5 columns = Phylum, Class, Oder, Family, Genus
###############################################################################
suppressPackageStartupMessages({ library(ggplot2); library(reshape2); library(scales) })
source("R/figure_helpers.R"); source("R/supplement_io.R")

CFG <- list(
  meta_file  = file.path("data", "20261019_GFREDIZ23_metadata.xlsx"),
  tax_file   = file.path("data", "FullTaxCounts_GF23_20250814OutMM.txt"),
  fig_dir    = file.path("figures", "fig2"),
  res_dir    = file.path("results", "fig2"),
  tax_cols   = c("Phylum", "Class", "Oder", "Family", "Genus"),
  rank_col   = "Oder",
  top_n      = 3,                              # top orders per category / habitat, cyano and NCD
  cat_levels = c("Zero", "Low", "Medium", "High")
)
dir.create(CFG$res_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Data -------------------------------------------------------------------
meta  <- read_full_metadata(CFG$meta_file)
taxa  <- read.delim(CFG$tax_file, row.names = 1, check.names = FALSE)
stopifnot(all(CFG$tax_cols %in% names(taxa)))
sample_cols <- intersect(setdiff(names(taxa), CFG$tax_cols), rownames(meta))
counts   <- as.matrix(taxa[, sample_cols])                 # ASV x sample
is_cyano <- taxa$Phylum == "Cyanobacteria"

# N2 fixation categories: zero + tertiles of the non-zero rates
n2 <- meta[!is.na(meta$N2FixationRates), , drop = FALSE]
q  <- quantile(n2$N2FixationRates[n2$N2FixationRates > 0], probs = c(1/3, 2/3))
cat_of <- function(x) factor(ifelse(x == 0, "Zero", ifelse(x <= q[1], "Low",
                       ifelse(x <= q[2], "Medium", "High"))), levels = CFG$cat_levels)
meta$Category <- NA; meta$Category[!is.na(meta$N2FixationRates)] <- as.character(cat_of(n2$N2FixationRates))
cat("N2 fixation tertile breaks (non-zero rates):", signif(q, 4), "\n")
sites_by_cat  <- split(rownames(meta), meta$Category)
sites_by_type <- split(rownames(meta), meta$Type)

# Group share per site (% of diazotroph reads, normalised to cyano + NCD)
tot <- meta$cyanoMM + meta$non_cyanoMM
meta$pct_cyano <- meta$cyanoMM / tot * 100
meta$pct_ncd   <- meta$non_cyanoMM / tot * 100

# ---- Order-level proportions WITHIN each group (cyano or NCD) -----------------
# Unclassified ("_") orders removed before normalising; NaN where a sample has no reads in the group
order_prop <- function(rows) {
  agg <- rowsum(counts[rows, , drop = FALSE], group = taxa[[CFG$rank_col]][rows])
  agg <- agg[rownames(agg) != "_", , drop = FALSE]
  out <- as.data.frame(sweep(agg, 2, colSums(agg, na.rm = TRUE), "/"), check.names = FALSE)
  cbind(Oder = rownames(out), out, row.names = NULL)
}
cyano_prop <- order_prop(is_cyano)
NCD_prop   <- order_prop(!is_cyano)
prop <- list(Cyano = cyano_prop, NCD = NCD_prop)

# Paired Wilcoxon: NCD vs cyano total proportion per site (sites with cyano reads only)
tot_df <- na.omit(data.frame(Cyano = colSums(cyano_prop[, -1], na.rm = TRUE) / 1,
                             NCD = colSums(NCD_prop[, -1], na.rm = TRUE) / 1))
tot_df <- tot_df[tot_df$Cyano > 0, ]
print(wilcox.test(tot_df$NCD, tot_df$Cyano, paired = TRUE, alternative = "greater"))

# Subset the proportion table to a set of sites (only sites present in the table)
sub_sites <- function(df, sites) df[, c("Oder", sites[sites %in% colnames(df)]), drop = FALSE]

# Top orders across the sites of a subset (mean proportion, NaN ignored)
get_top <- function(df, n = CFG$top_n) {
  avg <- rowMeans(df[, -1, drop = FALSE], na.rm = TRUE)
  df[order(avg, decreasing = TRUE)[seq_len(min(n, nrow(df)))], , drop = FALSE]
}
to_long <- function(df) {
  d <- melt(df, id.vars = "Oder", variable.name = "Site", value.name = "Proportion")
  d$Proportion[is.nan(d$Proportion)] <- NA
  d[!is.na(d$Proportion), ]                   # sites without reads in the group are dropped
}

# ---- Colours and pies ---------------------------------------------------------
col_cyano <- mean_hex(PAL_CYANO); col_ncd <- mean_hex(PAL_NCD)   # mean RGB of each palette

pie_plot <- function(values, colors, title) {
  d <- data.frame(Group = factor(names(values), levels = names(values)), Percent = as.numeric(values))
  ggplot(d, aes(x = "", y = Percent, fill = Group)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") + scale_fill_manual(values = colors) + theme_void() +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = "none")
}

# ---- 2B: by habitat type ------------------------------------------------------
tab2b <- NULL
for (h in c("Lentic", "Lotic")) {
  s <- sites_by_type[[h]]
  pc <- mean(meta[s, "pct_cyano"], na.rm = TRUE); pn <- mean(meta[s, "pct_ncd"], na.rm = TRUE)
  ggsave(out_path(CFG$fig_dir, paste0("Fig2B_pie_", h, ".png")),
         pie_plot(c(Cyanobacteria = pc, NCD = pn), c(Cyanobacteria = col_cyano, NCD = col_ncd),
                  paste("Habitat type:", h)), width = 5, height = 5, dpi = 600)
  # top-3 orders of each group in this habitat; mean within-group proportion, rescaled to the
  # group's share of the community (individual slices sum to the group total)
  sl <- list()
  for (g in c("Cyano", "NCD")) {
    top <- get_top(sub_sites(prop[[g]], s))
    avg <- aggregate(Proportion ~ Oder, data = to_long(top), FUN = mean, na.rm = TRUE)
    avg$Individual <- avg$Proportion / sum(avg$Proportion) * (if (g == "Cyano") pc else pn) / 100
    avg$Group <- g; avg$Habitat <- h
    sl[[g]] <- avg
  }
  sl <- do.call(rbind, sl); tab2b <- rbind(tab2b, sl)
  pal <- c(PAL_CYANO, PAL_NCD)
  ggsave(out_path(CFG$fig_dir, paste0("Fig2B_pie_orders_", h, ".png")),
         pie_plot(setNames(sl$Individual * 100, sl$Oder), pal[sl$Oder], paste("Habitat type:", h)),
         width = 5, height = 5, dpi = 600)
}
write.csv(tab2b, file.path(CFG$res_dir, "Fig2B_top3_orders_by_habitat.csv"), row.names = FALSE)
print(tab2b)

# ---- 2C-i: pies by N2 fixation category ----------------------------------------
cat_means <- do.call(rbind, lapply(CFG$cat_levels, function(cl) {
  s <- sites_by_cat[[cl]]
  data.frame(Category = cl, n = sum(!is.na(meta[s, "pct_cyano"])),
             pct_cyano = mean(meta[s, "pct_cyano"], na.rm = TRUE), pct_ncd = mean(meta[s, "pct_ncd"], na.rm = TRUE))
}))
write.csv(cat_means, file.path(CFG$res_dir, "Fig2Ci_category_means.csv"), row.names = FALSE)
print(cat_means)
for (i in seq_len(nrow(cat_means)))
  ggsave(out_path(CFG$fig_dir, paste0("Fig2Ci_pie_", cat_means$Category[i], ".png")),
         pie_plot(c(Cyanobacteria = cat_means$pct_cyano[i], NCD = cat_means$pct_ncd[i]),
                  c(Cyanobacteria = col_cyano, NCD = col_ncd),
                  paste("N\u2082 fixation category:", cat_means$Category[i])), width = 5, height = 5, dpi = 600)
for (v in c("pct_cyano", "pct_ncd")) {
  ok <- !is.na(meta$Category) & !is.na(meta[[v]])
  cat("\n", v, "\n"); print(kruskal.test(meta[[v]][ok] ~ factor(meta$Category[ok])))
  print(suppressWarnings(pairwise.wilcox.test(meta[[v]][ok], meta$Category[ok], p.adjust.method = "BH")))
}

# ---- 2C-ii / 2C-iii: top-3 orders per category ---------------------------------
box_plot <- function(d, palette) {
  lev <- names(sort(tapply(d$Proportion, d$Oder, mean), decreasing = TRUE))   # descending mean
  d$Oder <- factor(d$Oder, levels = lev)
  ggplot(d, aes(x = Oder, y = Proportion, fill = Oder)) +
    geom_boxplot(outlier.shape = NA, width = 0.8) +
    geom_jitter(shape = 21, fill = alpha("black", 0.4), color = "black", width = 0.2, size = 8, stroke = NA) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    scale_fill_manual(values = palette, drop = FALSE) +
    labs(y = "Relative abundance", x = "Order") +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
}
all_top <- list()
for (g in c("Cyano", "NCD")) for (cl in CFG$cat_levels) {
  d <- to_long(get_top(sub_sites(prop[[g]], sites_by_cat[[cl]])))
  d$Category <- cl; d$Group <- g; all_top[[paste(g, cl)]] <- d
  pal <- if (g == "Cyano") PAL_CYANO else PAL_NCD
  miss <- setdiff(unique(as.character(d$Oder)), names(pal))
  if (length(miss)) warning(g, "/", cl, ": orders without a palette colour: ", paste(miss, collapse = ", "))
  ggsave(out_path(CFG$fig_dir, sprintf("Fig2C_%s_top3_%s.png", g, cl)), box_plot(d, pal),
         width = 10, height = 14, dpi = 600)
}
write.csv(do.call(rbind, all_top), file.path(CFG$res_dir, "Fig2C_top3_orders_by_category.csv"), row.names = FALSE)
