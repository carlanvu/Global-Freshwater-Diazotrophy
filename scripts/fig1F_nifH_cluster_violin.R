###############################################################################
# Figure 1F: relative abundance of nifH clusters (1-4 = I-IV), lentic vs lotic
#            mirrored violin plot, lotic up (+) and lentic down (-), pseudo-log10 y axis
#
# Based on nifHclustersRelativeAbundance.R (final "p_log" version).
# Upstream (nifHclustersRelativeAbundance.R, not repeated here): ASV counts + taxonomy;
# ASVs without a cluster assignment removed; replicate 88C removed -> GF23CountsTaxaClusters_clean.csv.
#
# Run from the repo root:  source("scripts/fig1F_nifH_cluster_violin.R")
# Inputs
#   data/GF23CountsTaxaClusters_clean.csv   ASV x sample counts + AUID, Phylum, Class, Oder,
#                                           Family, Genus, Cluster (sample columns = site IDs)
#   data/20261019_GFREDIZ23_metadata.xlsx   full metadata (Type)
# Requires ggplot2 >= 3.5.0 (violin `bounds` argument)
###############################################################################
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(readxl); library(scales) })
source("R/figure_helpers.R")
stopifnot(packageVersion("ggplot2") >= "3.5.0")

CFG <- list(
  cluster_file = file.path("data", "GF23CountsTaxaClusters_clean.csv"),
  meta_file    = file.path("data", "20261019_GFREDIZ23_metadata.xlsx"),
  fig_dir      = file.path("figures", "fig1"),
  res_dir      = file.path("results", "fig1"),
  exclude_sites = character(0),            # none: sites are dropped only when a value is missing
  roman        = TRUE                      # label clusters I-IV (the original script used 1-4)
)
dir.create(CFG$res_dir, showWarnings = FALSE, recursive = TRUE)
source("R/supplement_io.R")

meta   <- read_full_metadata(CFG$meta_file)
counts <- read.csv(CFG$cluster_file, check.names = FALSE)
tax_cols <- c("AUID", "Phylum", "Class", "Oder", "Family", "Genus", "Cluster")
stopifnot(all(tax_cols %in% names(counts)))
sample_cols <- setdiff(intersect(names(counts), rownames(meta)), CFG$exclude_sites)
cat("Samples in the cluster table with metadata:", length(sample_cols), "\n")

# Per-sample relative abundance of each ASV (count / sample total, over ASVs with a cluster),
# then summed per cluster -> one value per sample and cluster
CLUSTER_REL <- counts[, c("Cluster", sample_cols)] %>%
  pivot_longer(-Cluster, names_to = "Sample", values_to = "Count") %>%
  group_by(Sample) %>% mutate(RelAbund = Count / sum(Count, na.rm = TRUE)) %>% ungroup() %>%
  group_by(Sample, Cluster) %>% summarise(RelAbund = sum(RelAbund), .groups = "drop") %>%
  filter(!is.na(RelAbund)) %>%                       # samples without reads in any clustered ASV
  mutate(Type = factor(meta[Sample, "Type"], levels = c("Lentic", "Lotic")),
         Cluster = factor(Cluster, levels = as.character(1:4))) %>%
  filter(!is.na(Type), !is.na(Cluster))

# Summary per cluster x habitat
SUMMARY_REL <- CLUSTER_REL %>%
  group_by(Type, Cluster) %>%
  summarise(Mean = mean(RelAbund), SD = sd(RelAbund), SE = SD / sqrt(dplyr::n()), .groups = "drop")

# Wilcoxon rank-sum test per cluster; FDR (Benjamini-Hochberg) adjustment
ClusterTests <- CLUSTER_REL %>%
  group_by(Cluster) %>%
  summarise(p_value = wilcox.test(RelAbund ~ Type)$p.value, .groups = "drop") %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))
write.csv(SUMMARY_REL, file.path(CFG$res_dir, "Fig1F_cluster_summary.csv"), row.names = FALSE)
write.csv(ClusterTests, file.path(CFG$res_dir, "Fig1F_cluster_wilcoxon_FDR.csv"), row.names = FALSE)
print(as.data.frame(SUMMARY_REL)); print(as.data.frame(ClusterTests))

# Mirrored data: lotic positive, lentic negative
lotic_raw  <- CLUSTER_REL %>% filter(Type == "Lotic")  %>% mutate(y =  RelAbund)
lentic_raw <- CLUSTER_REL %>% filter(Type == "Lentic") %>% mutate(y = -RelAbund)
lotic_sum  <- SUMMARY_REL %>% filter(Type == "Lotic")  %>% mutate(Mean_y =  Mean, ymin =  Mean - SE, ymax =  Mean + SE)
lentic_sum <- SUMMARY_REL %>% filter(Type == "Lentic") %>% mutate(Mean_y = -Mean, ymin = -(Mean + SE), ymax = -(Mean - SE))

col_lotic <- COL$lotic_vio; col_lentic <- COL$lentic_vio
y_max_data <- max(c(lotic_raw$y, abs(lentic_raw$y)), na.rm = TRUE)
y_breaks <- c(-0.8, -0.6, -0.4, -0.2, -0.05, -0.005, 0.005, 0.05, 0.2, 0.4, 0.6, 0.8)
y_breaks <- y_breaks[abs(y_breaks) <= y_max_data * 1.2]

half_violin <- function(data, col, bounds)
  geom_violin(data = data, aes(x = Cluster, y = y), fill = col, colour = col, alpha = 0.25,
              trim = TRUE, scale = "width", width = 0.7, bounds = bounds)   # hard boundary at zero
jit <- function(data, col)
  geom_jitter(data = data, aes(x = Cluster, y = y), colour = col, alpha = 0.55, size = 1.5,
              width = 0.08, height = 0)
xbar <- function(data, col)
  geom_crossbar(data = data, aes(x = Cluster, y = Mean_y, ymin = ymin, ymax = ymax),
                colour = col, fill = col, alpha = 0.7, width = 0.35, linewidth = 0.5, fatten = 2)

p <- ggplot() +
  half_violin(lotic_raw, col_lotic, c(0, Inf)) + half_violin(lentic_raw, col_lentic, c(-Inf, 0)) +
  jit(lotic_raw, col_lotic) + jit(lentic_raw, col_lentic) +
  xbar(lotic_sum, col_lotic) + xbar(lentic_sum, col_lentic) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.7) +
  scale_y_continuous(trans = pseudo_log_trans(sigma = 0.01, base = 10),
                     breaks = y_breaks, labels = as.character(abs(y_breaks)),
                     expand = expansion(mult = 0.05)) +
  annotate("text", x = 0.55, y =  y_max_data * 0.92, label = "Lotic",  hjust = 0,
           colour = col_lotic, fontface = "bold", size = 4.5) +
  annotate("text", x = 0.55, y = -y_max_data * 0.92, label = "Lentic", hjust = 0,
           colour = col_lentic, fontface = "bold", size = 4.5) +
  labs(title = "Relative Abundance of nifH Clusters",
       subtitle = "Lotic (upper) vs Lentic (lower) \u2014 pseudo-log10 scale",
       x = "nifH Cluster", y = "Relative Abundance") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, colour = "grey40", margin = margin(b = 8)),
        axis.title.x = element_text(margin = margin(t = 8)), axis.title.y = element_text(margin = margin(r = 8)),
        axis.text.x = element_text(size = 12, face = "bold"), axis.text.y = element_text(size = 10),
        panel.grid.major.y = element_line(colour = "grey95", linewidth = 0.4),
        panel.grid.major.x = element_blank(), panel.border = element_blank(), legend.position = "none")
if (CFG$roman) p <- p + scale_x_discrete(labels = c("I", "II", "III", "IV"))

ggsave(out_path(CFG$fig_dir, "Fig1F_nifH_MirrorViolin_pseudolog.png"), p, width = 9, height = 7, dpi = 600)
