###############################################################################
# Figure 1B, 1D (bubble plots) and boxplot insets (1B, 1D, 1E)
#   1B: DA:BA ratio (bubbles) + DA:BA inset
#   1D: N2 fixation per cell (fmol N cell-1 d-1) (bubbles) + inset with axis break
#   1E: Shannon-Wiener index inset
# Maps 1A, 1C, 1E were made in ArcGIS Pro (not reproduced here).
#
# Run from the repo root:  source("scripts/fig1_bubbles_and_insets.R")
# Input : data/20261019_GFREDIZ23_metadata.xlsx  (full metadata)
#         columns used: Type, latitude, longitude, DAtoBA, N2FixationCell_fmol, ShannonMM
###############################################################################
suppressPackageStartupMessages({ library(ggplot2); library(readxl) })
has_break <- requireNamespace("ggbreak", quietly = TRUE)
if (!has_break) warning("Package 'ggbreak' not installed: the axis break in the 1D inset is omitted.")
source("R/figure_helpers.R"); source("R/supplement_io.R")

CFG <- list(
  meta_file = file.path("data", "20261019_GFREDIZ23_metadata.xlsx"),
  fig_dir   = file.path("figures", "fig1"),
  res_dir   = file.path("results", "fig1"),
  font      = pick_font("Arial"),
  # No sites are excluded by ID. A site is left out of a panel only when the value
  # needed for that panel (Type, coordinates or the plotted variable) is missing.
  exclude_sites = character(0)
)
dir.create(CFG$res_dir, showWarnings = FALSE, recursive = TRUE)

d <- read_full_metadata(CFG$meta_file)
need <- c("Type", "latitude", "longitude", "DAtoBA", "N2FixationCell_fmol", "ShannonMM")
stopifnot(all(need %in% names(d)))
d <- d[!rownames(d) %in% CFG$exclude_sites, ]
d$Type <- factor(d$Type, levels = c("Lentic", "Lotic"))
d <- d[!is.na(d$Type), ]
cat("Sites plotted:", nrow(d), "\n")

fill_pt <- c(Lentic = COL$lentic_pt, Lotic = COL$lotic)
col_box <- c(Lentic = COL$lentic_box, Lotic = COL$lotic)

bubble <- function(data, size_var, legend_title, size_range, breaks = waiver()) {
  data <- data[!is.na(data[[size_var]]), ]
  ggplot(data, aes(x = longitude, y = latitude)) +
    geom_point(aes(size = .data[[size_var]], fill = Type),
               shape = 21, color = "black", alpha = 0.25, stroke = 0.5) +
    scale_size_continuous(range = size_range, breaks = breaks) +
    scale_fill_manual(values = fill_pt) +
    theme_classic(base_family = CFG$font) +
    xlab("Longitude") + ylab("Latitude") +
    theme(axis.text = element_text(size = 18, color = "black"),
          axis.title = element_text(size = 20, color = "black"),
          legend.title = element_text(color = "black", size = 18),
          legend.text = element_text(color = "black", size = 18)) +
    guides(size = guide_legend(title = legend_title, override.aes = list(stroke = 1.2)),
           fill = guide_legend(title = "Type", override.aes = list(size = 6, alpha = 0.6)))
}

p1d <- bubble(d, "N2FixationCell_fmol", "N\u2082 fixation\n(femtomol N cell\u207b\u00b9 d\u207b\u00b9)",
              c(2, 32), breaks = c(0, 0.5, 1, 1.5, 2))
ggsave(out_path(CFG$fig_dir, "Fig1D_bubble_N2fixCell.png"), p1d, width = 12, height = 6, dpi = 600)

p1b <- bubble(d, "DAtoBA", "DA:BA", c(2, 32))
ggsave(out_path(CFG$fig_dir, "Fig1B_bubble_DABA.png"), p1b, width = 16, height = 8, dpi = 600)

inset <- function(data, y, ylab) {
  data <- data[!is.na(data[[y]]), ]
  ggplot(data, aes(x = Type, y = .data[[y]], colour = Type)) +
    geom_boxplot(width = 0.8, position = position_dodge(width = 0.8),
                 outlier.shape = NA, linewidth = 1) +
    geom_jitter(width = 0.15, height = 0) +
    scale_colour_manual(values = col_box) +
    ylab(ylab) +
    theme_classic() +
    theme(axis.text = element_text(size = 15, color = "black"),
          axis.title.y = element_text(size = 20), axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "none") +
    coord_flip()
}

in1d <- inset(d, "N2FixationCell_fmol", "N\u2082 fixation\n(fmol N cell\u207b\u00b9 d\u207b\u00b9)") +
  scale_y_continuous(limits = c(0, 3.6), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 3.5, 3.6))
if (has_break) in1d <- in1d + ggbreak::scale_y_break(c(0.2, 3.5))   # two high-value sites above the break
ggsave(out_path(CFG$fig_dir, "Fig1D_inset_N2fixCell.png"), in1d, width = 10, height = 7, dpi = 600)

in1b <- inset(d, "DAtoBA", "DA:BA")
ggsave(out_path(CFG$fig_dir, "Fig1B_inset_DABA.png"), in1b, width = 10, height = 7, dpi = 600)

in1e <- inset(d, "ShannonMM", "Shannon-Wiener index") +
  theme(axis.text = element_text(size = 20, color = "black"))
ggsave(out_path(CFG$fig_dir, "Fig1E_inset_Shannon.png"), in1e, width = 10, height = 6, dpi = 600)

# Wilcoxon rank-sum tests, lentic vs lotic
tests <- do.call(rbind, lapply(c("DAtoBA", "ShannonMM", "N2FixationCell_fmol"), function(v) {
  w <- wilcox.test(d[[v]] ~ d$Type)
  data.frame(variable = v, n_lentic = sum(!is.na(d[[v]]) & d$Type == "Lentic"),
             n_lotic = sum(!is.na(d[[v]]) & d$Type == "Lotic"), W = unname(w$statistic), p = w$p.value)
}))
write.csv(tests, file.path(CFG$res_dir, "Fig1_lentic_vs_lotic_wilcoxon.csv"), row.names = FALSE)
print(tests)
