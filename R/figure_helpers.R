# Shared constants and helpers for the figure scripts.
# Source from the repo root:  source("R/figure_helpers.R")

COL <- list(
  above      = "#800026",   # dark red: above-threshold regime
  below      = "#001a66",   # dark blue: below-threshold regime
  lentic_pt  = "goldenrod1",       # bubble plots (#FFC125)
  lentic_box = "darkgoldenrod1",   # boxplot insets (#FFB90F), as in the original code
  lotic      = "magenta4",         # #8B008B
  lentic_vio = "#D4A017",          # Fig 1F
  lotic_vio  = "#6A0572"           # Fig 1F
)

# Latest palettes used for Fig 2 (earlier drafts in the original script were discarded)
PAL_CYANO <- c("Nostocales" = "#1B7745", "Oscillatoriales" = "#66C08F",
               "Chroococcales" = "#A8E6A3")
PAL_NCD <- c("Desulfuromonadales" = "#B22222", "Burkholderiales" = "#FF6347",
             "Pseudomonadales" = "#F08080", "Desulfovibrionales" = "#FFA07A",
             "Rhodocyclales" = "#FFDAB9")

# Mean RGB of a set of hex colours -> one hex colour
mean_hex <- function(cols) {
  m <- rowMeans(grDevices::col2rgb(cols))
  grDevices::rgb(m[1], m[2], m[3], maxColorValue = 255)
}

# Create output dirs and return a path builder
out_path <- function(dir, file) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  file.path(dir, file)
}

# Font: "Arial" if available, otherwise the default sans family
pick_font <- function(preferred = "Arial") {
  if (preferred %in% names(grDevices::pdfFonts()) ||
      requireNamespace("systemfonts", quietly = TRUE) &&
      preferred %in% systemfonts::system_fonts()$family) preferred else "sans"
}

# Read the analysis table. Site IDs are the row names (first column).
read_sites <- function(path) {
  d <- read.csv(path, row.names = 1, check.names = FALSE)
  if ("Type" %in% names(d)) d$Type <- factor(d$Type, levels = c("Lentic", "Lotic"))
  d
}
