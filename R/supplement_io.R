# Readers for the supplementary tables (Excel) and the mapping between the
# descriptive names used in the supplement/figures and the column names used in code.
# Requires readxl.   Source from the repo root:  source("R/supplement_io.R")

# descriptive name (supplement)  ->  column name in the full metadata table
VAR_MAP <- c(
  "Latitude" = "latitude", "latitude" = "latitude",
  "Phototrophs abundance" = "PhotoBA", "Heterotrophs abundance" = "Het", "Het" = "Het",
  "DA" = "DA", "DA:BA" = "DAtoBA", "C:N" = "DOCtoDON",
  "Ammonia" = "NH4", "Phosphate" = "P.PO4", "Nitrate" = "N.NO3", "N:P" = "NtoP",
  "Conductivity" = "Conductivity", "Electrical conductivity" = "Conductivity",
  "Dissolved oxygen" = "DO", "pH" = "pH", "Irradiance" = "GHI", "Temperature" = "Temp",
  "No. of wet days" = "Wdays",
  "Diazotrophs richness" = "RichnessMM", "Diazotrophs diversity" = "ShannonMM",
  "N2 Fixation Rates" = "N2FixationRates", "N2 fixation rates" = "N2FixationRates",
  "N2 Fixation Rates per cell" = "N2FixationCell_fmol", "N2 fix per cell" = "N2FixationCell_fmol",
  "cyano" = "cyanoMM", "Cyanodiazotrophs relative abundance" = "cyanoMM",
  "non_cyano" = "non_cyanoMM", "NCDs relative abundance" = "non_cyanoMM"
)

# Predictor set of the NCD model search (target = non_cyanoMM); Temp is both
# predictor and threshold variable.
NCD_PREDICTORS <- c("latitude", "PhotoBA", "Het", "DA", "DAtoBA", "DOCtoDON", "NH4", "P.PO4",
                    "N.NO3", "NtoP", "Conductivity", "DO", "pH", "GHI", "Temp", "Wdays",
                    "ShannonMM", "N2FixationRates", "N2FixationCell_fmol")

# Predictor set of the N2-fixation model search (target = N2FixationRates)
N2_PREDICTORS <- c("latitude", "PhotoBA", "Het", "DA", "DAtoBA", "DOCtoDON", "NH4", "P.PO4",
                   "N.NO3", "NtoP", "Conductivity", "DO", "pH", "GHI", "Temp", "Wdays",
                   "cyanoMM", "non_cyanoMM", "RichnessMM", "ShannonMM")

# column name -> descriptive label (first matching entry of VAR_MAP)
pretty_name <- function(var) {
  vapply(var, function(v) { k <- names(VAR_MAP)[VAR_MAP == v]; if (length(k)) k[length(k)] else v }, character(1))
}

to_var <- function(x) {
  v <- unname(VAR_MAP[x])
  v[is.na(v) & x %in% VAR_MAP] <- x[is.na(v) & x %in% VAR_MAP]   # already a column name (e.g. "NtoP")
  if (anyNA(v)) stop("No column mapping for: ", paste(x[is.na(v)], collapse = ", "),
                     ". Add it to VAR_MAP in R/supplement_io.R")
  v
}

# Full metadata table: one row per site, row names = site IDs
read_full_metadata <- function(path, sheet = 1) {
  if (!requireNamespace("readxl", quietly = TRUE)) stop("Install the 'readxl' package")
  d <- as.data.frame(readxl::read_excel(path, sheet = sheet))
  rownames(d) <- d$Site
  d
}

# Model-selection table (Fig 3 / Fig 4 supplementary sheets)
# columns: Temperature, Group, MLR predictors ("a+b+c"), RMSE, AdjR2, n, k, AIC [, p]
read_model_table <- function(path, sheet) {
  d <- as.data.frame(readxl::read_excel(path, sheet = sheet))
  stopifnot(all(c("Temperature", "Group", "MLR predictors", "RMSE", "AdjR2", "n", "k", "AIC") %in% names(d)))
  d <- d[!is.na(d$Temperature), ]
  d$ID    <- paste0(d$Group, d$Temperature)
  d$Temp  <- d$Temperature
  d$Group <- factor(d$Group, levels = c("Above", "Below"))
  d$RSS   <- d$RMSE^2 * d$n                       # derived
  d$Combination <- vapply(strsplit(d[["MLR predictors"]], "+", fixed = TRUE),
                          function(v) paste(to_var(v), collapse = "+"), character(1))
  d
}

# Variable-importance table (LMG, absolute shares of R2) from the supplement.
# Returns data.frame: Variable (descriptive), var (column name), Above, Below, plus totals.
read_importance <- function(path, sheet) {
  d <- as.data.frame(readxl::read_excel(path, sheet = sheet))
  names(d) <- c("Variable", "Above", "Below")
  tot <- d[d$Variable == "Total MLR R2", ]
  d <- d[d$Variable != "Total MLR R2", ]
  num <- function(x) suppressWarnings(as.numeric(x))        # "Not present" -> NA
  d$Above <- num(d$Above); d$Below <- num(d$Below)
  d$var <- to_var(d$Variable)
  attr(d, "total_R2") <- c(Above = num(tot$Above), Below = num(tot$Below))
  d
}

# Sites with a missing categorical value are dropped, not plotted as "NA"
drop_na_type <- function(d) d[!is.na(d$Type), ]
