#!/usr/bin/env Rscript
# coactivation_figures.R
# =============================================================================
# Figures 3 and S2: co-activation correlation heatmaps and companion panels.
#
# Accompanies:
#   Martinez M, Ozawa A, Van Zant D, Thornberry J, Toll L.
#   Whole Brain Cellular Activation After Mu and NOP Receptor Agonism
#   Identified Differential Regional and Network Consequences.
#   [JOURNAL] 2026. [PAPER_DOI]
#
# Code author:  Madeline Martinez
# Contact:      Lawrence Toll (corresponding author)
# Affiliation:  Toll Lab, Stiles-Nicholson Brain Institute,
#               Charles E. Schmidt College of Medicine,
#               Florida Atlantic University
# Repository:   https://github.com/toll-lab-code/cfos-opioid-brain-analyses
# Archived:     [ZENODO_DOI] (archive created at acceptance)
# License:      MIT (see LICENSE)
# Tested with:  R 4.4.3 (see sessionInfo.txt)
#
# Usage (from the repository root):
#   Rscript scripts/coactivation_figures.R
# Input and output locations are set in the Config section below; see README.
# =============================================================================
#
# Run network_analysis.R first: this script reads 01_Connectivity.xlsx, which
# that step produces. Heatmap style follows Kimbrough et al. 2020.
# =====================================================================
# Renders INDIVIDUAL pieces for hand-assembly in Illustrator:
#   * one heatmap per condition x grouping (Combined / Males / Females)
#     = square correlation matrix (jet colormap, -1..1) + a single
#       anatomy colour bar on the RIGHT (no title, no grid lines)
#   * two legends as separate files:
#       - colour key (anatomy divisions)
#       - correlation-strength colourbar (jet, horizontal)
#
# Combined  -> main Figure 3
# Males / Females -> supplemental figures
#
# All outputs saved as BOTH .tiff and .png into a "Figure 3" folder,
# matching the per-figure folder convention used by the other figure scripts.
# Base R graphics + optional showtext/Calibri.
# =====================================================================

# ---- Packages ----------------------------------------------------------------
required_pkgs <- c("readxl")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages. Please run:\n  install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}
# showtext and sysfonts are optional: if absent the script falls back to the
# device default font (see resolve_font() below).

# ---- Config ------------------------------------------------------------------
# Paths are relative to the repository root. Run from the repo root, e.g.:
#   Rscript scripts/coactivation_figures.R
# Override by editing the directories below or by setting the environment
# variables CFOS_OUTPUT_DIR / CFOS_SCRIPT_DIR before running.
DATA_DIR   <- path.expand(Sys.getenv("CFOS_DATA_DIR",   "data"))     # region_division_lookup.csv, read by divisions.R
OUTPUT_DIR <- path.expand(Sys.getenv("CFOS_OUTPUT_DIR", "results"))
SCRIPT_DIR <- path.expand(Sys.getenv("CFOS_SCRIPT_DIR", "scripts"))

# Panels are written to a per-figure subfolder, as in the other figure scripts. The
# subfolder name is set per figure set by apply_mode(); the value here is only
# a fallback for calling main() directly.
FIGURES_ROOT <- file.path(OUTPUT_DIR, "figures")
FIG_DIRNAME  <- "Figure3"

# ---- Figure sets -----------------------------------------------------
# Both figure sets are produced in a single run, each into its own folder:
#   "Figure3"  - Combined per-condition heatmaps plus the companion panels
#                (mean |r| by condition, region involvement, within-division
#                mean |r|), with condition titles baked into the panels.
#   "FigureS2" - the sex-stratified set (Males/Females per condition, condition
#                differences, sex differences) as small, titleless panels; row
#                and column labels are added at assembly rather than here.
# To render only one set, drop the other from RENDER_MODES. The size and panel
# overrides for each mode are applied by apply_mode() below.
RENDER_MODES <- c("Figure3", "FigureS2")
SUPP_MAT_MM  <- 30   # matrix side for the small supplement heatmaps

# Shared division scheme (13-division palette, abbreviations, region lookup) -
# identical across every figure so the anatomy grouping cannot drift.
source(file.path(SCRIPT_DIR, "divisions.R"))

# ---- Data source -----------------------------------------------------
# PRIMARY: the 12 correlation matrices are read from the network pipeline's
# 01_Connectivity.xlsx (sheets "Corr <trt> <sex>"), written by
# network_analysis.R.
CONNECTIVITY_XLSX <- file.path(OUTPUT_DIR, "01_Connectivity.xlsx")

# OPTIONAL fallback: a folder of "<stem>_correlation_matrix.csv" files
# (used only if the xlsx above is not found). Set to NULL to disable.
MATRICES_DIR <- NULL

# Sheet-name abbreviations used in 01_Connectivity.xlsx ("Corr <trt> <sex>").
TRT_ABBR <- c(Veh = "Vehicle", AcMor = "Acute Morphine",
              MorDep = "Morphine-Dependent", Ro = "Ro 64-6198")
SEX_ABBR <- c(C = "Combined", M = "Males", F = "Females")

# Condition stem (CSV fallback) -> display name + ordering.
CONDITION_DISPLAY <- c(
  "Vehicle" = "Vehicle",
  "Acute_Morphine" = "Acute Morphine", "Morphine" = "Acute Morphine",
  "Morphine-Dependent" = "Morphine-Dependent", "Chronic" = "Morphine-Dependent",
  "Ro" = "Ro 64-6198", "Ro_64-6198" = "Ro 64-6198"
)
CONDITION_ORDER <- c("Vehicle", "Acute Morphine", "Morphine-Dependent", "Ro 64-6198")
GROUPING_ORDER  <- c("Combined", "Males", "Females")

# ---- Sizes (mm) & resolution ----------------------------------------
# RENDER AT FINAL SIZE. Each panel is rendered at the size it will occupy
# in the assembled figure and placed at 100% in Illustrator (no rescaling),
# so every text element keeps its nominal point size on the printed page.
# NPP artwork spec: every text element >= 8 pt at final size, as in the other
# figure scripts. Set MAT_MM to that final size.
# Quick guide for a full-width (178 mm)
# Figure 3 with the two legends in a left column (~34 mm):
#   2x2 heatmaps  -> MAT_MM ~ 64  (TITLE_PT 8)   <- default
#   1x4 row       -> MAT_MM ~ 36  (TITLE_PT 6.5, tight -- prefer 2x2)
# For a 1.5-column (114 mm) 2x2: MAT_MM ~ 46 (TITLE_PT 7).
MAT_MM   <- 59      # final matrix side (square)
BAR_MM   <- 2.5     # right anatomy bar width
GAP_MM   <- 0.5     # gap between matrix and bar
BAKE_TITLES <- TRUE # condition name baked above each heatmap
TITLE_MM <- 6.0     # title strip height

# Apply Supplemental S2 overrides now that base sizes are defined
DPI      <- 1200    # line-art + raster; 1200 dpi for publication
MASK_DIAGONAL  <- FALSE   # TRUE -> blank the self-correlation diagonal
DRAW_BOUNDARIES<- FALSE   # TRUE -> thin black lines between divisions
LWD_FRAME <- 0.5

# ---- Type (points) ---------------------------------------------------
BASE_PS  <- 10
PT_TITLE     <- 10        # condition title above each heatmap (bold; bumped 8->10 pt to stand out vs the matrix; still NPP, not "large type")
PT_KEY_TITLE <- 8 
PT_KEY_ITEM  <- 8  
PT_CB_TITLE  <- 8 
PT_CB_TICK   <- 8
PT_CB_R      <- 8 
ptcex <- function(pt) pt / BASE_PS

# ====================================================================
# COMPANION PANELS - config (mean |r| line plot + dependence chord)
# ====================================================================
# Toggle which components this script generates:
DO_HEATMAPS <- TRUE
DO_DIFF     <- TRUE    # difference heatmaps: condition - condition (Delta r)
DO_SEXDIFF  <- TRUE    # sex-difference heatmaps: Male - Female, per condition
DO_MEANR    <- TRUE
DO_INVOLVE  <- TRUE    # ranked bar: regions by # dependence-strengthened edges
DO_DIVMAT   <- FALSE   # division x division matrix of strengthened edges (cut)
DO_DIVR     <- TRUE    # co-activation strength (mean |r|) by division x condition

# The values above define the main Figure 3 render. They are captured here so
# that apply_mode() can restore them when switching between figure sets.
MAIN_MAT_MM      <- MAT_MM
MAIN_BAKE_TITLES <- BAKE_TITLES
MAIN_TITLE_MM    <- TITLE_MM
MAIN_DO <- list(HEATMAPS = DO_HEATMAPS, DIFF = DO_DIFF, SEXDIFF = DO_SEXDIFF,
                MEANR    = DO_MEANR,    INVOLVE = DO_INVOLVE,
                DIVMAT   = DO_DIVMAT,   DIVR    = DO_DIVR)

# Configure the globals for one figure set. Called once per mode by the render
# loop at the end of the script. Supplement panels are small and titleless, and
# the Figure 3 companion panels (mean |r|, involvement, within-division) are
# skipped there because they are not part of S2.
apply_mode <- function(mode) {
  supplement  <- identical(mode, "FigureS2")
  FIG_DIRNAME <<- mode
  MAT_MM      <<- if (supplement) SUPP_MAT_MM else MAIN_MAT_MM
  BAKE_TITLES <<- if (supplement) FALSE       else MAIN_BAKE_TITLES
  TITLE_MM    <<- if (supplement) 0           else MAIN_TITLE_MM
  DO_HEATMAPS <<- MAIN_DO$HEATMAPS
  DO_DIFF     <<- MAIN_DO$DIFF
  DO_SEXDIFF  <<- MAIN_DO$SEXDIFF
  DO_DIVMAT   <<- MAIN_DO$DIVMAT
  DO_MEANR    <<- if (supplement) FALSE else MAIN_DO$MEANR
  DO_INVOLVE  <<- if (supplement) FALSE else MAIN_DO$INVOLVE
  DO_DIVR     <<- if (supplement) FALSE else MAIN_DO$DIVR
  CBAR_CORR_H <<- if (supplement) 44 else 65   # correlation colorbar height (mm): S2 matches the Dr bar (44) for consistency
  CBAR_DIFF_H <<- if (supplement) 44 else 65   # Dr colorbar height (mm): supplement bar sized to a single diff row
  invisible(mode)
}

# -- Panel 1: mean |r| by condition --
# By default this is computed from the correlation matrices (mean of
# |off-diagonal r| per condition x grouping); nothing is hardcoded. Set
# MEANR_FROM_DATA <- FALSE to use the reported fallback values instead. The
# computed table is printed so it can be checked against the reported values.
MEANR_FROM_DATA <- TRUE
MEANR_REPORTED <- data.frame(             # fallback / cross-check only
  Combined = c(0.45, 0.53, 0.74, 0.50),
  Male     = c(0.49, 0.53, 0.62, 0.52),
  Female   = c(0.55, 0.62, 0.83, 0.56),
  row.names = CONDITION_ORDER
)
COL_COMBINED <- "#222222"            # pooled (neutral dark)
COL_MALE     <- "#3f6fa3"            # muted blue  (matches Regional figs)
COL_FEMALE   <- "#7d2452"            # wine        (matches Regional figs)
COND_SHORT   <- c("Veh", "Mor", "MorDep", "Ro")   # x-axis tick labels (angled); abbreviations, matches Figs 2/3

# -- Panels 2 & 3: dependence-strengthened pairs (bar + division matrix) --
# Strengthened = pairs whose correlation is HIGHER in STRENGTH_HIGH than in
# STRENGTH_LOW. The loader finds the comparison in either orientation and
# fixes the sign automatically, so only the two conditions need naming.
STRENGTH_HIGH <- "Morphine-Dependent"   # condition with the increased coupling
STRENGTH_LOW  <- "Vehicle"              # baseline
STRENGTH_SEX  <- "Combined"
N_TOP         <- 20         # regions shown in the involvement bar
DIVMAT_METRIC <- "count"    # "count" = # strengthened edges; "sumdr" = summed Delta_r
STRENGTH_PALETTE <- colorRampPalette(
  c("#fff5eb","#fdd0a2","#fdae6b","#fd8d3c","#e6550d","#a63603"))(256)
DIVISION_ABBR <- c(division_abbr, "Other" = "Other")   # shared 13 abbrevs + Other fallback

# -- Panel 4: co-activation strength (mean |r|) by anatomical division x
# condition. Mean of |r| over within-division region pairs, one value per
# division per condition (Ardinger 2024 Fig 3 / Bloch 2024 precedent).
# Bars grouped by condition, coloured by division (decoded by the shared
# Color Key), drawn from a 0 baseline (field-standard). Same |r| metric as
# the whole-brain mean|r| panel, resolved by division.
DIVR_SEX <- "Combined"           # grouping to use: Combined / Males / Females

# -- Panel sizes (mm) & extra type sizes --
MR_W <- 88; MR_H <- 34          # mean|r| panel (bottom margin trimmed; taller plot)
IB_W <- 70; IB_H <- 87           # involvement bar (left margin trimmed)
DM_W <- 78; DM_H <- 72           # division x division matrix
DS_W <- 82; DS_H <- 52          # division-strength bar (bottom margin trimmed; taller plot)
LWD_AX  <- 0.6
PT_AXIS <- 8; PT_LAB <- 8; PT_LEG <- 8; PT_BAR <- 8; PT_CELL <- 8

# ---- Font ------------------------------------------------------------
# Arial throughout, matching every other figure script: NPP asks for a single
# sans-serif family (Arial / Helvetica) across all figures in a submission.
# Falls back to the device sans family (Helvetica on macOS) if Arial is absent.
FONT_NAME <- Sys.getenv("CFOS_FONT", "Arial")   # then Helvetica, then device sans
setup_font <- function(dpi = NA) {
  fam <- "sans"               # device sans = Helvetica on macOS (also approved)
  if (requireNamespace("showtext", quietly = TRUE) &&
      requireNamespace("sysfonts", quietly = TRUE)) {
    tryCatch({
      if (!(FONT_NAME %in% sysfonts::font_families())) {
        cands <- switch(FONT_NAME,
                        "Arial" = c("/Library/Fonts/Arial.ttf",
                                    "/System/Library/Fonts/Supplemental/Arial.ttf",
                                    "C:/Windows/Fonts/arial.ttf"),
                        "Helvetica" = c("/System/Library/Fonts/Helvetica.ttc",
                                        "/Library/Fonts/Helvetica.ttf"),
                        character(0))
        cands <- path.expand(cands); hit <- cands[file.exists(cands)]
        if (length(hit)) {
          b <- sub("(arial)\\.ttf$", "\\1bd.ttf", hit[1], ignore.case = TRUE)
          sysfonts::font_add(FONT_NAME, regular = hit[1],
                             bold = if (file.exists(b)) b else hit[1])
        }
      }
      if (FONT_NAME %in% sysfonts::font_families()) {
        showtext::showtext_auto()
        if (!is.na(dpi)) showtext::showtext_opts(dpi = dpi)
        fam <- FONT_NAME
      }
    }, error = function(e) invisible(NULL))
  }
  fam
}

# ---- Allen CCF groupings, colours, abbreviations ---------------------
ANATOMY_ORDER <- c(
  "Isocortex", "Olfactory", "Hippocampus",
  "Cortical Subplate", "Amygdala", "Striatum",
  "Pallidum/Septum", "Thalamus", "Hypothalamus",
  "Midbrain", "Pons", "Medulla",
  "Cerebellum"
)

ANATOMY_COLORS <- division_colors   # shared 13-division palette (Striatum = #ff5da2)

# Region -> division (Allen CCFv3 fine divisions; amygdalar nuclei kept as one
# Amygdala division and lateral septum grouped with Pallidum, matching the
# canonical region_division_lookup used by every other figure).
ALLEN_LOOKUP <- division_lookup   # shared canonical region -> division (region_division_lookup.csv)

# ---- Look: colormap + how r maps to colour ---------------------------
# Region order follows the field standard (Kimbrough 2020): Allen anatomical
# blocks, with regions reordered WITHIN each block by hierarchical clustering
# (see REORDER_* below). Colour preset:
#   "diverging" : blue-white-red, raw r, fixed -1..+1 (DEFAULT; the literature
#                 convention -- shows anticorrelations; Dependent stays red,
#                 which is the finding -- structure comes from the reordering).
#   "absolute"  : sequential white->red, scaled to the data range shared across
#                 the four conditions. Decompresses the warm end but clips the
#                 (few) anticorrelations to white. A labelled readability choice.
#   "relative"  : diverging, each matrix as deviation from grand-mean coupling.
#                 Maximises within-panel structure in every condition.
#   "jet"       : Kimbrough's literal rainbow, fixed -1..+1.
FIGURE_LOOK       <- "diverging"
SCALE_QUANTILES   <- c(0.02, 0.98)   # data range used by "absolute"
RELATIVE_HALFSPAN <- 0.40            # symmetric +/- range used by "relative"

# Within-block hierarchical reordering (Kimbrough 2020, Fig. 3): one region
# order is derived from a reference matrix and applied to ALL panels so they
# stay comparable. Values are NOT changed -- ordering only.
REORDER_WITHIN_GROUPS      <- TRUE
REORDER_REFERENCE_DISPLAY  <- "Morphine-Dependent"  # condition to derive order from
REORDER_REFERENCE_GROUPING <- "Combined"            # grouping whose matrix is used

.SEQ_ANCHORS <- c("#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a",
                  "#ef3b2c","#cb181d","#a50f15","#67000d")   # white -> dark red
.DIV_ANCHORS <- c("#2166ac","#4393c3","#92c5de","#f7f7f7","#f4a582",
                  "#d6604d","#b2182b")                        # blue -> white -> red
.JET_ANCHORS <- c("#0000FF","#00FFFF","#00FF00","#FFFF00","#FF0000")

LOOK <- switch(FIGURE_LOOK,
               diverging = list(anchors = .DIV_ANCHORS, mode = "fixed",  center = FALSE,
                                label = "Co-activation (Pearson r)"),
               absolute  = list(anchors = .SEQ_ANCHORS, mode = "shared", center = FALSE,
                                label = "Co-activation (Pearson r)"),
               relative  = list(anchors = .DIV_ANCHORS, mode = "center", center = TRUE,
                                label = "Co-activation (r \u2212 mean)"),
               jet       = list(anchors = .JET_ANCHORS, mode = "fixed",  center = FALSE,
                                label = "Correlation (Pearson r)"),
               stop("FIGURE_LOOK must be 'diverging', 'absolute', 'relative', or 'jet'"))
HEAT_PALETTE <- colorRampPalette(LOOK$anchors)(256)

# ---- Difference heatmaps (condition - condition) ---------------------
# Element-wise Delta r between two conditions, shown on a diverging
# blue-white-red scale centred at 0: red = stronger co-activation in the
# first condition, blue = weaker. Same region order and side-bars as the
# main matrices. A shared symmetric scale across all pairs keeps them
# comparable. Each pair is c(high, low) by DISPLAY name.
DIFF_PAIRS <- list(
  c("Acute Morphine",     "Vehicle"),
  c("Morphine-Dependent", "Vehicle"),
  c("Ro 64-6198",         "Vehicle"),
  c("Morphine-Dependent", "Acute Morphine"))   # escalation: red = gain under dependence
DIFF_SEX      <- "Combined"   # grouping to difference (Combined / Males / Females)
DIFF_LIMIT    <- 1            # fixed symmetric +/-1: shared Dr scale w/ sex-diffs, matches Fig 3's +/-1 correlation bar
DIFF_QUANTILE <- 0.98         # quantile of |Delta r| used when DIFF_LIMIT = "auto"
DIFF_PALETTE  <- colorRampPalette(.DIV_ANCHORS)(256)   # blue -> white -> red

# ---- Sex-difference heatmaps (Male - Female, per condition) ----------
# Within each condition, the difference between the male and female
# co-activation matrices. Red = stronger co-activation in the first sex.
# Own shared symmetric scale across the set (separate from the condition
# differences, since sex differences are typically smaller), own colourbar.
SEXDIFF_HIGH       <- "Males"     # first sex  (red where this sex is higher)
SEXDIFF_LOW        <- "Females"   # second sex
SEXDIFF_CONDITIONS <- CONDITION_ORDER   # which conditions (default: all four)
SEXDIFF_LIMIT      <- 1           # fixed +/-1: same shared Dr scale as condition-diffs -> one Dr colorbar for B and C
SEXDIFF_QUANTILE   <- 0.98

# Scale (lo, hi, center) for one grouping's matrices, per the chosen LOOK.
compute_scale <- function(mats) {
  offs <- unlist(lapply(mats, function(m) { d <- m; diag(d) <- NA; d[is.finite(d)] }))
  if (LOOK$mode == "fixed")  return(list(lo = -1, hi = 1, center = 0))
  if (LOOK$mode == "center") { gm <- mean(offs, na.rm = TRUE)
  return(list(lo = -RELATIVE_HALFSPAN, hi = RELATIVE_HALFSPAN, center = gm)) }
  q <- as.numeric(quantile(offs, SCALE_QUANTILES, na.rm = TRUE))   # "shared"
  list(lo = q[1], hi = q[2], center = 0)
}

# Aligned element-wise difference (high - low) for one grouping.
diff_matrix <- function(items, high, low, grouping) {
  hi_it <- Find(function(x) x$display == high && x$grouping == grouping, items)
  lo_it <- Find(function(x) x$display == low  && x$grouping == grouping, items)
  if (is.null(hi_it) || is.null(lo_it)) return(NULL)
  H <- as.matrix(hi_it$matrix); Lm <- as.matrix(lo_it$matrix)
  common <- intersect(rownames(H), rownames(Lm))
  if (length(common) < 2) return(NULL)
  H[common, common] - Lm[common, common]
}
# Aligned difference between two groupings (sexes) within one condition.
sex_diff_matrix <- function(items, cond, high, low) {
  hi_it <- Find(function(x) x$display == cond && x$grouping == high, items)
  lo_it <- Find(function(x) x$display == cond && x$grouping == low, items)
  if (is.null(hi_it) || is.null(lo_it)) return(NULL)
  H <- as.matrix(hi_it$matrix); Lm <- as.matrix(lo_it$matrix)
  common <- intersect(rownames(H), rownames(Lm))
  if (length(common) < 2) return(NULL)
  H[common, common] - Lm[common, common]
}
# Reference region order (within-block hclust) shared by all difference maps.
diff_region_order <- function(items) {
  if (!REORDER_WITHIN_GROUPS) return(NULL)
  refd <- Find(function(it) it$display == REORDER_REFERENCE_DISPLAY &&
                 it$grouping == REORDER_REFERENCE_GROUPING, items)
  if (is.null(refd)) refd <- items[[1]]
  compute_region_order(refd$matrix)
}
# Shared symmetric scale across a set of difference matrices.
diff_scale <- function(dmats, limit = DIFF_LIMIT, quantile_p = DIFF_QUANTILE) {
  if (is.numeric(limit)) {
    L <- abs(limit)
  } else {
    offs <- unlist(lapply(dmats, function(m) {
      d <- m; diag(d) <- NA; abs(d[is.finite(d)]) }))
    L <- as.numeric(quantile(offs, quantile_p, na.rm = TRUE))
  }
  L <- max(L, 1e-6)
  list(lo = -L, hi = L, center = 0)
}

# ---- Region ordering helpers -----------------------------------------
division_of <- function(regions) {
  d <- ALLEN_LOOKUP[regions]; d[is.na(d)] <- "Other"; unname(d)
}
sort_regions_by_anatomy <- function(regions) {
  cats <- division_of(regions)
  regions[order(factor(cats, levels = c(ANATOMY_ORDER, "Other")), regions)]
}
get_anatomy_colors_for_regions <- function(regions) {
  cols <- ANATOMY_COLORS[division_of(regions)]
  cols[is.na(cols)] <- "#999999"; unname(cols)
}
get_anatomy_blocks <- function(sorted_regions) {
  cats <- division_of(sorted_regions); r <- rle(cats)
  ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  data.frame(cat = r$values, start = starts, end = ends,
             stringsAsFactors = FALSE)
}

# Field-standard region order (Kimbrough 2020, Fig. 3): anatomical blocks in
# fixed order, regions reordered WITHIN each block by complete-linkage
# hierarchical clustering of their (Euclidean) correlation profiles. Derived
# from one reference matrix; values are never altered (ordering only).
compute_region_order <- function(ref_mat) {
  sr <- sort_regions_by_anatomy(rownames(ref_mat))     # anatomical base order
  if (!REORDER_WITHIN_GROUPS) return(sr)
  cm <- as.matrix(ref_mat[sr, sr])
  cm[!is.finite(cm)] <- 0; diag(cm) <- 1               # as in hclust_module_count
  blocks <- get_anatomy_blocks(sr)
  ordered <- character(0)
  for (k in seq_len(nrow(blocks))) {
    regs <- sr[blocks$start[k]:blocks$end[k]]
    if (length(regs) >= 3) {
      d  <- dist(cm[regs, , drop = FALSE], method = "euclidean")
      hc <- hclust(d, method = "complete")             # complete linkage, as in lineage
      regs <- regs[hc$order]
    }
    ordered <- c(ordered, regs)
  }
  ordered
}
# Choose black/white text for legibility on a coloured background.
contrast_text <- function(hex) {
  rgb <- col2rgb(hex) / 255
  lum <- 0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
  ifelse(lum > 0.55, "#000000", "#FFFFFF")
}

# ---- Renderers (draw into the current device) ------------------------
render_heatmap <- function(mat, title = NULL, family = "sans",
                           lo = -1, hi = 1, center = 0, order = NULL,
                           pal = HEAT_PALETTE) {
  if (is.null(order)) {
    sr <- sort_regions_by_anatomy(rownames(mat))
  } else {
    sr <- order[order %in% rownames(mat)]
    if (length(sr) != nrow(mat)) sr <- sort_regions_by_anatomy(rownames(mat))
  }
  m   <- as.matrix(mat[sr, sr]) - center
  n   <- length(sr)
  rc  <- get_anatomy_colors_for_regions(sr)
  if (MASK_DIAGONAL) diag(m) <- NA
  m[m < lo] <- lo; m[m > hi] <- hi          # clamp into the displayed range
  
  has_title <- BAKE_TITLES && !is.null(title)
  Wmm <- MAT_MM + GAP_MM + BAR_MM
  Hmm <- MAT_MM + (if (has_title) TITLE_MM else 0)
  nb  <- function(x0, x1, y0, y1) c(x0 / Wmm, x1 / Wmm, y0 / Hmm, y1 / Hmm)
  
  # matrix (occupies the lower MAT_MM of the canvas)
  z <- t(m)[, n:1, drop = FALSE]
  par(fig = nb(0, MAT_MM, 0, MAT_MM), new = TRUE, mar = c(0, 0, 0, 0))
  image(x = 0:n, y = 0:n, z = z, zlim = c(lo, hi), col = pal,
        axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
  if (DRAW_BOUNDARIES) {
    b <- get_anatomy_blocks(sr)
    for (k in seq_len(nrow(b))) if (b$start[k] > 1) {
      abline(h = n - b$start[k] + 1, col = "black", lwd = 0.3)
      abline(v = b$start[k] - 1,     col = "black", lwd = 0.3)
    }
  }
  box(lwd = LWD_FRAME)
  
  # right anatomy bar (region 1 at top), aligned to the matrix rows
  par(fig = nb(MAT_MM + GAP_MM, Wmm, 0, MAT_MM), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new(); plot.window(c(0, 1), c(0, n), xaxs = "i", yaxs = "i")
  rect(0, (n - 1):0, 1, n:1, col = rc, border = NA)
  rect(0, 0, 1, n, border = "black", lwd = LWD_FRAME)
  
  # title centred over the matrix, auto-shrunk if wider than the matrix
  if (has_title) {
    par(fig = nb(0, MAT_MM, MAT_MM, Hmm), new = TRUE, mar = c(0, 0, 0, 0))
    plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
    cx <- ptcex(PT_TITLE)
    w  <- strwidth(title, cex = cx, family = family)
    if (w > 0.96) cx <- cx * 0.96 / w
    text(0.5, 0.42, title, font = 2, cex = cx, family = family,
         adj = c(0.5, 0.5))
  }
}

render_colorkey <- function(family) {
  divs <- ANATOMY_ORDER; m <- length(divs)
  par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = FALSE); plot.new()
  plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  text(0.5, 0.985, "Color Key", font = 2, cex = ptcex(PT_KEY_TITLE),
       family = family, adj = c(0.5, 1))
  top <- 0.92; bot <- 0.00; gap <- 0.012   # 8pt title gets a gap above the swatches; Cerebellum bottom = matrix bottom (flush at assembly)
  rh  <- (top - bot - gap * (m - 1)) / m
  for (i in seq_len(m)) {
    y1 <- top - (i - 1) * (rh + gap); y0 <- y1 - rh
    col <- ANATOMY_COLORS[divs[i]]
    rect(0.06, y0, 0.94, y1, col = col, border = NA)
    text(0.5, (y0 + y1) / 2, divs[i], font = 2, cex = ptcex(PT_KEY_ITEM),
         family = family, col = contrast_text(col), adj = c(0.5, 0.5))
  }
}

render_corrbar <- function(family, lo = -1, hi = 1, label = "Correlation\nStrength",
                           pal = HEAT_PALETTE, sym = "r") {
  # VERTICAL colour bar: bar on the left, tick values to its right,
  # title rotated along the right edge, italic symbol below the bar.
  par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = FALSE); plot.new()
  plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  bx0 <- 0.22; bx1 <- 0.56; by0 <- 0.046; by1 <- 0.954   # bar = 0.908*65mm = 59mm (matches new matrix); ~3mm pad each end for +/-1 labels
  np  <- length(pal)
  ys  <- seq(by0, by1, length.out = np + 1)
  rect(bx0, ys[-(np + 1)], bx1, ys[-1], col = pal, border = NA)   # lo bottom -> hi top
  rect(bx0, by0, bx1, by1, border = "black", lwd = LWD_FRAME)
  ticks <- pretty(c(lo, hi), n = 4); ticks <- ticks[ticks >= lo & ticks <= hi]
  ty <- by0 + (ticks - lo) / (hi - lo) * (by1 - by0)
  segments(bx1, ty, bx1 + 0.07, ty, lwd = LWD_FRAME)
  text(bx1 + 0.12, ty, formatC(ticks, format = "g"),
       cex = ptcex(PT_CB_TICK), family = family, adj = c(0, 0.5))
  ttl <- if (grepl("Pearson", label))
    expression(bold("Co-activation (Pearson ") * bolditalic(r) * bold(")"))
  else if (grepl("\u0394", sym))
    expression(bold("Co-activation (") * bold(Delta) * bolditalic(r) * bold(")"))
  else gsub("\n", " ", label)
  text(0.10, 0.5, ttl, font = 2, cex = ptcex(PT_CB_TITLE),
       family = family, srt = 90, adj = c(0.5, 0.5))
}

# ---- Save one renderer to both .tiff and .png ------------------------
save_both <- function(render_fn, path_noext, w_mm, h_mm) {
  fam <- setup_font(dpi = DPI)
  for (ext in c("tiff", "png")) {
    f <- paste0(path_noext, ".", ext)
    if (ext == "tiff") {
      tiff(f, width = w_mm, height = h_mm, units = "mm", res = DPI,
           pointsize = BASE_PS, compression = "lzw", type = "cairo", bg = "white")
    } else {
      png(f, width = w_mm, height = h_mm, units = "mm", res = DPI,
          pointsize = BASE_PS, type = "cairo", bg = "white")
    }
    par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0)); plot.new()
    render_fn(fam)
    dev.off()
  }
  message("  ", basename(path_noext), " (.tiff/.png)")
}

# ---- Data loading ----------------------------------------------------
# Parse "<condition>_<Grouping>" stems; returns list of
# list(grouping, cond_key, display, matrix).
discover_csv_matrices <- function(matrices_dir) {
  if (!dir.exists(matrices_dir)) stop("Matrices dir not found: ", matrices_dir)
  files <- list.files(matrices_dir, pattern = "_correlation_matrix\\.csv$")
  out <- list()
  for (f in files) {
    stem <- sub("_correlation_matrix\\.csv$", "", f)
    g <- NA
    for (grp in GROUPING_ORDER)               # case-sensitive: "Males" != "Females"
      if (grepl(grp, stem, fixed = TRUE)) { g <- grp; break }
    if (is.na(g)) { message("  (skip, no grouping token) ", f); next }
    cond_key <- gsub(g, "", stem, fixed = TRUE)          # remove token anywhere
    cond_key <- trimws(gsub("(^[_ -]+)|([_ -]+$)", "", cond_key))  # tidy seps
    disp <- CONDITION_DISPLAY[cond_key]; if (is.na(disp)) disp <- gsub("_", " ", cond_key)
    df <- read.csv(file.path(matrices_dir, f), row.names = 1, check.names = FALSE)
    rownames(df) <- trimws(rownames(df)); colnames(df) <- trimws(colnames(df))
    out[[length(out) + 1]] <- list(grouping = g, cond_key = cond_key,
                                   display = unname(disp), matrix = as.matrix(df))
  }
  out
}

# Parse a "Corr <trt> <sex>" sheet name -> grouping / cond_key / display.
parse_corr_sheet <- function(sheet) {
  s <- sub("^Corr\\s+", "", sheet)
  toks <- strsplit(trimws(s), "\\s+")[[1]]
  if (length(toks) < 2) return(NULL)
  sex_tok <- toks[length(toks)]
  trt_tok <- paste(toks[-length(toks)], collapse = " ")
  grp <- SEX_ABBR[sex_tok]; if (is.na(grp)) grp <- sex_tok
  dsp <- TRT_ABBR[trt_tok]; if (is.na(dsp)) dsp <- trt_tok
  list(grouping = unname(grp), cond_key = trt_tok, display = unname(dsp))
}

# PRIMARY loader: auto-discover every "Corr ..." sheet in 01_Connectivity.xlsx.
load_connectivity_xlsx <- function(xlsx_path) {
  if (!requireNamespace("readxl", quietly = TRUE))
    stop("Package 'readxl' is required to read ", basename(xlsx_path),
         ".\n  Run: install.packages('readxl')")
  sheets <- readxl::excel_sheets(xlsx_path)
  corr   <- sheets[grepl("^Corr\\b", sheets)]
  if (length(corr) == 0)
    stop("No 'Corr ...' sheets found in ", xlsx_path)
  out <- list()
  for (sh in corr) {
    meta <- parse_corr_sheet(sh); if (is.null(meta)) next
    raw <- as.data.frame(readxl::read_excel(xlsx_path, sheet = sh))
    rn  <- as.character(raw[[1]])                      # first col = region names
    mat <- as.matrix(raw[, -1, drop = FALSE]); storage.mode(mat) <- "double"
    rownames(mat) <- trimws(rn); colnames(mat) <- trimws(colnames(mat))
    out[[length(out) + 1]] <- c(meta, list(matrix = mat))
    message(sprintf("  %-14s -> %-18s (%s)", sh, meta$display, meta$grouping))
  }
  out
}

slug <- function(s) gsub("[^A-Za-z0-9]+", "_", s)

# =====================================================================
# COMPANION PANEL RENDERERS
# =====================================================================
# ---- Panel 1: mean |r| by condition ---------------------------------
# ---- Panel 1: mean |r| by condition ---------------------------------
# Mean of |off-diagonal r| for one matrix.
mean_abs_r <- function(mat) {
  d <- as.matrix(mat); diag(d) <- NA
  mean(abs(d[is.finite(d)]))
}
# Build the Combined/Male/Female x condition table from loaded Corr matrices.
compute_meanr <- function(items) {
  gmap <- c(Combined = "Combined", Male = "Males", Female = "Females")
  tab <- data.frame(row.names = CONDITION_ORDER)
  for (col in names(gmap)) {
    vals <- vapply(CONDITION_ORDER, function(cond) {
      it <- Find(function(x) x$display == cond && x$grouping == gmap[[col]], items)
      if (is.null(it)) NA_real_ else mean_abs_r(it$matrix)
    }, numeric(1))
    tab[[col]] <- round(unname(vals), 3)
  }
  tab
}

render_meanr <- function(family) {
  y <- .MEANR
  rng  <- range(unlist(y), na.rm = TRUE)
  lo   <- floor((rng[1] - 0.01) / 0.1) * 0.1
  hi   <- ceiling((rng[2] + 0.01) / 0.1) * 0.1
  yt   <- seq(lo, hi, 0.1)
  ylim <- c(lo, hi)
  xL <- 0.70; xR <- 4.30
  par(family = family, ps = BASE_PS,
      mar = c(2.1, 3.0, 0.8, 2.4), mgp = c(2.1, 0.5, 0), tcl = -0.25)
  plot(NA, xlim = c(xL, xR), ylim = ylim, axes = FALSE, xlab = "", ylab = "")
  # clean L-shaped axes that meet at the bottom-left corner
  segments(xL, ylim[1], xR, ylim[1], lwd = LWD_AX)            # x baseline
  segments(xL, ylim[1], xL, max(yt),  lwd = LWD_AX)           # y axis
  axis(2, at = yt, pos = xL, lwd = 0, lwd.ticks = LWD_AX,
       cex.axis = ptcex(PT_AXIS), las = 1)
  segments(1:4, ylim[1], 1:4, ylim[1] - 0.008, lwd = LWD_AX, xpd = NA)  # x ticks
  mtext(expression("Mean |" * italic(r) * "|"), side = 2, line = 1.2, cex = ptcex(PT_LAB))
  text(1:4, ylim[1] - 0.045, COND_SHORT, srt = 30, adj = 1, xpd = NA,
       cex = ptcex(PT_AXIS))
  lines(1:4, y$Male,    col = COL_MALE,   lty = 2, lwd = 1.3)
  points(1:4, y$Male,   col = COL_MALE,   pch = 1, cex = 0.8, lwd = 1.2)
  lines(1:4, y$Female,  col = COL_FEMALE, lty = 3, lwd = 1.3)
  points(1:4, y$Female, col = COL_FEMALE, pch = 2, cex = 0.8, lwd = 1.2)
  lines(1:4, y$Combined,  col = COL_COMBINED, lwd = 2.4)
  points(1:4, y$Combined, col = COL_COMBINED, pch = 19, cex = 1.0)
  # manual legend: top-right, first row dead-level with the 0.9 tick, pushed right
  .ll <- c("Combined","Male","Female")
  .lc <- c(COL_COMBINED, COL_MALE, COL_FEMALE)
  .lt <- c(1,2,3); .lp <- c(19,1,2); .lw <- c(2.4,1.3,1.3)
  .kL <- xR - 0.40; .kR <- xR - 0.12; .tx <- xR - 0.06   # sits at the plot edge / start of margin; +/- to all three shifts it right/left
  .y0 <- hi; .dy <- 0.08                                  # .y0 = 0.9 (dead level); row pitch widened for 8pt
  for (i in 1:3) {
    yy <- .y0 - (i - 1) * .dy
    segments(.kL, yy, .kR, yy, col = .lc[i], lty = .lt[i], lwd = .lw[i], xpd = NA)
    points((.kL + .kR) / 2, yy, col = .lc[i], pch = .lp[i], cex = 0.8, lwd = 1.2, xpd = NA)
    text(.tx, yy, .ll[i], adj = c(0, 0.5), cex = ptcex(PT_LEG), xpd = NA)
  }
}

# ---- Panel 4: co-activation strength by division x condition --------
# Mean |r| over within-division region pairs for one matrix.
div_strength_vec <- function(mat) {
  regs <- rownames(mat); dv <- division_of(regs)
  d <- as.matrix(mat); diag(d) <- NA
  out <- setNames(rep(NA_real_, length(ANATOMY_ORDER)), ANATOMY_ORDER)
  for (dn in ANATOMY_ORDER) {
    idx <- which(dv == dn)
    # Need >=3 regions (>=3 within-division pairs) for a stable mean; this
    # excludes Cortical Subplate (2 regions = 1 pair) from the strength panels
    # only. It still appears, coloured, in the heatmap anatomy key.
    if (length(idx) < 3) next
    sub <- d[idx, idx, drop = FALSE]
    vals <- abs(sub[upper.tri(sub)])
    out[dn] <- if (any(is.finite(vals))) mean(vals, na.rm = TRUE) else NA_real_
  }
  out
}
# divisions (rows) x conditions (cols) matrix of within-division mean |r|.
compute_div_strength <- function(items, grouping) {
  M <- matrix(NA_real_, length(ANATOMY_ORDER), length(CONDITION_ORDER),
              dimnames = list(ANATOMY_ORDER, CONDITION_ORDER))
  for (cond in CONDITION_ORDER) {
    it <- Find(function(x) x$display == cond && x$grouping == grouping, items)
    if (!is.null(it)) M[, cond] <- div_strength_vec(it$matrix)
  }
  M
}

render_divr <- function(family) {
  M <- .DIVR
  M <- M[rowSums(is.finite(M)) > 0, , drop = FALSE]   # drop divisions w/ no within-division estimate (Cortical Subplate: 2 regions) -> NOTE in caption
  divs <- rownames(M); conds <- colnames(M)
  nd <- length(divs); nc <- length(conds)
  cols <- ANATOMY_COLORS[divs]
  ymax <- max(M, na.rm = TRUE)
  ylim <- c(0, ymax * 1.06)
  bw <- 1; gap_in <- 0.08; gap_grp <- 5.0           # bar / intra-group / group gaps
  grp_w <- nd * bw + (nd - 1) * gap_in
  xstep <- grp_w + gap_grp
  xmax  <- nc * grp_w + (nc - 1) * gap_grp
  par(family = family, ps = BASE_PS, mar = c(2.1, 3.0, 0.7, 0.6),
      mgp = c(2.1, 0.5, 0), tcl = -0.25)
  plot(NA, xlim = c(0, xmax), ylim = ylim, axes = FALSE, xlab = "", ylab = "")
  yt <- pretty(c(0, ymax), n = 5); yt <- yt[yt <= ylim[2]]
  segments(0, 0, xmax, 0, lwd = LWD_AX)                       # x baseline
  segments(0, 0, 0, max(yt), lwd = LWD_AX)                    # y axis
  axis(2, at = yt, pos = 0, lwd = 0, lwd.ticks = LWD_AX,
       cex.axis = ptcex(PT_AXIS), las = 1)
  mtext(expression("Mean |" * italic(r) * "| (within division)"), side = 2, line = 1.2, cex = ptcex(PT_LAB))
  for (j in seq_len(nc)) {
    x0 <- (j - 1) * xstep
    for (i in seq_len(nd)) {
      v <- M[i, j]; if (!is.finite(v)) next
      xl <- x0 + (i - 1) * (bw + gap_in)
      rect(xl, 0, xl + bw, v, col = cols[i], border = NA)
    }
    text(x0 + grp_w / 2, -ymax * 0.045, COND_SHORT[j], srt = 30,
         adj = c(1, 1), xpd = NA, cex = ptcex(PT_AXIS))
  }
}

# ---- Panel 2: chord of dependence-strengthened pairs ----------------
load_strengthened <- function(xlsx, high, low, sex) {
  if (!requireNamespace("readxl", quietly = TRUE))
    stop("readxl required to read ", basename(xlsx))
  fz <- as.data.frame(readxl::read_excel(xlsx, sheet = "FisherZ_SigPairs"))
  message("  Comparisons available: ", paste(unique(fz$Comparison), collapse = " | "))
  lab1 <- paste(high, "vs", low); lab2 <- paste(low, "vs", high)
  d <- fz[fz$Sex == sex & fz$Comparison %in% c(lab1, lab2), ]
  if (nrow(d) == 0) {
    message("  no rows for '", lab1, "' or '", lab2, "' (Sex=", sex,
            ") -- check STRENGTH_HIGH/LOW/SEX against the list above.")
    return(data.frame(Region_A = character(0), Region_B = character(0),
                      Delta_r = numeric(0)))
  }
  # orient so positive = increase in `high`; keep only strengthened pairs
  inc <- ifelse(d$Comparison == lab2, -d$Delta_r, d$Delta_r)
  d <- d[inc > 0, , drop = FALSE]; inc <- inc[inc > 0]
  data.frame(Region_A = trimws(d$Region_A), Region_B = trimws(d$Region_B),
             Delta_r = inc, stringsAsFactors = FALSE)
}
# Region involvement: # strengthened edges each region participates in.
region_involvement <- function(edges) {
  tb <- sort(table(c(edges$Region_A, edges$Region_B)), decreasing = TRUE)
  data.frame(region = names(tb), n = as.integer(tb), stringsAsFactors = FALSE)
}
# Division x division matrix of strengthened edges (count or summed Delta_r).
build_divmat <- function(edges, metric = "count") {
  da <- division_of(edges$Region_A); db <- division_of(edges$Region_B)
  present <- unique(c(da, db))
  divs <- c(ANATOMY_ORDER[ANATOMY_ORDER %in% present],
            if ("Other" %in% present) "Other")
  M <- matrix(0, length(divs), length(divs), dimnames = list(divs, divs))
  for (k in seq_len(nrow(edges))) {
    i <- da[k]; j <- db[k]
    v <- if (metric == "sumdr") edges$Delta_r[k] else 1
    M[i, j] <- M[i, j] + v
    if (i != j) M[j, i] <- M[j, i] + v
  }
  M
}

# ---- Panel 2: ranked involvement bar --------------------------------
render_involvement <- function(family) {
  d <- .INVO[order(.INVO$n), ]                 # ascending -> largest at top
  cols <- ANATOMY_COLORS[division_of(d$region)]; cols[is.na(cols)] <- "#999999"
  xmax <- max(d$n)
  par(family = family, ps = BASE_PS, mar = c(3.4, 2.7, 0.2, 1.0),
      mgp = c(2.0, 0.5, 0), yaxs = "i")
  bp <- barplot(d$n, horiz = TRUE, col = cols, border = NA, names.arg = d$region,
                las = 1, cex.names = ptcex(PT_BAR), axes = FALSE,
                ylim = c(0.2, 1.2 * length(d$n)),
                xlim = c(0, xmax * 1.10), xpd = FALSE)
  axis(1, lwd = LWD_AX, cex.axis = ptcex(PT_AXIS))
  mtext("Dependence-strengthened edges", side = 1, line = 1.45, cex = ptcex(PT_LAB))
  text(d$n + xmax * 0.015, bp, d$n, adj = c(0, 0.5),
       cex = ptcex(PT_BAR), family = family, xpd = NA)
  # division colours are decoded by the shared Fig3_legend_ColorKey
}

# ---- Panel 3: division x division matrix ----------------------------
render_divmat <- function(family) {
  M <- .DIVMAT; n <- nrow(M); divs <- rownames(M)
  ab <- DIVISION_ABBR[divs]; ab[is.na(ab)] <- divs
  lab <- if (DIVMAT_METRIC == "sumdr") "\u03a3 \u0394r" else "Edges"
  vmax <- max(M); if (vmax <= 0) vmax <- 1
  par(family = family, ps = BASE_PS, mar = c(0.4, 3.0, 3.0, 4.2))
  plot.new(); plot.window(xlim = c(0, n), ylim = c(0, n), asp = 1,
                          xaxs = "i", yaxs = "i")
  for (i in seq_len(n)) for (j in seq_len(n)) {
    v <- M[i, j]
    col <- if (v <= 0) "#ffffff" else STRENGTH_PALETTE[max(1, round(v / vmax * 255) + 1)]
    rect(j - 1, n - i, j, n - i + 1, col = col, border = "white", lwd = 0.5)
    if (v > 0)
      text(j - 0.5, n - i + 0.5, if (DIVMAT_METRIC == "sumdr") sprintf("%.1f", v) else v,
           cex = ptcex(PT_CELL), family = family,
           col = if (v > 0.62 * vmax) "white" else "black")
  }
  rect(0, 0, n, n, border = "black", lwd = LWD_FRAME)
  for (k in seq_len(n)) {
    text(k - 0.5, n + 0.12, ab[k], srt = 45, adj = c(0, 0.5), xpd = NA,
         cex = ptcex(PT_AXIS), col = ANATOMY_COLORS[divs[k]], family = family)
    text(-0.12, n - k + 0.5, ab[k], adj = c(1, 0.5), xpd = NA,
         cex = ptcex(PT_AXIS), col = ANATOMY_COLORS[divs[k]], family = family)
  }
  # compact value colourbar on the right
  npl <- length(STRENGTH_PALETTE)
  bx0 <- n + 0.30; bx1 <- n + 0.70; ys <- seq(0, n, length.out = npl + 1)
  rect(bx0, ys[-(npl + 1)], bx1, ys[-1], col = STRENGTH_PALETTE, border = NA, xpd = NA)
  rect(bx0, 0, bx1, n, border = "black", lwd = 0.5, xpd = NA)
  text(bx1 + 0.08, c(0, n), c("0", if (DIVMAT_METRIC == "sumdr") sprintf("%.0f", vmax) else vmax),
       adj = c(0, 0.5), cex = ptcex(PT_AXIS), family = family, xpd = NA)
  text((bx0 + bx1) / 2, n + 0.62, lab, adj = c(0.5, 0), cex = ptcex(PT_LAB),
       family = family, xpd = NA)
}


# ---- Main ------------------------------------------------------------
main <- function() {
  fig_dir <- file.path(FIGURES_ROOT, FIG_DIRNAME)
  if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
  
  # ---------- Load Corr matrices once (heatmaps + computed mean|r|) ----------
  items <- NULL
  if (DO_HEATMAPS || DO_DIFF || DO_SEXDIFF || (DO_MEANR && MEANR_FROM_DATA) || DO_DIVR) {
    if (!is.null(CONNECTIVITY_XLSX) && file.exists(CONNECTIVITY_XLSX)) {
      message("Reading matrices from: ", CONNECTIVITY_XLSX)
      items <- load_connectivity_xlsx(CONNECTIVITY_XLSX)
    } else if (!is.null(MATRICES_DIR) && dir.exists(MATRICES_DIR)) {
      message("Reading matrices from CSV folder: ", MATRICES_DIR)
      items <- discover_csv_matrices(MATRICES_DIR)
    } else if (DO_HEATMAPS) {
      stop("No data source found. Set one of these at the top of the script:\n",
           "  CONNECTIVITY_XLSX -> full path to your 01_Connectivity.xlsx\n",
           "      (tried: ", CONNECTIVITY_XLSX, ")\n",
           "  MATRICES_DIR      -> folder of *_correlation_matrix.csv files\n",
           "      (", if (is.null(MATRICES_DIR)) "disabled" else MATRICES_DIR, ")")
    }
    if (!is.null(items)) {
      if (length(items) == 0) stop("Data source found but no matrices loaded.")
      ord <- order(match(vapply(items, `[[`, "", "grouping"), GROUPING_ORDER),
                   match(vapply(items, `[[`, "", "display"), CONDITION_ORDER))
      items <- items[ord]
    }
  }
  
  # ---------- Heatmaps + colourbars + colour key ----------
  if (DO_HEATMAPS) {
    if (is.null(items)) stop("Heatmaps need a data source (CONNECTIVITY_XLSX/MATRICES_DIR).")
    Wmm <- MAT_MM + GAP_MM + BAR_MM
    Hmm <- MAT_MM + (if (BAKE_TITLES) TITLE_MM else 0)
    groupings <- intersect(GROUPING_ORDER,
                           unique(vapply(items, `[[`, "", "grouping")))
    
    region_order <- NULL
    if (REORDER_WITHIN_GROUPS) {
      ref <- Find(function(it) it$display == REORDER_REFERENCE_DISPLAY &&
                    it$grouping == REORDER_REFERENCE_GROUPING, items)
      if (is.null(ref)) ref <- items[[1]]
      region_order <- compute_region_order(ref$matrix)
      message("Region order: within-block hclust from ", ref$display,
              " (", ref$grouping, ")")
    }
    
    message("Heatmaps -> ", fig_dir, "   (look: ", FIGURE_LOOK, ")")
    for (g in groupings) {
      g_items <- Filter(function(it) it$grouping == g, items)
      sc <- compute_scale(lapply(g_items, `[[`, "matrix"))
      message(sprintf("  %-9s scale [%.2f, %.2f]%s", g, sc$lo, sc$hi,
                      if (sc$center != 0) sprintf(" centered on %.2f", sc$center) else ""))
      for (it in g_items) {
        base <- file.path(fig_dir, sprintf("Fig3_%s_%s_heatmap",
                                           it$grouping, slug(it$display)))
        local({
          m <- it$matrix; ttl <- it$display; s <- sc; ordr <- region_order
          save_both(function(fam) render_heatmap(m, ttl, fam, s$lo, s$hi,
                                                 s$center, ordr),
                    base, Wmm, Hmm)
        })
      }
      local({
        s <- sc
        save_both(function(fam) render_corrbar(fam, s$lo, s$hi, LOOK$label),
                  file.path(fig_dir, paste0("Fig3_legend_CorrelationStrength_", g)),
                  18, CBAR_CORR_H)
      })
    }
    save_both(render_colorkey, file.path(fig_dir, "Fig3_legend_ColorKey"), 26, 59)
    message("  ", length(items), " heatmaps + ", length(groupings),
            " colourbars + 1 colour key")
  }
  
  # ---------- Difference heatmaps (condition - condition) ----------
  if (DO_DIFF) {
    if (is.null(items)) {
      message("Difference heatmaps skipped: need Corr matrices (CONNECTIVITY_XLSX).")
    } else {
      ord_d <- diff_region_order(items)
      dmats <- list(); titles <- character(0)
      for (pr in DIFF_PAIRS) {
        dm <- diff_matrix(items, pr[1], pr[2], DIFF_SEX)
        if (is.null(dm)) {
          message("  diff skipped (missing matrix): ", pr[1], " \u2212 ", pr[2],
                  " [", DIFF_SEX, "]")
          next
        }
        dmats[[length(dmats) + 1]] <- dm
        titles <- c(titles, paste0(pr[1], " \u2212 ", pr[2]))
      }
      if (length(dmats) == 0) {
        message("Difference heatmaps: no pairs available for grouping '",
                DIFF_SEX, "'. Check DIFF_PAIRS / DIFF_SEX.")
      } else {
        sc <- diff_scale(dmats)
        message(sprintf("Difference heatmaps (%s) -> %s   scale [%.2f, %.2f] (%s)",
                        DIFF_SEX, fig_dir, sc$lo, sc$hi,
                        if (is.numeric(DIFF_LIMIT)) "fixed" else
                          sprintf("auto q%.2f", DIFF_QUANTILE)))
        Wmm <- MAT_MM + GAP_MM + BAR_MM
        Hmm <- MAT_MM + (if (BAKE_TITLES) TITLE_MM else 0)
        for (i in seq_along(dmats)) {
          local({
            m <- dmats[[i]]; ttl <- titles[i]; s <- sc; ordr <- ord_d
            base <- file.path(fig_dir, sprintf("Fig3_Diff_%s_heatmap", slug(ttl)))
            save_both(function(fam) render_heatmap(m, ttl, fam, s$lo, s$hi,
                                                   s$center, ordr, pal = DIFF_PALETTE), base, Wmm, Hmm)
          })
        }
        save_both(function(fam) render_corrbar(fam, sc$lo, sc$hi,
                                               "Co-activation change", pal = DIFF_PALETTE, sym = "\u0394r"),
                  file.path(fig_dir, "Fig3_legend_DiffStrength"), 18, CBAR_DIFF_H)
        message("  ", length(dmats), " difference heatmaps + 1 colourbar")
      }
    }
  }
  
  # ---------- Sex-difference heatmaps (Male - Female, per condition) ----------
  if (DO_SEXDIFF) {
    if (is.null(items)) {
      message("Sex-difference heatmaps skipped: need Corr matrices (CONNECTIVITY_XLSX).")
    } else {
      sex_lab <- c(Males = "Male", Females = "Female", Combined = "Combined")
      hlab <- if (!is.na(sex_lab[SEXDIFF_HIGH])) sex_lab[SEXDIFF_HIGH] else SEXDIFF_HIGH
      llab <- if (!is.na(sex_lab[SEXDIFF_LOW]))  sex_lab[SEXDIFF_LOW]  else SEXDIFF_LOW
      ord_s <- diff_region_order(items)
      sdmats <- list(); sconds <- character(0)
      for (cond in SEXDIFF_CONDITIONS) {
        dm <- sex_diff_matrix(items, cond, SEXDIFF_HIGH, SEXDIFF_LOW)
        if (is.null(dm)) {
          message("  sex-diff skipped (missing ", SEXDIFF_HIGH, "/", SEXDIFF_LOW,
                  " matrix): ", cond)
          next
        }
        sdmats[[length(sdmats) + 1]] <- dm
        sconds <- c(sconds, cond)
      }
      if (length(sdmats) == 0) {
        message("Sex-difference heatmaps: no conditions have both ",
                SEXDIFF_HIGH, " and ", SEXDIFF_LOW, " matrices.")
      } else {
        sc <- diff_scale(sdmats, SEXDIFF_LIMIT, SEXDIFF_QUANTILE)
        message(sprintf("Sex-difference heatmaps (%s \u2212 %s) -> %s   scale [%.2f, %.2f] (%s)",
                        hlab, llab, fig_dir, sc$lo, sc$hi,
                        if (is.numeric(SEXDIFF_LIMIT)) "fixed" else
                          sprintf("auto q%.2f", SEXDIFF_QUANTILE)))
        Wmm <- MAT_MM + GAP_MM + BAR_MM
        Hmm <- MAT_MM + (if (BAKE_TITLES) TITLE_MM else 0)
        for (i in seq_along(sdmats)) {
          local({
            m <- sdmats[[i]]; ttl <- sprintf("%s (%s \u2212 %s)", sconds[i], hlab, llab)
            s <- sc; ordr <- ord_s
            base <- file.path(fig_dir, sprintf("Fig3_SexDiff_%s_heatmap", slug(sconds[i])))
            save_both(function(fam) render_heatmap(m, ttl, fam, s$lo, s$hi,
                                                   s$center, ordr, pal = DIFF_PALETTE), base, Wmm, Hmm)
          })
        }
        save_both(function(fam) render_corrbar(fam, sc$lo, sc$hi,
                                               sprintf("%s \u2212 %s", hlab, llab),
                                               pal = DIFF_PALETTE, sym = "\u0394r"),
                  file.path(fig_dir, "Fig3_legend_SexDiff"), 18, CBAR_DIFF_H)
        message("  ", length(sdmats), " sex-difference heatmaps + 1 colourbar")
      }
    }
  }
  
  # ---------- Panel 1: mean |r| by condition ----------
  if (DO_MEANR) {
    if (MEANR_FROM_DATA && !is.null(items)) {
      tab <- compute_meanr(items)
      message("Mean |r| computed from Corr matrices:")
      print(tab)
    } else {
      if (MEANR_FROM_DATA) message("Mean |r|: no matrices loaded; using reported values.")
      tab <- MEANR_REPORTED
    }
    .MEANR <<- tab
    message("Mean |r| panel -> ", fig_dir)
    save_both(render_meanr, file.path(fig_dir, "Fig3_MeanAbsR_byCondition"), MR_W, MR_H)
  }
  
  # ---------- Panel 4: co-activation strength by division x condition ----------
  if (DO_DIVR) {
    if (is.null(items)) {
      message("Division-strength panel skipped: needs Corr matrices (CONNECTIVITY_XLSX).")
    } else {
      M <- compute_div_strength(items, DIVR_SEX)
      if (all(is.na(M))) {
        message("Division-strength: no '", DIVR_SEX, "' Corr matrices found; ",
                "set DIVR_SEX to Combined / Males / Females.")
      } else {
        .DIVR <<- M
        message("Co-activation strength by division (", DIVR_SEX, "), mean |r|:")
        print(round(M, 3))
        save_both(render_divr,
                  file.path(fig_dir, "Fig3_Strength_byDivision"), DS_W, DS_H)
      }
    }
  }
  
  # ---------- Panels 2 & 3: dependence-strengthened pairs ----------
  if (DO_INVOLVE || DO_DIVMAT) {
    if (is.null(CONNECTIVITY_XLSX) || !file.exists(CONNECTIVITY_XLSX)) {
      message("Strengthened-pair panels skipped: need CONNECTIVITY_XLSX ",
              "(FisherZ_SigPairs sheet).")
    } else {
      message("Strengthened-pair panels -> ", fig_dir)
      edges <- load_strengthened(CONNECTIVITY_XLSX, STRENGTH_HIGH, STRENGTH_LOW,
                                 STRENGTH_SEX)
      message(sprintf("  strengthened pairs: %d", nrow(edges)))
      if (nrow(edges) < 1) {
        message("  none found -- check STRENGTH_HIGH / STRENGTH_LOW / ",
                "STRENGTH_SEX against the comparison list above.")
      } else {
        if (DO_INVOLVE) {
          invo <- region_involvement(edges)
          .INVO <<- head(invo, N_TOP)
          message(sprintf("  involvement bar: top %d of %d involved regions",
                          nrow(.INVO), nrow(invo)))
          save_both(render_involvement,
                    file.path(fig_dir, "Fig3_Involvement_byRegion"), IB_W, IB_H)
        }
        if (DO_DIVMAT) {
          .DIVMAT <<- build_divmat(edges, DIVMAT_METRIC)
          message(sprintf("  division matrix: %d divisions (%s)",
                          nrow(.DIVMAT), DIVMAT_METRIC))
          save_both(render_divmat,
                    file.path(fig_dir, "Fig3_DivisionMatrix_Strengthened"), DM_W, DM_H)
        }
      }
    }
  }
  
  message("Done.")
}

# Runs when sourced (RStudio "Source") or via Rscript. Comment this out to load
# the functions without generating figures.
#
# Each mode in RENDER_MODES is rendered in turn into its own folder. Figure 3 is
# rendered first so that the computed tables it caches (.MEANR, .DIVR, .INVO)
# are still available for inspection after the run.
for (.mode in RENDER_MODES) {
  message("\n", strrep("=", 60))
  message("Rendering ", .mode)
  message(strrep("=", 60))
  apply_mode(.mode)
  main()
}