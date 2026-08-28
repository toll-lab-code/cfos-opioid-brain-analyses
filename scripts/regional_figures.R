# regional_figures.R
# =============================================================================
# Figure 2 (regional activation, sex-resolved) and Supplementary Figure S1
# (whole-brain overview). Reworked layout, Aug 2026:
#   - Old Figure 2 (lollipop / counts / whole-brain / exemplars) -> Supp. Fig S1.
#   - Old Figure 3 becomes the new main-text Figure 2, recombined with the
#     whole-brain total (now sex-stratified) and a sex-resolved exemplar grid.
#   - Old Figure 3 scatter and DMH/NOD exemplars are dropped.
# All text is set to the NPP >= 8 pt floor via the FS_* constants below.
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
#   Rscript scripts/regional_figures.R
# Input and output locations are set in the Config section below; see README.
# =============================================================================
#
# Run regional_analysis.R first: this script reads the statistics workbooks that
# it produces. Panels are rendered individually and assembled into the final
# figures separately.
#
# PANEL MANIFEST
# --- Figure 2 (new main-text figure; PIL-assembled) --------------------------
# Figure 2, Panel A: Diverging bar chart - males vs. females, Veh vs. Acute
#                    Morphine (formerly Fig 3A; unchanged content).
# Figure 2, Panel B: Sex-stratified whole-brain total cFos+ cells - Males | Females
#                    facets, all 4 conditions, per-sex one-way ANOVA (sex-split of
#                    the former Fig 2C whole-brain panel).
# Figure 2, Panel C: Sex-resolved exemplar grid - SNr / AVPV / PRE / PRM (columns)
#                    x Males / Females (rows); per-region free y shared across sex;
#                    within-sex FDR brackets (formerly Fig 3B recipe, regions swapped).
# Figure 2, Panel D: Representative light-sheet images (MBO) - composited into the
#                    figure separately, not produced by this script.
# --- Supplementary Figure S1 (former Figure 2; PIL-assembled) ----------------
# Fig S1, Panel A:   Lollipop - Veh vs Acute Morphine; FDR-significant regions,
#                    top 30 by Cohen's |d| (read live from the stats). Table S2.
# Fig S1, Panel B:   Significant-region count summary - uncorrected & FDR counts
#                    per drug condition (combined sex).
# Fig S1, Panel C:   Whole-brain total cFos+ cells, by condition (combined sex;
#                    one-way ANOVA).
# Fig S1, Panel D:   Faceted exemplar bars - MBO, ASO, TRS, MEA; all 4 conditions;
#                    combined sex.
#
# Output: 1200 dpi LZW TIFF + matching PNG per panel, for PIL assembly
# =============================================================================

# ---- Packages ----------------------------------------------------------------
required_pkgs <- c("tidyverse", "readxl", "showtext", "ggbeeswarm", "scales", "ggh4x")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages. Please run:\n  install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}
library(tidyverse)
library(readxl)
library(showtext)
library(ggbeeswarm)
library(scales)
library(ggh4x)      # facet_grid2(): per-region free y with regions across, sex down

# ---- Font sizes (NPP compliance: every text element >= 8 pt at final size) ----
# TARGET JOURNAL: Neuropsychopharmacology. FIGURE_SPEC.md v2.1 still describes a
# PNAS submission (178 mm double column, 7 pt ticks / 6 pt legend). NPP is the
# live target -- the 8 pt floor below is correct; the spec file is stale.
# Neuropsychopharmacology artwork spec floors all figure text at 8 pt on the
# printed page. Every panel is rendered at its final slot size and placed 1:1
# with no scaling, so the pt values below map EXACTLY to the printed page.
# Retarget another journal by editing these four values in ONE place; never
# scatter raw sizes through the panels and never rescale the assembled artwork.
FS_TITLE   <- 8    # axis titles + facet/strip titles (pt); strip titles stay bold
FS_TICK    <- 8    # axis tick labels, region labels (pt)
FS_LEGEND  <- 8    # legend title + text (pt); legend title stays bold
# geom_text()/annotate() take size in MILLIMETRES, not points: mm = pt / .pt
FS_MARK_MM <- 8 / ggplot2::.pt   # ~2.81 mm = 8 pt: significance stars & in-plot text

# ---- Font --------------------------------------------------------------------
# Figures use Arial. The default search path is the macOS system font directory;
# set the ARIAL_DIR environment variable to point somewhere else. If the font
# files are not found, the script falls back to the device default sans family,
# which alters glyph metrics slightly but nothing else about the figures.
arial_dir   <- Sys.getenv("ARIAL_DIR", "/System/Library/Fonts/Supplemental")
arial_files <- file.path(arial_dir, c("Arial.ttf", "Arial Bold.ttf",
                                      "Arial Italic.ttf", "Arial Bold Italic.ttf"))
FONT_FAMILY <- "sans"
if (all(file.exists(arial_files))) {
  font_add("Arial", regular = arial_files[1], bold = arial_files[2],
           italic = arial_files[3], bolditalic = arial_files[4])
  FONT_FAMILY <- "Arial"
} else {
  message("Arial not found in ", arial_dir, ". Falling back to the default sans ",
          "family; set ARIAL_DIR to the directory containing Arial.ttf to ",
          "reproduce the published figures exactly.")
}
showtext_auto()
showtext_opts(dpi = 1200)

# ---- Config ------------------------------------------------------------------
# Paths are relative to the repository root. Run from the repo root, e.g.:
#   Rscript scripts/regional_figures.R
# Every location is overridable with an environment variable, so a working copy
# laid out differently needs no edit to this file, e.g.
#   CFOS_DATA_DIR=~/Data/Drug CFOS_OUTPUT_DIR=~/Data/Drug CFOS_SCRIPT_DIR=~/Data/Drug \
#   CFOS_VOL_DIR=~/Data CFOS_RAW_XLSX="Drug Composite.xlsx" \
#     Rscript scripts/regional_figures.R
DATA_DIR   <- path.expand(Sys.getenv("CFOS_DATA_DIR",   "data"))     # Drug_Composite.xlsx, region_division_lookup.csv
OUTPUT_DIR <- path.expand(Sys.getenv("CFOS_OUTPUT_DIR", "results"))  # stats xlsx in; figures out under <OUTPUT_DIR>/figures/
SCRIPT_DIR <- path.expand(Sys.getenv("CFOS_SCRIPT_DIR", "scripts"))  # divisions.R
VOL_DIR    <- path.expand(Sys.getenv("CFOS_VOL_DIR",    DATA_DIR))   # ccfv3_volumes.xlsx
RAW_XLSX   <- Sys.getenv("CFOS_RAW_XLSX", "Drug_Composite.xlsx")
VOL_XLSX   <- Sys.getenv("CFOS_VOL_XLSX", "ccfv3_volumes.xlsx")

raw_file   <- file.path(DATA_DIR, RAW_XLSX)
vol_file   <- file.path(VOL_DIR,  VOL_XLSX)

# This first block builds the Supplementary Figure S1 panels (the former main-
# text Figure 2). output_dir is reassigned to figures/Figure2 before the new
# main-text Figure 2 section further down.
output_dir <- file.path(OUTPUT_DIR, "figures", "FigureS1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -- Colors --------------------------------------------------------------------
# The division palette and the region-to-division lookup are both defined in
# divisions.R, so the 13-division scheme is identical across every figure.
source(file.path(SCRIPT_DIR, "divisions.R"))
structure_colors <- division_colors            # 13 divisions, canonical palette

# Canonical region -> division lookup (Allen CCFv3 divisions, with the lab's two
# functional groupings: amygdalar nuclei kept as one Amygdala division; lateral
# septum grouped with Pallidum). Every panel resolves `structure` from THIS table
# by acronym, so the grouping cannot drift between figures. Generated from the
# Node_Roles Structure column (region_division_lookup.csv).
division_lookup_file <- file.path(DATA_DIR, "region_division_lookup.csv")
.divtab <- read.csv(division_lookup_file, stringsAsFactors = FALSE)
div_map <- setNames(.divtab$Structure, .divtab$Region)

# Short display labels (acronym -> abbreviated name) for the two label-heavy
# panels (lollipop, diverging bars). The DATA stays live from the workbooks;
# only the printed region name is shortened so the labels fit the panel width at
# 8 pt. Any acronym not listed falls back to the full region_name.
short_label <- c(
  AHN  = "Anterior hypothalamic nucleus",          ASO  = "Accessory supraoptic group",
  AT   = "Anterior tegmental nucleus",             ATN  = "Anterior group, dorsal thalamus",
  AVPV = "Anteroventral periventricular nucleus",  BAC  = "Bed nucleus of anterior commissure",
  CP   = "Caudoputamen",                           GPe  = "Globus pallidus, external",
  LHA  = "Lateral hypothalamic area",              LPO  = "Lateral preoptic area",
  MBO  = "Mammillary body",                        MEA  = "Medial amygdalar nucleus",
  MEPO = "Median preoptic nucleus",                MPN  = "Medial preoptic nucleus",
  MRN  = "Midbrain reticular nucleus",             NB   = "Nucleus of brachium of inf. colliculus",
  NDB  = "Diagonal band nucleus",                  NLL  = "Nucleus of the lateral lemniscus",
  NOD  = "Nodulus",                                PGRN = "Paragigantocellular reticular nucleus",
  PMd  = "Dorsal premammillary",                   POST = "Postsubiculum",
  PPN  = "Pedunculopontine nucleus",               PRE  = "Presubiculum",
  PRM  = "Paramedian lobule",                      PVH  = "Paraventricular hypothalamic nucleus",
  PVa  = "PVN, anterior part",                     PVi  = "PVN, intermediate part",
  PVpo = "PVN, preoptic part",                     PeF  = "Perifornical nucleus",
  RAmb = "Midbrain raphe nuclei",                  RR   = "Retrorubral area",
  RT   = "Reticular nucleus of thalamus",          SAG  = "Nucleus sagulum",
  SBPV = "Subparaventricular zone",                SCm  = "Superior colliculus, motor",
  SCs  = "Superior colliculus, sensory",           SI   = "Substantia innominata",
  SNr  = "Substantia nigra, reticular part",       TRS  = "Triangular nucleus of septum",
  TT   = "Taenia tecta",                           VTA  = "Ventral tegmental area",
  VeCB = "Vestibulocerebellar nucleus"
)

condition_colors <- c(
  "Vehicle"            = "#9e9e9e",   # grey (control)
  "Acute Morphine"     = "#c98ba0",   # dusty rose
  "Morphine-Dependent" = "#8f72b3",   # dusty mauve-purple
  "Ro 64-6198"         = "#4fae9f"    # soft teal-green
)

condition_order <- c("Vehicle", "Acute Morphine", "Morphine-Dependent", "Ro 64-6198")

sex_colors <- c("Males" = "#3f6fa3", "Females" = "#7d2452")   # muted blue / wine

# -- Statistics workbooks ------------------------------------------------------
# Defined up front so every panel (lollipop, counts, diverging, exemplars) can
# read its region set / effect sizes / significance live from the same files and
# cannot drift from the statistics. Produced by the analysis script.
stat_primary_file <- file.path(OUTPUT_DIR, "Drug_Statistical_Results_Primary.xlsx")
stat_bysex_file   <- file.path(OUTPUT_DIR, "Drug_Statistical_Results_BySex.xlsx")


# -- Save helper -------------------------------------------------------------
# Writes each panel twice at 1200 dpi: a TIFF (LZW) for archival/submission and
# a matching PNG. PNG is lossless and carries the dpi, so panels land at their
# true physical size when composited (no JPEG artifacts, no mis-scaling).
# Pass the .tiff path; the .png path is derived from it.
save_panel <- function(plot, path_tiff, width, height) {
  ggsave(path_tiff, plot, width = width, height = height, units = "mm",
         dpi = 1200, compression = "lzw", device = "tiff")
  ggsave(sub("\\.tiff?$", ".png", path_tiff), plot,
         width = width, height = height, units = "mm",
         dpi = 1200, device = "png")
}


# =============================================================================
# PANEL A - Lollipop (S1): FDR-significant regions (q<0.05, Vehicle vs. Acute
#           Morphine, combined sex), top 30 by Cohen's |d|. Points are filled
#           with the structure color. Region set and effect sizes are read at
#           run time from Drug_Statistical_Results_Primary.xlsx so the panel
#           cannot drift from the statistics. Full list in Table S2.
# =============================================================================

.figS1a_raw <- suppressMessages(read_excel(stat_primary_file, sheet = "Veh vs Mor"))

figS1a_data <- .figS1a_raw %>%
  transmute(
    acronym,
    label   = dplyr::coalesce(unname(short_label[acronym]), region_name),
    abs_d   = abs(cohens_d),
    fdr_sig = p_fdr < 0.05
  ) %>%
  filter(!is.na(abs_d), fdr_sig %in% TRUE) %>%   # FDR-significant (q<0.05)
  arrange(desc(abs_d)) %>%
  slice_head(n = 30) %>%                          # top 30 by |d|
  mutate(
    label      = factor(label, levels = rev(label)),
    structure  = factor(unname(div_map[acronym]), levels = names(structure_colors)),
    point_fill = unname(structure_colors[as.character(structure)])
  )

message(sprintf("Fig S1 Panel A: %d FDR-significant regions shown (top |d|)",
                nrow(figS1a_data)))

p_figS1a <- ggplot(figS1a_data, aes(x = abs_d, y = label, color = structure)) +
  
  geom_segment(aes(x = 0, xend = abs_d, yend = label),
               linewidth = 0.55, lineend = "round") +
  
  # All plotted regions are FDR-significant; point filled with the structure color.
  geom_point(aes(fill = point_fill),
             size = 2.4, stroke = 0.6, shape = 21,
             show.legend = c(colour = FALSE, fill = FALSE)) +

  scale_color_manual(values = structure_colors, name = "Structure",
                     guide = guide_legend(order = 1)) +
  scale_fill_identity() +

  # Limit is derived from the data, not hard-coded: with `expand = c(0,0)` a
  # hard limit silently clips any bar that exceeds it on a re-run.
  scale_x_continuous(
    name   = expression(paste("Cohen's |", italic(d), "|")),
    limits = c(0, ceiling(max(figS1a_data$abs_d) * 4) / 4 + 0.05),
    breaks = scales::pretty_breaks(n = 5),
    expand = c(0, 0)
  ) +
  
  theme_classic(base_family = FONT_FAMILY, base_size = 8) +
  theme(
    axis.title.x        = element_text(size = FS_TITLE, margin = margin(t = 4)),
    axis.title.y        = element_blank(),
    axis.text.x         = element_text(size = FS_TICK, color = "black"),
    axis.text.y         = element_text(size = FS_TICK, color = "black", hjust = 1),
    axis.line.x         = element_line(linewidth = 0.4, color = "#333333"),
    axis.line.y         = element_blank(),
    axis.ticks.y        = element_blank(),
    axis.ticks.x        = element_line(linewidth = 0.3, color = "#333333"),
    axis.ticks.length.x = unit(2, "pt"),
    panel.grid.major.x  = element_line(linewidth = 0.2, color = "#eeeeee"),
    panel.grid.major.y  = element_blank(),
    panel.grid.minor    = element_blank(),
    legend.title        = element_text(size = FS_LEGEND, face = "bold"),
    legend.text         = element_text(size = FS_LEGEND),
    legend.key.size     = unit(8, "pt"),
    legend.spacing.y    = unit(2, "pt"),
    legend.position     = "right",
    legend.box.spacing  = unit(2, "pt"),   # pull legend in close to the plot
    legend.background   = element_rect(fill = "white", color = NA),
    legend.margin       = margin(2, 2, 2, 2),
    plot.margin         = margin(8, 2, 6, 4),
    plot.background     = element_rect(fill = "white", color = NA),
    plot.caption        = element_text(size = FS_LEGEND, color = "#555555", hjust = 0)
  ) +
  
  labs(tag = NULL)

save_panel(p_figS1a, file.path(output_dir, "FigureS1_PanelA.tiff"),
           width = 114, height = 110)   # composite slot (Fig S1 left column)
cat("Fig S1 Panel A saved (lollipop).\n")


# =============================================================================
# PANEL D (S1) - Faceted exemplar bars: MBO, ASO, TRS, MEA; all 4 conditions;
#               combined sex. (Relettered from E: S1 has no representative-image
#               strip - that panel moved to main-text Figure 2.)
# =============================================================================

# -- Load CCFv3 volumes --------------------------------------------------------
vol_raw <- read_excel(vol_file, skip = 1)
vol_map <- setNames(
  as.numeric(vol_raw[["Mean Volume (m)"]]),
  trimws(as.character(vol_raw[["abbreviation"]]))
)

# -- Load and normalize per-animal data ---------------------------------------
load_and_normalize <- function(region_acronym, sheet_name, sex_label) {
  raw <- read_excel(raw_file, sheet = sheet_name)
  
  acr_col    <- grep("acronym", names(raw), ignore.case = TRUE)[1]
  region_row <- which(trimws(as.character(raw[[acr_col]])) == region_acronym)
  if (length(region_row) == 0)
    stop(paste("Region not found:", region_acronym, "in", sheet_name))
  
  region_vol <- vol_map[region_acronym]
  if (is.na(region_vol) || region_vol == 0)
    stop(paste("Volume not found:", region_acronym))
  
  all_cols  <- names(raw)
  data_cols <- grep("^(Vehicle|Morphine|Chronic|Ro)", all_cols, value = TRUE)
  
  cond_map <- list(
    "Vehicle"            = grep("^Vehicle",  data_cols, value = TRUE),
    "Acute Morphine"     = grep("^Morphine", data_cols, value = TRUE),
    "Morphine-Dependent" = grep("^Chronic",  data_cols, value = TRUE),
    "Ro 64-6198"         = grep("^Ro",       data_cols, value = TRUE)
  )
  
  rows_list <- list()
  for (cond in names(cond_map)) {
    cols <- cond_map[[cond]]
    if (length(cols) == 0) next
    vals <- suppressWarnings(as.numeric(raw[region_row, cols]))
    vals <- vals[!is.na(vals) & vals > 0]
    if (length(vals) > 0) {
      rows_list[[cond]] <- tibble(
        condition = cond,
        sex       = sex_label,
        value     = log10(vals / region_vol)
      )
    }
  }
  bind_rows(rows_list)
}

# -- Significance lookup from statistical result files ------------------------
# Reads sig_fdr (Benjamini-Hochberg corrected) directly from the ANOVA output so
# brackets always reflect the current statistics and are consistent with the FDR
# encoding used in the lollipop and diverging chart. Rerunning the analysis flows
# straight through here. (stat_primary_file / stat_bysex_file are defined up top.)

# Map comparison name -> x-axis positions (Veh=1, Mor=2, MorDep=3, Ro=4)
comparison_positions <- list(
  "Veh vs Mor"    = c(1, 2),
  "Veh vs MorDep" = c(1, 3),
  "Veh vs Ro"     = c(1, 4),
  "Mor vs MorDep" = c(2, 3),
  "Mor vs Ro"     = c(2, 4)
)

# Look up sig_fdr for one region in one sheet; returns "ns" if absent
get_sig <- function(stat_file, sheet, region_acronym) {
  df <- suppressMessages(read_excel(stat_file, sheet = sheet))
  acr_col <- grep("acronym", names(df), ignore.case = TRUE)[1]
  row <- df[trimws(as.character(df[[acr_col]])) == region_acronym, ]
  if (nrow(row) == 0) return("ns")
  sig <- as.character(row[["sig_fdr"]][1])
  if (is.na(sig)) return("ns")
  sig
}

# Build bracket list (combined sex - Primary file). Comparisons listed narrowest
# first so they stack without overlapping. Only FDR q<0.05 and below (*, **, ***)
# produce a bracket; FDR trends q<0.10 (dagger) and ns are not drawn.
build_brackets_combined <- function(region_acronym, comparisons) {
  brackets <- list()
  for (comp in comparisons) {
    sig <- get_sig(stat_primary_file, comp, region_acronym)
    if (sig %in% c("*", "**", "***")) {
      pos <- comparison_positions[[comp]]
      brackets[[length(brackets) + 1]] <- list(x1 = pos[1], x2 = pos[2], sig = sig)
    }
  }
  brackets
}

# Build bracket list within a sex (BySex file). Same p<0.05 threshold.
build_brackets_sex <- function(region_acronym, sex_label, comparisons) {
  sex_suffix <- if (sex_label == "Males") "Males" else "Fem"
  brackets <- list()
  for (comp in comparisons) {
    sheet <- paste0(comp, " - ", sex_suffix)
    sig <- get_sig(stat_bysex_file, sheet, region_acronym)
    if (sig %in% c("*", "**", "***")) {
      pos <- comparison_positions[[comp]]
      brackets[[length(brackets) + 1]] <- list(x1 = pos[1], x2 = pos[2], sig = sig)
    }
  }
  brackets
}

# Region full names for titles
region_labels <- c(
  MBO = "Mammillary body",
  ASO = "Accessory supraoptic group",
  TRS = "Triangular nucleus of septum",
  MEA = "Medial amygdalar nucleus"
)

# Significance annotations - read from ANOVA output at runtime (vs Vehicle)
# Combined sex: Veh vs Mor (1-2), Veh vs MorDep (1-3), Veh vs Ro (1-4)
figS1d_comparisons <- c("Veh vs Mor", "Veh vs MorDep", "Veh vs Ro")

# -- Faceted exemplar bars (one panel: MBO, ASO, TRS, MEA) ---------------------
# Single faceted panel with a shared x-axis and per-region (free) y-scales,
# replacing the four standalone plots. This removes the repeated y-axis and
# x-axis ink and reads as one unit. Significance brackets are added per facet
# from the same FDR lookup used everywhere else.
cat("Loading exemplar-bar data...\n")
bar_regions <- c("MBO","ASO","TRS","MEA")

bar_df <- purrr::map_dfr(bar_regions, function(acr) {
  m <- load_and_normalize(acr, "Clustered Males",   "Males")
  f <- load_and_normalize(acr, "Clustered Females", "Females")
  # NB: no `filter(value > 0)` here. `value` is log10(density), so that filter
  # silently dropped any animal with density < 1 cell/mm3 from the FIGURE while
  # the statistics kept it. Zero counts are already removed inside
  # load_and_normalize() (`vals > 0`, on raw counts), which is the correct place.
  bind_rows(m, f) %>% dplyr::select(-sex) %>%
    dplyr::mutate(region = acr)
}) %>%
  dplyr::mutate(condition = factor(condition, levels = condition_order),
                region    = factor(region,    levels = bar_regions))

# Per-region y-limits and bracket positions (facets use free y-scales)
brk_list <- list(); lim_list <- list()
for (acr in bar_regions) {
  d    <- dplyr::filter(bar_df, region == acr)
  ymax <- max(d$value, na.rm = TRUE); ymin <- min(d$value, na.rm = TRUE)
  brks <- Filter(function(b) b$sig %in% c("*","**","***"),
                 build_brackets_combined(acr, figS1d_comparisons))
  step <- max((ymax - ymin) * 0.14, 0.13)
  if (length(brks)) {
    for (i in seq_along(brks)) {
      b <- brks[[i]]
      brk_list[[length(brk_list) + 1]] <- data.frame(
        region = acr, x1 = b$x1, x2 = b$x2, xm = (b$x1 + b$x2) / 2,
        y = ymax + 0.10 + step * (i - 1), sig = b$sig,
        stringsAsFactors = FALSE)
    }
    ytop <- ymax + 0.10 + step * (length(brks) - 1) + step * 0.9
  } else {
    ytop <- ymax + 0.15
  }
  lim_list[[length(lim_list) + 1]] <- data.frame(
    region = acr, ylo = ymin - (ymax - ymin) * 0.05, yhi = ytop)
}
brk_df <- if (length(brk_list)) do.call(rbind, brk_list) else
  data.frame(region = character(), x1 = double(), x2 = double(),
             xm = double(), y = double(), sig = character())
lim_df <- do.call(rbind, lim_list)
brk_df$region <- factor(brk_df$region, levels = bar_regions)
lim_df$region <- factor(lim_df$region, levels = bar_regions)

# Invisible points pin each facet's y-range; tick marks for the brackets
lim_long <- tidyr::pivot_longer(lim_df, c(ylo, yhi), values_to = "value") %>%
  dplyr::mutate(condition = factor("Vehicle", levels = condition_order))
brk_tick <- if (nrow(brk_df)) {
  tidyr::pivot_longer(brk_df, c(x1, x2), values_to = "x") %>%
    dplyr::mutate(yend = y - 0.05)
} else brk_df

p_bars <- ggplot(bar_df, aes(x = condition, y = value)) +
  
  geom_beeswarm(aes(color = condition),
                size = 0.8, shape = 16, cex = 2.6, alpha = 0.80) +
  
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.25, linewidth = 0.45, color = "black") +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.40, color = "black", fill = NA, linewidth = 0.45) +
  
  geom_blank(data = lim_long, aes(x = condition, y = value)) +
  
  geom_segment(data = brk_df, aes(x = x1, xend = x2, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.3) +
  geom_segment(data = brk_tick, aes(x = x, xend = x, y = y, yend = yend),
               inherit.aes = FALSE, linewidth = 0.3) +
  geom_text(data = brk_df, aes(x = xm, y = y + 0.04, label = sig),
            inherit.aes = FALSE, family = FONT_FAMILY, size = FS_MARK_MM, vjust = 0) +
  
  scale_color_manual(values = condition_colors, guide = "none") +
  scale_x_discrete(labels = c("Vehicle"="Veh","Acute Morphine"="Mor",
                              "Morphine-Dependent"="MorDep","Ro 64-6198"="Ro")) +
  scale_y_continuous(name = expression(log[10](cFos^"+" ~ cells/mm^3)),
                     breaks = pretty_breaks(n = 4),
                     labels = function(x) sprintf("%.1f", x),   # uniform 1-decimal ticks across facets
                     expand = expansion(mult = c(0.02, 0.02))) +
  
  facet_wrap(~ region, nrow = 1, scales = "free_y") +
  coord_cartesian(clip = "off") +
  
  theme_classic(base_family = FONT_FAMILY, base_size = 8) +
  theme(
    strip.text        = element_text(size = FS_TITLE, face = "bold", margin = margin(b = 2)),
    strip.background  = element_blank(),
    axis.title.x      = element_blank(),
    axis.title.y      = element_text(size = FS_TITLE, margin = margin(r = 2)),
    axis.text.x       = element_text(size = FS_TICK, color = "black",
                                     angle = 30, hjust = 1, vjust = 1),
    axis.text.y       = element_text(size = FS_TICK, color = "black"),
    axis.line         = element_line(linewidth = 0.4, color = "#333333"),
    axis.ticks        = element_line(linewidth = 0.3, color = "#333333"),
    axis.ticks.length = unit(2, "pt"),
    panel.spacing     = unit(5, "pt"),
    legend.position   = "none",
    plot.margin       = margin(4, 6, 2, 2),
    plot.background   = element_rect(fill = "white", color = NA)
  )

save_panel(p_bars, file.path(output_dir, "FigureS1_PanelD_bars.tiff"),
           width = 169, height = 46)   # composite slot (Fig S1 bottom row, full width)
cat("Fig S1 Panel D (faceted exemplar bars) saved.\n")


# =============================================================================
# PANEL B - Significant-region count summary (combined sex)
# Number of regions significant vs. Vehicle for each drug condition, at two
# thresholds: uncorrected p<0.05 (faded) and FDR q<0.05 (solid, a subset of the
# uncorrected count). Counts are computed at run time from
# Drug_Statistical_Results_Primary.xlsx so the panel cannot drift from the
# statistics. The faded->solid collapse at Ro shows the restricted NOP-agonist
# profile (few regions, none surviving FDR). Full counts in Table S2.
# =============================================================================

count_levels <- c("Acute Morphine", "Morphine-Dependent", "Ro 64-6198")

# For each Vehicle-vs-condition sheet: regions at uncorrected p_tukey < 0.05 and
# at FDR q < 0.05 (p_fdr column).
.count_sheets <- c("Acute Morphine"     = "Veh vs Mor",
                   "Morphine-Dependent" = "Veh vs MorDep",
                   "Ro 64-6198"         = "Veh vs Ro")
.sig_counts <- function(sheet) {
  d <- suppressMessages(read_excel(stat_primary_file, sheet = sheet))
  c(unc = sum(d$p_tukey < 0.05, na.rm = TRUE),
    fdr = sum(d$p_fdr   < 0.05, na.rm = TRUE))
}
.counts <- vapply(.count_sheets, .sig_counts, numeric(2))   # rows: unc, fdr

count_data <- tibble(
  condition = factor(rep(count_levels, 2), levels = count_levels),
  tier      = factor(rep(c("Tukey p<0.05", "FDR q<0.05"), each = 3),
                     levels = c("Tukey p<0.05", "FDR q<0.05")),
  count     = c(.counts["unc", ], .counts["fdr", ])
) %>% arrange(tier)   # uncorrected drawn first (behind); FDR overlaid on top
message(sprintf("Fig S1 Panel B counts (unc/fdr): Mor %d/%d, MorDep %d/%d, Ro %d/%d",
                .counts["unc","Acute Morphine"], .counts["fdr","Acute Morphine"],
                .counts["unc","Morphine-Dependent"], .counts["fdr","Morphine-Dependent"],
                .counts["unc","Ro 64-6198"], .counts["fdr","Ro 64-6198"]))

p_counts <- ggplot(count_data, aes(condition, count, fill = condition, alpha = tier)) +
  
  geom_col(position = "identity", width = 0.68) +
  
  # Count labels: uncorrected total above each bar, FDR count atop its solid part
  geom_text(data = subset(count_data, tier == "Tukey p<0.05"),
            aes(label = count), vjust = -0.4, size = FS_MARK_MM, fontface = "bold",
            color = "#555555", family = FONT_FAMILY, show.legend = FALSE) +
  geom_text(data = subset(count_data, tier == "FDR q<0.05" & count > 0),
            aes(label = count), vjust = -0.4, size = FS_MARK_MM, fontface = "bold",
            color = "black", family = FONT_FAMILY, show.legend = FALSE) +
  
  scale_fill_manual(values = condition_colors, guide = "none") +
  scale_alpha_manual(values = c("Tukey p<0.05" = 0.32, "FDR q<0.05" = 1.0),
                     name = NULL) +
  
  scale_y_continuous(name = "Significant regions\n(vs. Vehicle)",   # 2 lines so it fits the panel at 8 pt
                     # top adapts to the live counts (was hardcoded 82) so tall
                     # bars and their labels are never clipped
                     limits = c(0, max(90, ceiling(max(count_data$count) * 1.08 / 10) * 10)),
                     breaks = seq(0, 80, 20),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_x_discrete(labels = c("Acute Morphine"     = "Mor",
                              "Morphine-Dependent" = "MorDep",
                              "Ro 64-6198"         = "Ro")) +
  
  guides(alpha = guide_legend(override.aes = list(fill = "#777777"))) +
  
  theme_classic(base_family = FONT_FAMILY, base_size = 8) +
  theme(
    axis.title.x         = element_blank(),
    axis.title.y         = element_text(size = FS_TITLE, margin = margin(r = 3)),
    axis.text.x          = element_text(size = FS_TICK, color = "black",
                                        angle = 30, hjust = 1, vjust = 1),
    axis.text.y          = element_text(size = FS_TICK, color = "black"),
    axis.line            = element_line(linewidth = 0.4, color = "#333333"),
    axis.ticks           = element_line(linewidth = 0.3, color = "#333333"),
    axis.ticks.length    = unit(2, "pt"),
    legend.position      = "top",              # out of the plot -> no overlap with bars
    legend.direction     = "horizontal",    # single row, right-aligned
    legend.justification = "right",
    legend.text          = element_text(size = FS_LEGEND),
    legend.key.size      = unit(7, "pt"),
    legend.key.spacing.y = unit(1, "pt"),
    legend.margin        = margin(0, 0, 1, 0),
    legend.background    = element_rect(fill = "white", color = NA),
    plot.margin          = margin(4, 4, 4, 4),
    plot.background      = element_rect(fill = "white", color = NA)
  )

save_panel(p_counts, file.path(output_dir, "FigureS1_PanelB_CountSummary.tiff"),
           width = 50, height = 47)   # composite slot (Fig S1 right sidebar, top)
cat("Fig S1 Panel B (count summary) saved.\n")


# =============================================================================
# PANEL C - Whole-brain total cFos+ cells per animal, by condition (combined sex)
# Per-animal total = sum of per-region counts across all 198 analysis regions,
# pooled across sex. Supports the non-significant whole-brain total: a one-way
# ANOVA omnibus p is computed at runtime and annotated as a bare p-value (no
# "n.s." marker). Raw counts (zeros retained - they are real here, unlike in the
# log-density panels, where log10(0) is dropped).
# NOTE: assumes the 198 rows are the non-overlapping analysis region set (no
# parent/child double-counting); the absolute total scales with that set but the
# across-condition comparison is unaffected.
# =============================================================================

read_totals <- function(sheet, sex_label) {
  raw       <- read_excel(raw_file, sheet = sheet)
  data_cols <- grep("^(Vehicle|Morphine|Chronic|Ro)", names(raw), value = TRUE)
  cond_map  <- list(
    "Vehicle"            = grep("^Vehicle",  data_cols, value = TRUE),
    "Acute Morphine"     = grep("^Morphine", data_cols, value = TRUE),
    "Morphine-Dependent" = grep("^Chronic",  data_cols, value = TRUE),
    "Ro 64-6198"         = grep("^Ro",       data_cols, value = TRUE)
  )
  rows_list <- list()
  for (cond in names(cond_map)) {
    cols <- cond_map[[cond]]
    if (length(cols) == 0) next
    totals <- vapply(cols, function(cc)
      sum(suppressWarnings(as.numeric(raw[[cc]])), na.rm = TRUE), numeric(1))
    rows_list[[cond]] <- tibble(condition = cond, sex = sex_label,
                                animal = cols, total = as.numeric(totals))
  }
  bind_rows(rows_list)
}

total_data <- bind_rows(
  read_totals("Clustered Males",   "Males"),
  read_totals("Clustered Females", "Females")
) %>%
  mutate(condition = factor(condition, levels = condition_order))

# Omnibus one-way ANOVA on whole-brain totals (combined sex)
.aov_e <- aov(total ~ condition, data = total_data)
.p_e   <- summary(.aov_e)[[1]][["Pr(>F)"]][1]
.p_lab <- sprintf("ANOVA p = %.2f", .p_e)
cat(sprintf("Panel C - whole-brain total ANOVA: p = %.4f\n", .p_e))

p_e <- ggplot(total_data, aes(x = condition, y = total)) +
  
  geom_beeswarm(aes(color = condition),
                size = 0.8, shape = 16, cex = 2.6, alpha = 0.80) +
  
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.25, linewidth = 0.5, color = "black") +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.40, color = "black", fill = NA, linewidth = 0.5) +
  
  scale_color_manual(values = condition_colors, guide = "none") +
  
  scale_y_continuous(
    name   = expression("Total cFos"^"+" ~ "cells (whole brain)"),
    labels = scales::comma,
    expand = expansion(mult = c(0.03, 0.12))
  ) +
  
  scale_x_discrete(
    labels = c("Vehicle"            = "Veh",
               "Acute Morphine"     = "Mor",
               "Morphine-Dependent" = "MorDep",
               "Ro 64-6198"         = "Ro")
  ) +
  
  annotate("text", x = 2.5, y = Inf, vjust = 1.6, label = .p_lab,
           size = FS_MARK_MM, family = FONT_FAMILY, color = "#555555") +
  
  coord_cartesian(clip = "off") +
  
  theme_classic(base_family = FONT_FAMILY, base_size = 8) +
  theme(
    axis.title.x      = element_blank(),
    axis.title.y      = element_text(size = FS_TITLE, margin = margin(r = 3)),
    axis.text.x       = element_text(size = FS_TICK, color = "black",
                                     angle = 30, hjust = 1, vjust = 1),
    axis.text.y       = element_text(size = FS_TICK, color = "black"),
    axis.line         = element_line(linewidth = 0.4, color = "#333333"),
    axis.ticks        = element_line(linewidth = 0.3, color = "#333333"),
    axis.ticks.length = unit(2, "pt"),
    legend.position   = "none",
    plot.margin       = margin(6, 8, 4, 4),
    plot.background   = element_rect(fill = "white", color = NA)
  )

save_panel(p_e, file.path(output_dir, "FigureS1_PanelC_TotalCells.tiff"),
           width = 50, height = 58)   # composite slot (Fig S1 right sidebar, bottom)
cat("Fig S1 Panel C (whole-brain total cFos+ cells) saved.\n")


cat("\nAll Supplementary Figure S1 panels saved to:", output_dir, "\n")


# =============================================================================
# FIGURE 2 (new main text) - redirect output to Figure2 subdirectory
# Panels: A diverging bars | B sex-split whole-brain totals |
#         C sex-resolved exemplar grid | D representative images (composited separately).
# total_data (whole-brain per-animal totals, with sex) was built in the S1
# section above and is reused here for Panel B.
# =============================================================================
output_dir <- file.path(OUTPUT_DIR, "figures", "Figure2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# PANEL A - Diverging bar chart: males vs. females, Veh vs. Acute Morphine
# Males extend right (+), females extend left (-)
#
# INCLUSION RULE (principled, stated): a region is plotted if it survives FDR
# (q<0.05) in EITHER sex, OR is significant at uncorrected p<0.05 in FEMALES.
# Under acute morphine no female region survives FDR, so this resolves to the
# 20 male FDR-significant regions + the 6 additional female uncorrected-p<0.05
# regions (the female footprint) = 26 regions. Values are read from
# Drug_Statistical_Results_BySex.xlsx ("Veh vs Mor - Males" / "- Fem").
#
# Opacity: full = FDR q<0.05, mid = uncorrected p<0.05 only, faint = ns.
# Tip dot = FDR-significant. Male-uncorrected-only regions (PAG, PPN, RN, DG,
# CA, GPi, etc.) are intentionally excluded by the rule; they remain in Table S2.
# =============================================================================

# Region set, effect sizes, and significance flags are read at run time from
# Drug_Statistical_Results_BySex.xlsx ("Veh vs Mor - Males" / "- Fem") so the
# panel cannot drift from the statistics.
.fig2a_m <- suppressMessages(read_excel(stat_bysex_file, sheet = "Veh vs Mor - Males")) %>%
  transmute(acronym,
            label = dplyr::coalesce(unname(short_label[acronym]), region_name),
            d_m = abs(cohens_d), sig_m = p_tukey < 0.05, fdr_m = p_fdr < 0.05)
.fig2a_f <- suppressMessages(read_excel(stat_bysex_file, sheet = "Veh vs Mor - Fem")) %>%
  transmute(acronym,
            d_f = abs(cohens_d), sig_f = p_tukey < 0.05, fdr_f = p_fdr < 0.05)

panel_a_data <- .fig2a_m %>%
  left_join(.fig2a_f, by = "acronym") %>%
  mutate(across(c(sig_m, sig_f, fdr_m, fdr_f), ~ tidyr::replace_na(., FALSE)),
         d_f = tidyr::replace_na(d_f, 0)) %>%
  # inclusion rule: FDR in either sex, or uncorrected p<0.05 in females
  filter(fdr_m | fdr_f | sig_f) %>%
  # Sort by male |d| descending - females mirror on left
  arrange(desc(d_m)) %>%
  mutate(
    label     = factor(label, levels = rev(label)),
    structure = factor(unname(div_map[acronym]), levels = names(structure_colors)),
    # Opacity tiers: full = FDR q<0.05, mid = uncorrected p<0.05 only, faint = ns.
    alpha_m   = dplyr::case_when(fdr_m ~ 1.0, sig_m ~ 0.55, TRUE ~ 0.22),
    alpha_f   = dplyr::case_when(fdr_f ~ 1.0, sig_f ~ 0.55, TRUE ~ 0.22)
  )

message(sprintf("Fig 2 Panel A: %d regions (%d male-FDR, %d female uncorrected-only)",
                nrow(panel_a_data),
                sum(panel_a_data$fdr_m, na.rm = TRUE),
                sum(panel_a_data$sig_f & !panel_a_data$fdr_m, na.rm = TRUE)))

# ---- Panel A axis limits and label offsets ----------------------------------
# Limits are derived from the data (a hard-coded limit with expand = c(0,0)
# silently clips any bar that exceeds it on a re-run).
.p2a_lo <- -(ceiling(max(panel_a_data$d_f, na.rm = TRUE) * 2) / 2 + 0.35)
.p2a_hi <-   ceiling(max(panel_a_data$d_m, na.rm = TRUE) * 2) / 2 + 0.35
# Offset of the "Females" / "Males" labels from the zero line, as a FRACTION of
# the axis range rather than an absolute number of data units. Panel width is
# fixed, so mm-per-data-unit scales inversely with the range: a fixed offset in
# data units changes its printed size whenever the limits move, which is how the
# labels ended up 0.02 mm from the axis line (touching it) at the previous
# limits. A fraction of the range keeps the printed gap constant. 0.068 gives
# ~1.2 mm of clear space on the tighter (Males) side at 8 pt.
.p2a_lab_off <- 0.068 * (.p2a_hi - .p2a_lo)

# Build diverging chart - males positive, females negative
# Reshape to long format for easier plotting
df_long_a <- bind_rows(
  panel_a_data %>%
    transmute(label, acronym, structure,
              value = d_m, side = "Males",
              alpha = alpha_m, sig = sig_m, fdr = fdr_m),
  panel_a_data %>%
    transmute(label, acronym, structure,
              value = -d_f, side = "Females",
              alpha = alpha_f, sig = sig_f, fdr = fdr_f)
) %>%
  # `sig` is p_tukey, which is ALREADY adjusted across the 6 within-region
  # pairwise comparisons -- it is not an uncorrected p. Labelled accordingly so
  # the figure, the table columns and the Results text all agree.
  dplyr::mutate(sig_tier = factor(dplyr::case_when(
    fdr ~ "FDR q<0.05", sig ~ "Tukey p<0.05", TRUE ~ "ns"),
    levels = c("FDR q<0.05", "Tukey p<0.05", "ns")))

p_a <- ggplot(df_long_a,
              aes(x = value, y = label, fill = structure, alpha = sig_tier)) +
  
  # Bars
  geom_col(width = 0.75, color = NA) +
  
  # Zero line
  geom_vline(xintercept = 0, linewidth = 0.5, color = "#333333") +
  
  # Solid dot at bar tip marks FDR-significant regions (q<0.05)
  geom_point(data = df_long_a %>% filter(fdr %in% TRUE),
             aes(x = value + sign(value) * 0.04, color = structure),
             size = 1.1, show.legend = FALSE) +
  
  scale_fill_manual(values  = structure_colors, name = "Structure",
                    guide = guide_legend(order = 1)) +
  scale_color_manual(values = structure_colors, name = "Structure", guide = "none") +
  scale_alpha_manual(
    name   = "vs. Vehicle",
    values = c("FDR q<0.05" = 1.0, "Tukey p<0.05" = 0.55, "ns" = 0.22),
    labels = c("FDR q<0.05  \u25CF", "Tukey p<0.05", "ns"),
    guide  = guide_legend(order = 2, override.aes = list(fill = "grey30"))) +
  
  # Symmetric-per-side limits derived from the data (see note on Fig S1 Panel A).
  # Females plot as -|d|, so the axis is mirrored: labels are absolute values.
  scale_x_continuous(
    name   = expression(paste("Cohen's |", italic(d), "|")),
    limits = c(.p2a_lo, .p2a_hi),
    breaks = scales::breaks_width(1),
    labels = function(x) as.character(abs(x)),
    expand = c(0, 0)
  ) +
  
  # Headroom above the top bar so the sex labels are not clipped
  scale_y_discrete(expand = expansion(add = c(0.6, 2.0))) +
  
  # Female / Male axis labels (placed in the headroom)
  # Anchor each label to its own side of the zero line with a small gap so they
  # never touch: Females right-aligned just left of 0, Males left-aligned just right.
  annotate("text", x = -.p2a_lab_off, y = length(levels(df_long_a$label)) + 1.3,
           label = "\u2190 Females", hjust = 1,
           size = FS_MARK_MM, family = FONT_FAMILY, fontface = "bold",
           color = sex_colors[["Females"]]) +
  annotate("text", x = .p2a_lab_off, y = length(levels(df_long_a$label)) + 1.3,
           label = "Males \u2192", hjust = 0,
           size = FS_MARK_MM, family = FONT_FAMILY, fontface = "bold",
           color = sex_colors[["Males"]]) +
  
  coord_cartesian(clip = "off") +
  
  theme_classic(base_family = FONT_FAMILY, base_size = 8) +
  theme(
    axis.title.x   = element_text(size = FS_TITLE, margin = margin(t = 4)),
    axis.title.y   = element_blank(),
    axis.text.x    = element_text(size = FS_TICK, color = "black"),
    axis.text.y    = element_text(size = FS_TICK, color = "black", hjust = 1),
    axis.line.x    = element_line(linewidth = 0.4, color = "#333333"),
    axis.line.y    = element_blank(),
    axis.ticks.y   = element_blank(),
    axis.ticks.x   = element_line(linewidth = 0.3, color = "#333333"),
    axis.ticks.length.x = unit(2, "pt"),
    panel.grid.major.x  = element_line(linewidth = 0.2, color = "#eeeeee"),
    panel.grid.major.y  = element_blank(),
    legend.title        = element_text(size = FS_LEGEND, face = "bold"),
    legend.text         = element_text(size = FS_LEGEND),
    legend.key.size     = unit(8, "pt"),
    legend.position     = "right",
    legend.background   = element_rect(fill = "white", color = NA),
    plot.margin         = margin(18, 4, 6, 4),
    plot.background     = element_rect(fill = "white", color = NA),
    plot.caption        = element_text(size = FS_LEGEND, color = "#555555", hjust = 0)
  ) +
  
  labs(tag = NULL)

save_panel(p_a, file.path(output_dir, "Figure2_PanelA.tiff"),
           width = 106, height = 110)   # composite slot (new Fig 2, top-left, diverging bars)
cat("Fig 2 Panel A saved (diverging bars).\n")


# =============================================================================
# PANEL B - Sex-stratified whole-brain total cFos+ cells (Males | Females)
# Sex-split of the S1 whole-brain panel: same per-animal totals, faceted by sex,
# with a per-sex one-way ANOVA omnibus p annotated in each facet (bare p, no
# "n.s." marker). Shared y across the two sex facets so Males and Females are
# directly comparable. Reuses total_data built in the S1 section above.
# =============================================================================
total_data_sex <- total_data %>%
  dplyr::mutate(sex = factor(sex, levels = c("Males", "Females")))

# Per-sex omnibus one-way ANOVA on whole-brain totals. Two-line label so it fits
# centered under each sex strip without colliding with the y-axis.
sex_totals_p <- total_data_sex %>%
  dplyr::group_by(sex) %>%
  dplyr::summarise(p = summary(aov(total ~ condition))[[1]][["Pr(>F)"]][1],
                   .groups = "drop") %>%
  dplyr::mutate(lab = sprintf("ANOVA\np = %.2f", p))
cat(sprintf("Fig 2 Panel B - per-sex whole-brain ANOVA: %s\n",
            paste(sprintf("%s p=%.4f", sex_totals_p$sex, sex_totals_p$p),
                  collapse = "; ")))

# Significance brackets vs Vehicle. Whole-brain totals are summed in-script and
# are NOT in the per-region stat workbooks, so these are computed at runtime:
# gated on the per-sex omnibus ANOVA (a sex only gets brackets if its ANOVA is
# p < 0.05), then pairwise Welch t-tests of each drug vs Vehicle, BH-adjusted
# across the three comparisons. Stars: * q<0.05, ** q<0.01, *** q<0.001.
star_of <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}
whole_vs_veh <- list(c("Acute Morphine", 2), c("Morphine-Dependent", 3), c("Ro 64-6198", 4))
y_gmax <- max(total_data_sex$total, na.rm = TRUE)
y_gmin <- min(total_data_sex$total, na.rm = TRUE)
b_step <- (y_gmax - y_gmin) * 0.09

totB_list <- list()
for (sx in levels(total_data_sex$sex)) {
  om_p <- sex_totals_p$p[sex_totals_p$sex == sx]
  if (length(om_p) == 0 || is.na(om_p) || om_p >= 0.05) next   # gate on omnibus
  d_sx <- dplyr::filter(total_data_sex, sex == sx)
  veh  <- d_sx$total[d_sx$condition == "Vehicle"]
  y_sx_max <- max(d_sx$total, na.rm = TRUE)
  praw <- vapply(whole_vs_veh, function(cc) {
    oth <- d_sx$total[d_sx$condition == cc[[1]]]
    if (length(veh) < 2 || length(oth) < 2) return(NA_real_)
    stats::t.test(veh, oth)$p.value
  }, numeric(1))
  padj <- p.adjust(praw, method = "BH")
  k <- 0
  for (i in seq_along(whole_vs_veh)) {
    if (star_of(padj[i]) == "ns") next
    k  <- k + 1
    x2 <- as.numeric(whole_vs_veh[[i]][[2]])
    totB_list[[length(totB_list) + 1]] <- data.frame(
      sex = sx, x1 = 1, x2 = x2, xm = (1 + x2) / 2,
      y = y_sx_max + b_step * (0.6 + (k - 1)), sig = star_of(padj[i]),
      stringsAsFactors = FALSE)
  }
}
totB_df <- if (length(totB_list)) do.call(rbind, totB_list) else
  data.frame(sex = character(), x1 = double(), x2 = double(),
             xm = double(), y = double(), sig = character())
if (nrow(totB_df)) totB_df$sex <- factor(totB_df$sex, levels = levels(total_data_sex$sex))
totB_tick <- if (nrow(totB_df))
  tidyr::pivot_longer(totB_df, c(x1, x2), values_to = "x") %>%
    dplyr::mutate(yend = y - b_step * 0.28) else totB_df

# Pin a shared y-range that clears the global max, the tallest bracket, and the
# two-line ANOVA label above it (identical in both facets so the axes match).
y_ceiling <- max(y_gmax, if (nrow(totB_df)) max(totB_df$y) else -Inf)
y_top_pin <- y_ceiling + b_step * 1.4   # just enough for the two-line ANOVA label; keeps the axis short
pinB <- tidyr::crossing(
  sex       = factor(levels(total_data_sex$sex), levels = levels(total_data_sex$sex)),
  total     = c(y_gmin, y_top_pin),
  condition = factor("Vehicle", levels = condition_order))

p_fig2b <- ggplot(total_data_sex, aes(x = condition, y = total)) +

  geom_beeswarm(aes(color = condition),
                size = 0.8, shape = 16, cex = 2.6, alpha = 0.80) +

  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.25, linewidth = 0.5, color = "black") +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.40, color = "black", fill = NA, linewidth = 0.5) +

  geom_blank(data = pinB, aes(x = condition, y = total)) +

  # Significance brackets (only a sex whose omnibus ANOVA is significant)
  geom_segment(data = totB_df, aes(x = x1, xend = x2, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.3) +
  geom_segment(data = totB_tick, aes(x = x, xend = x, y = y, yend = yend),
               inherit.aes = FALSE, linewidth = 0.3) +
  geom_text(data = totB_df, aes(x = xm, y = y + b_step * 0.15, label = sig),
            inherit.aes = FALSE, family = FONT_FAMILY, size = FS_MARK_MM, vjust = 0) +

  scale_color_manual(values = condition_colors, guide = "none") +

  scale_y_continuous(
    name   = expression("Total cFos"^"+" ~ "cells (whole brain)"),
    labels = scales::comma,
    expand = expansion(mult = c(0.03, 0.02))
  ) +

  scale_x_discrete(
    labels = c("Vehicle"            = "Veh",
               "Acute Morphine"     = "Mor",
               "Morphine-Dependent" = "MorDep",
               "Ro 64-6198"         = "Ro")
  ) +

  # Two-line ANOVA label, centered under each sex strip
  geom_text(data = sex_totals_p, aes(x = 2.5, y = Inf, label = lab),
            inherit.aes = FALSE, vjust = 1.2, lineheight = 0.95,
            size = FS_MARK_MM, family = FONT_FAMILY, color = "#555555") +

  facet_wrap(~ sex, nrow = 1) +
  coord_cartesian(clip = "off") +

  theme_classic(base_family = FONT_FAMILY, base_size = 8) +
  theme(
    strip.text        = element_text(size = FS_TITLE, face = "bold", margin = margin(b = 2)),
    strip.background  = element_blank(),
    axis.title.x      = element_blank(),
    axis.title.y      = element_text(size = FS_TITLE, margin = margin(r = 3)),
    axis.text.x       = element_text(size = FS_TICK, color = "black",
                                     angle = 30, hjust = 1, vjust = 1),
    axis.text.y       = element_text(size = FS_TICK, color = "black"),
    axis.line         = element_line(linewidth = 0.4, color = "#333333"),
    axis.ticks        = element_line(linewidth = 0.3, color = "#333333"),
    axis.ticks.length = unit(2, "pt"),
    panel.spacing     = unit(6, "pt"),
    legend.position   = "none",
    plot.margin       = margin(6, 8, 4, 4),
    plot.background   = element_rect(fill = "white", color = NA)
  )

save_panel(p_fig2b, file.path(output_dir, "Figure2_PanelB_SexTotals.tiff"),
           width = 66, height = 104)   # composite slot (new Fig 2, top-right; right-aligned to band in assembler)
cat("Fig 2 Panel B saved (sex-split whole-brain totals).\n")


# =============================================================================
# PANEL C - Sex-resolved exemplar grid: SNr / AVPV / PRE / PRM x Males / Females
# The former Fig 3B recipe, regions swapped to SNr/AVPV/PRE (the old sex-bar set)
# plus PRM (paramedian lobule), and transposed so regions run ACROSS the columns and Males/Females
# stack as the two rows - the same horizontal exemplar-strip look as the S1
# MBO/ASO/TRS/MEA panel, now split by sex. All 4 conditions; individual animals
# overlaid. Y is free per region but SHARED across the two sex rows within a
# region (so Males vs Females are directly comparable) - achieved with
# ggh4x::facet_grid2(independent = "y") plus identical per-region y-limits pinned
# into both sex cells. Within-sex FDR brackets from the BySex ANOVA output.
# =============================================================================

# Significance comparisons read from ANOVA output at runtime. Within each sex:
# Veh vs Mor (1-2), Veh vs MorDep (1-3), Veh vs Ro (1-4), Mor vs MorDep (2-3).
# Brackets are drawn only where sig_fdr != "ns".
fig2c_comparisons <- c("Veh vs Mor", "Veh vs MorDep", "Veh vs Ro", "Mor vs MorDep")

cat("Loading Fig 2 Panel C data (sex-resolved exemplar grid)...\n")
bsex_regions <- c("SNr","AVPV","PRE","PRM")
sex_levels   <- c("Males","Females")

panelB_df <- purrr::map_dfr(bsex_regions, function(acr) {
  bind_rows(
    load_and_normalize(acr, "Clustered Males",   "Males"),
    load_and_normalize(acr, "Clustered Females", "Females")
  ) %>% dplyr::mutate(region = acr)   # see note in bar_df: value is log10 density
}) %>%
  dplyr::mutate(condition = factor(condition, levels = condition_order),
                region    = factor(region, levels = bsex_regions),
                sex       = factor(sex,    levels = sex_levels))

# Per-region y-limits (shared across sexes); per-(region,sex) bracket positions
brkB_list <- list(); limB_list <- list()
for (acr in bsex_regions) {
  d_reg <- dplyr::filter(panelB_df, region == acr)
  ymax  <- max(d_reg$value, na.rm = TRUE); ymin <- min(d_reg$value, na.rm = TRUE)
  step  <- max((ymax - ymin) * 0.14, 0.13)
  n_sig_max <- 0
  for (sx in sex_levels) {
    brks <- Filter(function(b) b$sig != "ns",
                   build_brackets_sex(acr, sx, fig2c_comparisons))
    if (length(brks)) {
      for (i in seq_along(brks)) {
        b <- brks[[i]]
        brkB_list[[length(brkB_list) + 1]] <- data.frame(
          region = acr, sex = sx, x1 = b$x1, x2 = b$x2,
          xm = (b$x1 + b$x2) / 2, y = ymax + 0.10 + step * (i - 1),
          sig = b$sig, stringsAsFactors = FALSE)
      }
      n_sig_max <- max(n_sig_max, length(brks))
    }
  }
  ytop <- if (n_sig_max > 0) ymax + 0.10 + step * (n_sig_max - 1) + step * 0.85
  else ymax + 0.15
  limB_list[[length(limB_list) + 1]] <- data.frame(
    region = acr, ylo = ymin - (ymax - ymin) * 0.05, yhi = ytop)
}
brkB_df <- if (length(brkB_list)) do.call(rbind, brkB_list) else
  data.frame(region = character(), sex = character(), x1 = double(),
             x2 = double(), xm = double(), y = double(), sig = character())
limB_df <- do.call(rbind, limB_list)
brkB_df$region <- factor(brkB_df$region, levels = bsex_regions)
limB_df$region <- factor(limB_df$region, levels = bsex_regions)
if (nrow(brkB_df)) brkB_df$sex <- factor(brkB_df$sex, levels = sex_levels)

# Pin each region row's y-range in both sex columns; bracket tick marks
limB_long <- limB_df %>%
  tidyr::pivot_longer(c(ylo, yhi), values_to = "value") %>%
  tidyr::crossing(sex = factor(sex_levels, levels = sex_levels)) %>%
  dplyr::mutate(condition = factor("Vehicle", levels = condition_order))
brkB_tick <- if (nrow(brkB_df)) {
  tidyr::pivot_longer(brkB_df, c(x1, x2), values_to = "x") %>%
    dplyr::mutate(yend = y - 0.05)
} else brkB_df

# Per-facet x-axis baseline - facet_grid2 + theme_classic drops the axis line on
# interior facets, so draw it at each facet's lower limit in both sex rows.
baseB <- limB_df %>%
  tidyr::crossing(sex = factor(sex_levels, levels = sex_levels))

p_bsex <- ggplot(panelB_df, aes(x = condition, y = value)) +
  
  geom_beeswarm(aes(color = condition),
                size = 0.8, shape = 16, cex = 2.6, alpha = 0.80) +
  
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.25, linewidth = 0.4, color = "black") +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.40, color = "black", fill = NA, linewidth = 0.4) +
  
  geom_blank(data = limB_long, aes(x = condition, y = value)) +
  
  geom_segment(data = baseB, aes(x = 0.5, xend = 4.5, y = ylo, yend = ylo),
               inherit.aes = FALSE, linewidth = 0.4, color = "#333333") +
  
  geom_segment(data = brkB_df, aes(x = x1, xend = x2, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.28) +
  geom_segment(data = brkB_tick, aes(x = x, xend = x, y = y, yend = yend),
               inherit.aes = FALSE, linewidth = 0.28) +
  geom_text(data = brkB_df, aes(x = xm, y = y + 0.04, label = sig),
            inherit.aes = FALSE, family = FONT_FAMILY, size = FS_MARK_MM, vjust = 0) +

  scale_color_manual(values = condition_colors, guide = "none") +
  scale_x_discrete(labels = c("Vehicle"="Veh","Acute Morphine"="Mor",
                              "Morphine-Dependent"="MorDep","Ro 64-6198"="Ro")) +
  scale_y_continuous(name = expression(log[10](cFos^"+" ~ cells/mm^3)),
                     breaks = pretty_breaks(n = 4),
                     labels = function(x) sprintf("%.1f", x),   # uniform 1-decimal ticks across facets
                     expand = expansion(mult = c(0, 0.04))) +
  
  # Regions ACROSS the columns, Males/Females as the two rows. independent = "y"
  # gives each panel its own y-scale; because both sex cells of a region are
  # pinned (via limB_long) to that region's shared [ylo, yhi], the two rows of a
  # column match while different regions keep their own scale.
  ggh4x::facet_grid2(sex ~ region, scales = "free_y", independent = "y") +
  coord_cartesian(clip = "off") +

  theme_classic(base_family = FONT_FAMILY, base_size = 8) +
  theme(
    strip.text.x      = element_text(size = FS_TITLE, face = "bold", margin = margin(b = 2)),
    strip.text.y      = element_text(size = FS_TITLE, face = "bold", margin = margin(l = 2)),
    strip.background  = element_blank(),
    axis.title.x      = element_blank(),
    axis.title.y      = element_text(size = FS_TITLE, margin = margin(r = 2)),
    axis.text.x       = element_text(size = FS_TICK, color = "black",
                                     angle = 30, hjust = 1, vjust = 1),
    axis.text.y       = element_text(size = FS_TICK, color = "black"),
    axis.line.y       = element_line(linewidth = 0.4, color = "#333333"),
    axis.line.x       = element_blank(),
    axis.ticks        = element_line(linewidth = 0.3, color = "#333333"),
    axis.ticks.length = unit(2, "pt"),
    panel.spacing.x   = unit(4, "pt"),
    panel.spacing.y   = unit(5, "pt"),
    legend.position   = "none",
    plot.margin       = margin(4, 6, 2, 2),
    plot.background   = element_rect(fill = "white", color = NA)
  )

save_panel(p_bsex, file.path(output_dir, "Figure2_PanelC_SexExemplars.tiff"),
           width = 169, height = 70)   # composite slot (new Fig 2, full-width band, sex x region)
cat("Fig 2 Panel C saved (sex-resolved exemplar grid).\n")


# =============================================================================
# PANEL D - Representative light-sheet images (Veh / Mor / MorDep / Ro)
# Not produced by this script: a representative MBO image strip composited into
# the figure separately, placed full-width at the very bottom of Figure 2.
# 500 um scale bar. Kept here as a manifest note so the panel inventory in this
# file stays complete.
# =============================================================================


cat("\nAll new Figure 2 panels saved to:", output_dir, "\n")
