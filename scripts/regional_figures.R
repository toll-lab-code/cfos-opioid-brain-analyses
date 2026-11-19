# regional_figures.R
# =============================================================================
# Figure 2 (all 198 regions x 4 conditions), Figure 3 (pooled regional
# detail, combined sexes) and Figure 4 (sex differences), numbered as in the
# published paper. All text is set to the journal's >= 8 pt floor via the FS_*
# constants below.
#
# Accompanies:
#   Martinez M, Ozawa A, Van Zant D, Thornberry J, Toll L.
#   Whole Brain Cellular Activation After Mu and NOP Receptor Agonism:
#   Shared Regional Targets with Divergent Network Consequences.
#   Submitted to Psychopharmacology, 2026.
#
# Code author:  Madeline Martinez
# Contact:      Lawrence Toll (corresponding author)
# Affiliation:  Toll Lab, Stiles-Nicholson Brain Institute,
#               Charles E. Schmidt College of Medicine,
#               Florida Atlantic University
# Repository:   https://github.com/toll-lab-code/cfos-opioid-brain-analyses
# Archived:     Zenodo; DOIs for this release are given in readme.md
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
# --- Figure 2 (comprehensive regional map; single panel, no assembly) ---------
# Figure 2:          All 198 analysed regions x 4 conditions. Left block, mean
#                    log10 density per condition; right block, Cohen's d vs
#                    vehicle with FDR marks. Drawn whole in base graphics at
#                    174 x ~179 mm; see the FIGURE 2 section at the end of this
#                    file. Answers R1 #3 and puts the pooled regional analysis
#                    in the main text (R3 major 1).
# --- Figure 4 (sex differences; PIL-assembled) -------------------------------
# Figure 4, Panel A: Diverging bar chart - males vs. females, Veh vs. Acute
#                    Morphine. Axis names the contrast (R1 #4); caption states
#                    that no region shows a condition-by-sex interaction (R3).
# Figure 4, Panel B: Whole-brain activation - Combined | Males | Females facets,
#                    all 4 conditions, omnibus ANOVA per stratum. The Combined
#                    facet carries the pooled whole-brain result.
#                    Caption gives the vehicle female-vs-male p.
# Figure 4, Panel C: Sex-resolved exemplar grid - SNr / AVPV / PRE / PRM (columns)
#                    x Males / Females (rows); per-region free y shared across sex;
#                    within-sex FDR brackets.
# Figure 3, Panel D: Representative light-sheet images (MBO) - built by
#                    panel_mbo_strip.py and composited separately.
# --- Figure 3 (pooled regional detail; PIL-assembled) ----------------------
# Figure 3, Panel A:   Lollipop - Veh vs Acute Morphine; FDR-significant regions,
#                    top 30 by Cohen's |d| (read live from the stats). Table S2.
# Figure 3, Panel B:   Significant-region count summary - uncorrected & FDR counts
#                    per drug condition (combined sex).
# Figure 3, Panel C:   Faceted exemplar bars - MBO, ASO, TRS, MEA; all 4 conditions;
#                    combined sex. (Pooled whole-brain activation is shown in
#                    Figure 4 Panel B, Combined facet.)
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
library(viridisLite)   # Figure 2 heat-map density scale
library(ggh4x)      # facet_grid2(): per-region free y with regions across, sex down

# ---- Font sizes (every text element >= 8 pt at final size) ----
# All figure text is floored at 8 pt on the printed page (174 mm full width,
# Psychopharmacology). Every panel is rendered at its final slot size and placed 1:1
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

# This first block builds the Figure 3 panels. output_dir is reassigned to
# figures/Figure4 for the Figure 4 section, and Figure 2 (at the end of this
# file) writes to figures/Figure2.
output_dir <- file.path(OUTPUT_DIR, "figures", "Figure3")
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
stat_global_file  <- file.path(OUTPUT_DIR, "Drug_Global_Activation.xlsx")

# ---- Whole-brain (global) activation statistics ------------------------------
# Computed once in regional_analysis.R (SUPPLEMENTARY 4) and READ here, never
# recomputed, so the whole-brain panels cannot drift from the statistics and
# every global number in the manuscript comes out of a deposited workbook.
#
# That section computes three measures, and THEY DO NOT AGREE:
#   total_counts   per-animal sum of raw counts   (condition effect weak)
#   offset_all     per-animal mean log10 density  (condition effect clear)
#   offset_common  as offset_all, common region set across all animals
# The raw total is carried by a handful of large structures; the offset weights
# every region equally. Whichever is set below is what the whole-brain panels
# plot AND annotate, and it must be the measure the Results text cites.
# Set to "offset_all" following the convention in this literature: the
# comparator studies report and analyse c-Fos DENSITY, log-transformed, and where
# they state a global figure it is a mean density rather than a summed count
# (Vasylieva et al. 2026 report global c-Fos as density; Kimbrough et al. 2020
# log10-transform explicitly "to bring brain regions with high numbers and lower
# numbers of Fos counts to a similar scale"). A raw summed whole-brain total is
# the one global measure none of the comparators tests, because it is carried by
# the largest structures.
GLOBAL_MEASURE <- "offset_all"

.global_pooled  <- suppressMessages(read_excel(stat_global_file, sheet = "Pooled"))
.global_bysex   <- suppressMessages(read_excel(stat_global_file, sheet = "By sex"))
.global_animals <- suppressMessages(read_excel(stat_global_file, sheet = "Per-animal values"))

# Contrast names as TukeyHSD writes them in the analysis script, keyed by the
# display condition labels used in the figures.
GLOBAL_CONTRAST <- list("Acute Morphine"     = "Morphine-Vehicle",
                        "Morphine-Dependent" = "Morphine-Dependent-Vehicle",
                        "Ro 64-6198"         = "Ro-Vehicle")

# p value for one measure / stratum / contrast. stratum is "Pooled", "Male" or
# "Female"; contrast is "OMNIBUS" or a GLOBAL_CONTRAST value.
global_stat <- function(measure, stratum, contrast = "OMNIBUS") {
  d <- if (stratum %in% c("Male", "Female")) .global_bysex else .global_pooled
  r <- d[d$measure == measure & grepl(paste0("^", stratum), d$stratum) &
           d$contrast == contrast, ]
  if (!nrow(r)) NA_real_ else as.numeric(r$p_value[1])
}

# Never print a small p as "0.00".
format_p <- function(p) ifelse(is.na(p), "NA",
                        ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p)))

# One wording for every density axis: "c-Fos+ density (log10 cells/mm3)" on the
# per-region panels, with "Mean" in front where the value is a mean over regions
# (Figure 4b) or over animals (Figure 2 colour bar).
DENS_Y_LAB <- expression("c-Fos"^"+" ~ "density (log"[10] ~ "cells/mm"^3 * ")")
GLOBAL_Y_LAB <- if (GLOBAL_MEASURE == "total_counts")
  expression("Total c-Fos"^"+" ~ "cells (whole brain)") else
  expression("Mean c-Fos"^"+" ~ "density (log"[10] ~ "cells/mm"^3 * ")")
GLOBAL_Y_FMT <- if (GLOBAL_MEASURE == "total_counts") scales::comma else
  function(x) sprintf("%.1f", x)          # one decimal, as on every density axis
GLOBAL_Y_BRK <- if (GLOBAL_MEASURE == "total_counts") waiver() else
  scales::breaks_width(0.5)               # 0.5 steps, as in Figure 3c
# second line of the two-line ANOVA label, in units of its own text height
ANOVA_P_VJUST <- 2.45


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
# PANEL A - Lollipop (Figure 3): FDR-significant regions (q<0.05, Vehicle vs. Acute
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

message(sprintf("Figure 3 Panel A: %d FDR-significant regions shown (top |d|)",
                nrow(figS1a_data)))

p_figS1a <- ggplot(figS1a_data, aes(x = abs_d, y = label, color = structure)) +
  
  geom_segment(aes(x = 0, xend = abs_d, yend = label),
               linewidth = 0.55, lineend = "round",
               key_glyph = key_square) +   # shared division-key swatch (divisions.R)
  
  # All plotted regions are FDR-significant; point filled with the structure color.
  geom_point(aes(fill = point_fill),
             size = 2.4, stroke = 0.6, shape = 21,
             show.legend = c(colour = FALSE, fill = FALSE)) +

  scale_color_manual(values = structure_colors, name = "CCFv3 division",
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
    panel.grid.major.x  = element_blank(),   # no gridlines in any panel
    panel.grid.major.y  = element_blank(),
    panel.grid.minor    = element_blank(),
    legend.key.width    = unit(KEY_SW_MM, "mm"),     # division-key spec (divisions.R)
    legend.key.height   = unit(KEY_PITCH_MM, "mm"),  # row pitch = key cell height
    legend.key.spacing.y = unit(0, "mm"),
    legend.text         = element_text(size = FS_LEGEND, margin = margin(l = KEY_GAP_MM, unit = "mm")),
    legend.title        = element_text(size = FS_LEGEND, face = "bold",
                                       margin = margin(b = KEY_TITLE_GAP_MM, unit = "mm")),
    legend.spacing.y    = unit(2, "pt"),
    legend.position     = "right",
    legend.box.spacing  = unit(2, "pt"),   # pull legend in close to the plot
    legend.background   = element_rect(fill = "white", color = NA),
    legend.margin       = margin(2, 2, 2, 2),
    plot.margin         = margin(8, 2, 6, 4),
    plot.background     = element_rect(fill = "white", color = NA),
    plot.caption        = element_text(size = FS_LEGEND, color = "black", hjust = 0)
  ) +
  
  labs(tag = NULL)

# The division key is saved as its own small panel and placed in the right
# column under panel b, which otherwise left ~60 mm of white space there; the
# lollipop gets the full 114 mm.
save_panel(p_figS1a + theme(legend.position = "none"),
           file.path(output_dir, "Figure3_PanelA.tiff"),
           width = 114, height = 110)   # composite slot (Figure 3 left column)
local({
  # build the grob on a real (showtext) device so text is measured in Arial
  grDevices::png(tf <- tempfile(fileext = ".png"), width = 50, height = 44,
                 units = "mm", res = 1200)
  g   <- ggplotGrob(p_figS1a)
  grDevices::dev.off(); unlink(tf)
  leg <- g$grobs[[grep("^guide-box", g$layout$name)[1]]]
  # top-left anchored, so the key lines up with panel b's left edge
  leg$vp <- grid::viewport(x = 0, y = 1, just = c(0, 1),
                           width = sum(leg$widths), height = sum(leg$heights))
  leg <- grid::grobTree(grid::rectGrob(gp = grid::gpar(fill = "white", col = NA)), leg)
  save_panel(leg, file.path(output_dir, "Figure3_PanelA_legend.tiff"),
             width = 50, height = 44)   # composite slot (Figure 3 right column, under b)
})
cat("Figure 3 Panel A saved (lollipop).\n")


# =============================================================================
# PANEL C (Figure 3) - Faceted exemplar bars: MBO, ASO, TRS, MEA; all 4 conditions;
#               combined sex. The MBO light-sheet images sit directly below as
#               Figure 3 Panel D (panel_mbo_strip.py).
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
  scale_y_continuous(name = DENS_Y_LAB,
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

save_panel(p_bars, file.path(output_dir, "Figure3_PanelC_bars.tiff"),
           width = 165, height = 46)   # composite slot (Figure 3 bottom row, full width)
                                      # 165 = 174 mm page - 2 x 4.5 mm margin.
                                      # Was 169, the band width at the old 178 mm page.
cat("Figure 3 Panel C (faceted exemplar bars) saved.\n")


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
message(sprintf("Figure 3 Panel B counts (unc/fdr): Mor %d/%d, MorDep %d/%d, Ro %d/%d",
                .counts["unc","Acute Morphine"], .counts["fdr","Acute Morphine"],
                .counts["unc","Morphine-Dependent"], .counts["fdr","Morphine-Dependent"],
                .counts["unc","Ro 64-6198"], .counts["fdr","Ro 64-6198"]))

p_counts <- ggplot(count_data, aes(condition, count, fill = condition, alpha = tier)) +
  
  geom_col(position = "identity", width = 0.68, key_glyph = key_square) +
  
  # Count labels: uncorrected total above each bar, FDR count atop its solid part
  geom_text(data = subset(count_data, tier == "Tukey p<0.05"),
            aes(label = count), vjust = -0.4, size = FS_MARK_MM, fontface = "bold",
            color = "black", alpha = 1,   # labels stay solid black, not tier-faded
            family = FONT_FAMILY, show.legend = FALSE) +
  geom_text(data = subset(count_data, tier == "FDR q<0.05" & count > 0),
            aes(label = count), vjust = -0.4, size = FS_MARK_MM, fontface = "bold",
            color = "black", family = FONT_FAMILY, show.legend = FALSE) +
  
  scale_fill_manual(values = condition_colors, guide = "none") +
  scale_alpha_manual(values = c("Tukey p<0.05" = 0.32, "FDR q<0.05" = 1.0),
                     labels = c("Tukey p<0.05" = expression("Tukey" ~ italic(p) < 0.05),
                                "FDR q<0.05"   = expression("FDR" ~ italic(q) < 0.05)),
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
    legend.text          = element_text(size = FS_LEGEND, margin = margin(l = KEY_GAP_MM, unit = "mm")),
    legend.key.width     = unit(KEY_SW_MM, "mm"),    # same swatch as the division keys
    legend.key.height    = unit(KEY_PITCH_MM, "mm"),
    legend.key.spacing.x = unit(3, "mm"),
    legend.key.spacing.y = unit(1, "pt"),
    legend.margin        = margin(0, 0, 1, 0),
    legend.background    = element_rect(fill = "white", color = NA),
    plot.margin          = margin(4, 4, 4, 4),
    plot.background      = element_rect(fill = "white", color = NA)
  )

save_panel(p_counts, file.path(output_dir, "Figure3_PanelB_CountSummary.tiff"),
           width = 50, height = 62)   # composite slot (Figure 3 right column, top; key below)
cat("Figure 3 Panel B (count summary) saved.\n")


# =============================================================================
# SHARED DATA - per-animal whole-brain values (used by Figure 4 Panel B)
# The pooled global-activation panel that used to live here has MOVED to
# Figure 4 Panel B as its "Combined" facet, so that the paper's first major
# finding appears in a main-text figure rather than only in the supplement.
# Only the data preparation remains, because Figure 4 Panel B builds on it.
#
# Per-animal whole-brain values and the omnibus one-way ANOVA are READ from
# Drug_Global_Activation.xlsx (regional_analysis.R, SUPPLEMENTARY 4), not
# recomputed here. Which measure is used is set by GLOBAL_MEASURE at the top.
#
# NB: on the raw whole-brain total the condition effect is weak (pooled omnibus
# p ~ 0.08) while on the per-animal global offset it is clear (p ~ 0.002). These
# are different measures of the same brains, not a discrepancy -- but the panel
# and the Results text must cite the SAME one.
# =============================================================================

# Per-animal values from the analysis workbook. `total` carries whichever
# measure GLOBAL_MEASURE names, so both whole-brain panels and their statistics
# are guaranteed to describe the same quantity.
total_data <- tibble(
  condition = factor(dplyr::recode(.global_animals$condition,
                                   "Morphine" = "Acute Morphine",
                                   "Ro"       = "Ro 64-6198"),
                     levels = condition_order),
  sex       = ifelse(.global_animals$sex == "Male", "Males", "Females"),
  animal    = .global_animals$animal,
  total     = as.numeric(.global_animals[[GLOBAL_MEASURE]]))
stopifnot(!all(is.na(total_data$total)), !any(is.na(total_data$condition)))

# Omnibus one-way ANOVA (combined sex), read from the analysis workbook.


cat("\nAll Figure 3 panels saved to:", output_dir, "\n")


# =============================================================================
# FIGURE 4 (sex differences) - redirect output to Figure4 subdirectory
# Panels: A diverging bars | B sex-split whole-brain totals |
#         C sex-resolved exemplar grid | D representative images (composited separately).
# total_data (whole-brain per-animal totals, with sex) was built in the Figure 3
# section above and is reused here for Panel B.
# =============================================================================
output_dir <- file.path(OUTPUT_DIR, "figures", "Figure4")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# PANEL A - Diverging bar chart: males vs. females, Veh vs. Acute Morphine
# Males extend right (+), females extend left (-)
#
# INCLUSION RULE: a region is plotted if it survives FDR (q<0.05) in EITHER sex.
# The same criterion is applied to both sexes. Under acute morphine no female
# region survives FDR, so this resolves to the 20 male FDR-significant regions.
# Values are read from Drug_Statistical_Results_BySex.xlsx ("Veh vs Mor - Males"
# / "- Fem").
#
# The previous rule additionally admitted regions significant at uncorrected
# p<0.05 in FEMALES ONLY. That was asymmetric: it admitted 8 such regions while
# excluding 57 male regions meeting the identical uncorrected standard, and the
# 8 admitted regions have comparable male and female effect sizes (e.g. LPO
# |d| 1.50 vs 1.50; NDB 1.73 vs 1.77), so they were not female-specific.
# Changed in response to reviewer comment; see Table S2 for all 198 regions.
#
# Opacity: full = FDR q<0.05, mid = uncorrected p<0.05 only, faint = ns.
# Tip dot = FDR-significant.
# =============================================================================

# Region set, effect sizes, and significance flags are read at run time from
# Drug_Statistical_Results_BySex.xlsx ("Veh vs Mor - Males" / "- Fem") so the
# panel cannot drift from the statistics.
.fig2a_m <- suppressMessages(read_excel(stat_bysex_file, sheet = "Veh vs Mor - Males")) %>%
  transmute(acronym,
            # Acronyms, as in Figure 2. Full region names consumed ~45 mm of this
              # panel, which no longer fits beside the three-facet Panel B in
              # the 165 mm band. Decoded in Figure 2's key and Table S1.
              label = acronym,
            d_m = abs(cohens_d), sig_m = p_tukey < 0.05, fdr_m = p_fdr < 0.05)
.fig2a_f <- suppressMessages(read_excel(stat_bysex_file, sheet = "Veh vs Mor - Fem")) %>%
  transmute(acronym,
            d_f = abs(cohens_d), sig_f = p_tukey < 0.05, fdr_f = p_fdr < 0.05)

panel_a_data <- .fig2a_m %>%
  left_join(.fig2a_f, by = "acronym") %>%
  mutate(across(c(sig_m, sig_f, fdr_m, fdr_f), ~ tidyr::replace_na(., FALSE)),
         d_f = tidyr::replace_na(d_f, 0)) %>%
  # inclusion rule: FDR in either sex (same criterion applied to both sexes)
  filter(fdr_m | fdr_f) %>%
  # Sort by male |d| descending - females mirror on left
  arrange(desc(d_m)) %>%
  mutate(
    label     = factor(label, levels = rev(label)),
    structure = factor(unname(div_map[acronym]), levels = names(structure_colors)),
    # Opacity tiers: full = FDR q<0.05, mid = uncorrected p<0.05 only, faint = ns.
    alpha_m   = dplyr::case_when(fdr_m ~ 1.0, sig_m ~ 0.55, TRUE ~ 0.22),
    alpha_f   = dplyr::case_when(fdr_f ~ 1.0, sig_f ~ 0.55, TRUE ~ 0.22)
  )

message(sprintf("Fig 2 Panel A: %d regions (%d male-FDR, %d female-FDR)",
                nrow(panel_a_data),
                sum(panel_a_data$fdr_m, na.rm = TRUE),
                sum(panel_a_data$fdr_f, na.rm = TRUE)))

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
    # Tier NAMES only. Spelling the thresholds out drove the legend to ~30 mm of
    # a 75 mm panel; "FDR q < 0.05 / Tukey p < 0.05 / n.s." is stated in the
    # manuscript figure legend instead.
    fdr ~ "FDR", sig ~ "Tukey", TRUE ~ "n.s."),
    levels = c("FDR", "Tukey", "n.s.")))

p_a <- ggplot(df_long_a,
              aes(x = value, y = label, fill = structure, alpha = sig_tier)) +
  
  # Bars
  geom_col(width = 0.75, color = NA, key_glyph = key_square) +   # shared division-key swatch (divisions.R)
  
  # Zero line
  geom_vline(xintercept = 0, linewidth = 0.5, color = "#333333") +
  
  # Solid dot at bar tip marks FDR-significant regions (q<0.05)
  geom_point(data = df_long_a %>% filter(fdr %in% TRUE),
             aes(x = value + sign(value) * 0.04, color = structure),
             size = 1.1, show.legend = FALSE) +
  
  # Full division names, as in every other division key (Figs 3a, 5, 6, S1, S4).
  scale_fill_manual(values  = structure_colors, name = "CCFv3 division",
                    guide = guide_legend(order = 1)) +
  scale_color_manual(values = structure_colors, name = "CCFv3 division", guide = "none") +
  scale_alpha_manual(
    name   = "vs. Vehicle",
    values = c("FDR" = 1.0, "Tukey" = 0.55, "n.s." = 0.22),
    labels = c(expression("FDR" ~ italic(q) < 0.05 ~ "\u25CF"),
               expression("Tukey" ~ italic(p) < 0.05), "n.s."),
    guide  = guide_legend(order = 2, override.aes = list(fill = "grey30"))) +
  
  # Symmetric-per-side limits derived from the data (see note on Figure 3 Panel A).
  # Females plot as -|d|, so the axis is mirrored: labels are absolute values.
  scale_x_continuous(
    name   = expression(paste("Cohen's |", italic(d),
                              "|, acute morphine vs. vehicle")),
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
           label = "\u2190 Female", hjust = 1,
           size = FS_MARK_MM, family = FONT_FAMILY, fontface = "bold",
           color = sex_colors[["Females"]]) +
  annotate("text", x = .p2a_lab_off, y = length(levels(df_long_a$label)) + 1.3,
           label = "Male \u2192", hjust = 0,
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
    panel.grid.major.x  = element_blank(),   # no gridlines in any panel
    panel.grid.major.y  = element_blank(),
    legend.key.width    = unit(KEY_SW_MM, "mm"),     # division-key spec (divisions.R)
    legend.key.height   = unit(KEY_PITCH_MM, "mm"),  # row pitch = key cell height
    legend.key.spacing.y = unit(0, "mm"),
    legend.text         = element_text(size = FS_LEGEND, margin = margin(l = KEY_GAP_MM, unit = "mm")),
    legend.title        = element_text(size = FS_LEGEND, face = "bold",
                                       margin = margin(b = KEY_TITLE_GAP_MM, unit = "mm")),
    legend.position     = "right",
    legend.box.spacing  = unit(3, "pt"),
    legend.margin       = margin(0, 0, 0, 2),
    legend.background   = element_rect(fill = "white", color = NA),
    # Top margin is dead space above the panel: the sex labels live INSIDE the
    # panel (y = n + 1.3, with matching expand), so nothing here needs reserving,
    # and the assembler sets the panel letter above the image. 18 pt left this
    # panel's content sitting ~4 mm lower than panel b's in the same row.
    plot.margin         = margin(7, 4, 6, 4),
    plot.background     = element_rect(fill = "white", color = NA),
    plot.caption        = element_text(size = FS_LEGEND, color = "black", hjust = 0),
      # hjust = 0 aligns the caption to the PANEL. This panel's y axis carries
      # full region names, so the panel starts ~45% across and a caption of this
      # length runs off the right edge. Align it to the whole plot instead.
      plot.caption.position = "plot"
  ) +
  
  labs(tag = NULL,
       # ASCII "-by-" rather than the multiplication sign: the glyph is absent
       # from some fallback fonts and silently renders as "..". The manuscript
       # text keeps the conventional "condition x sex" notation.
       caption = expression("No condition-by-sex interaction (all" ~ italic(q) > "0.71)"))

save_panel(p_a, file.path(output_dir, "Figure4_PanelA.tiff"),
           width = 75, height = 106)    # composite slot (Fig 4, top-left, diverging bars)
                                          # 75 + 2 gap + 88 (Panel B) = 165 mm band
cat("Fig 4 Panel A saved (diverging bars).\n")


# =============================================================================
# PANEL B - Sex-stratified whole-brain total cFos+ cells (Males | Females)
# Sex-split of the Figure 3 whole-brain panel: same per-animal totals, faceted by sex,
# with a per-sex one-way ANOVA omnibus p annotated in each facet (bare p, no
# "n.s." marker). Shared y across the two sex facets so Males and Females are
# directly comparable. Reuses total_data built in the Figure 3 section above.
# =============================================================================
# A Combined facet is included alongside Males and Females so that the POOLED
# global result -- the paper's first major finding, one-way ANOVA p = 0.002 --
# appears in a main-text figure. It previously existed only as a supplementary
# panel while the sex-split view led here, the same inversion R3 objected
# to for the regional analysis. Showing the pooled result and then its
# decomposition also makes the argument directly: the overall rise is carried by
# males, who start lower and travel further to the same endpoint.
total_data_sex <- dplyr::bind_rows(
    dplyr::mutate(total_data, sex = "Combined"), total_data) %>%
  dplyr::mutate(sex = factor(sex, levels = c("Combined", "Males", "Females")))

# Facet label -> stratum name as used in Drug_Global_Activation.xlsx.
.stratum_of <- function(sx)
  if (identical(sx, "Combined")) "Pooled" else sub("s$", "", sx)

# Per-sex omnibus one-way ANOVA, read from the analysis workbook. Two-line label
# so it fits centered under each sex strip without colliding with the y-axis.
sex_totals_p <- tibble(
  sex = factor(levels(total_data_sex$sex), levels = levels(total_data_sex$sex)),
  p   = vapply(levels(total_data_sex$sex),
               function(sx) global_stat(GLOBAL_MEASURE, .stratum_of(sx)),
               numeric(1))) %>%
  dplyr::mutate(lab   = "ANOVA",
                lab_p = ifelse(!is.na(p) & p < 0.001, 'italic(p) < "0.001"',
                               sprintf('italic(p) == "%s"', format_p(p))))
cat(sprintf("Fig 4 Panel B - whole-brain %s ANOVA by stratum: %s\n", GLOBAL_MEASURE,
            paste(sprintf("%s p=%.4g", sex_totals_p$sex, sex_totals_p$p),
                  collapse = "; ")))

# Significance brackets vs Vehicle, read from Drug_Global_Activation.xlsx: the
# per-sex Tukey HSD of each drug against Vehicle, adjusted within the family of
# pairwise condition comparisons -- the same convention as every other panel, so
# "Tukey p" means one thing throughout the paper. (This previously used Welch
# t-tests with a BH adjustment, a different family from the Tukey p values the
# text quotes.) Brackets are gated on the per-sex omnibus ANOVA: a sex only gets
# brackets if its ANOVA is p < 0.05. Stars: * p<0.05, ** p<0.01, *** p<0.001.
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
  padj <- vapply(whole_vs_veh, function(cc)
    global_stat(GLOBAL_MEASURE, .stratum_of(sx), GLOBAL_CONTRAST[[cc[[1]]]]),
    numeric(1))
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
# Headroom above the tallest bracket for the two-line ANOVA label. 1.4 steps
# sufficed when only one facet carried a bracket; the Combined facet stacks
# two (Mor and MorDep both differ from vehicle when sexes are pooled), so the
# upper bracket ran into the label. Scale the allowance with the number of
# brackets in the busiest facet rather than hard-coding it.
n_stack   <- if (nrow(totB_df)) max(table(totB_df$sex)) else 0
y_top_pin <- y_ceiling + b_step * (1.4 + 0.9 * max(0, n_stack - 1))
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
    name   = GLOBAL_Y_LAB,
    breaks = GLOBAL_Y_BRK,
    labels = GLOBAL_Y_FMT,
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
            inherit.aes = FALSE, vjust = 1.2,
            size = FS_MARK_MM, family = FONT_FAMILY, color = "black") +
  geom_text(data = sex_totals_p, aes(x = 2.5, y = Inf, label = lab_p),
            inherit.aes = FALSE, vjust = ANOVA_P_VJUST, parse = TRUE,
            size = FS_MARK_MM, family = FONT_FAMILY, color = "black") +

  facet_wrap(~ sex, nrow = 1,
             labeller = as_labeller(c(Combined = "Combined", Males = "Male", Females = "Female"))) +
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
    plot.background   = element_rect(fill = "white", color = NA),
    plot.caption      = element_text(size = FS_LEGEND, color = "black", hjust = 0),
    plot.caption.position = "plot"
  ) +

  labs(caption = expression("Vehicle, females vs. males:" ~ italic(p) == "0.29"))

save_panel(p_fig2b, file.path(output_dir, "Figure4_PanelB_SexTotals.tiff"),
           # Height is set so panel b's x-axis line lands at the SAME distance from
           # the panel top as panel a's, since the two sit side by side in row 1.
           # Measured on the rendered panels: a's axis sits 91.31 mm down, b's sat
           # 90.01 mm down at height 104. Everything below the axis (tick labels,
           # the vehicle sex-difference note) is fixed height, so added height goes
           # entirely into the plot region and moves the axis down 1:1 - hence
           # 104 + 1.30. Panel a is still the taller of the two, so the row height
           # and the page total are unchanged.
           width = 88, height = 105.3) # composite slot (Fig 4, top-right; x-axis aligned to panel a)
cat("Fig 4 Panel B saved (sex-split whole-brain totals).\n")


# =============================================================================
# PANEL C - Sex-resolved exemplar grid: SNr / AVPV / PRE / PRM x Males / Females
# The earlier sex-bar recipe, regions swapped to SNr/AVPV/PRE (the old sex-bar set)
# plus PRM (paramedian lobule), and transposed so regions run ACROSS the columns and Males/Females
# stack as the two rows - the same horizontal exemplar-strip look as the Figure 3
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

cat("Loading Fig 4 Panel C data (sex-resolved exemplar grid)...\n")
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
# The baseline has to span the panel's FULL x range, not the category range.
# scale_x_discrete expands by add = 0.6 either side, so with 4 conditions the
# panel runs 0.4 -> 4.6 while the segment was drawn 0.5 -> 4.5. That left a gap
# at each end: on the left the baseline stopped short of the y-axis line, so the
# two axes did not meet at the origin corner. Expansion and segment now derive
# from ONE constant, so they cannot drift apart again.
X_EXPAND_ADD <- 0.6
n_cond_B     <- length(unique(as.character(panelB_df$condition)))
x_axis_lo    <- 1 - X_EXPAND_ADD
x_axis_hi    <- n_cond_B + X_EXPAND_ADD

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
  
  geom_segment(data = baseB, aes(x = x_axis_lo, xend = x_axis_hi, y = ylo, yend = ylo),
               inherit.aes = FALSE, linewidth = 0.4, color = "#333333") +
  
  geom_segment(data = brkB_df, aes(x = x1, xend = x2, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.28) +
  geom_segment(data = brkB_tick, aes(x = x, xend = x, y = y, yend = yend),
               inherit.aes = FALSE, linewidth = 0.28) +
  geom_text(data = brkB_df, aes(x = xm, y = y + 0.04, label = sig),
            inherit.aes = FALSE, family = FONT_FAMILY, size = FS_MARK_MM, vjust = 0) +

  scale_color_manual(values = condition_colors, guide = "none") +
  scale_x_discrete(labels = c("Vehicle"="Veh","Acute Morphine"="Mor",
                              "Morphine-Dependent"="MorDep","Ro 64-6198"="Ro"),
                   expand = expansion(add = X_EXPAND_ADD)) +
  scale_y_continuous(name = DENS_Y_LAB,
                     breaks = pretty_breaks(n = 4),
                     labels = function(x) sprintf("%.1f", x),   # uniform 1-decimal ticks across facets
                     expand = expansion(mult = c(0, 0.04))) +
  
  # Regions ACROSS the columns, Males/Females as the two rows. independent = "y"
  # gives each panel its own y-scale; because both sex cells of a region are
  # pinned (via limB_long) to that region's shared [ylo, yhi], the two rows of a
  # column match while different regions keep their own scale.
  ggh4x::facet_grid2(sex ~ region, scales = "free_y", independent = "y",
                     labeller = labeller(sex = c(Males = "Male", Females = "Female"))) +
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

save_panel(p_bsex, file.path(output_dir, "Figure4_PanelC_SexExemplars.tiff"),
           width = 165, height = 70)   # composite slot (Fig 4, full-width band, sex x region)
                                      # 165 = 174 mm page - 2 x 4.5 mm margin.
                                      # Was 169, the band width at the old 178 mm page.
cat("Fig 4 Panel C saved (sex-resolved exemplar grid).\n")


# =============================================================================
# (Figure 4 has three panels. The representative MBO image strip that was its
# panel d is now Figure 3 Panel D, below the pooled MBO bars; see above.)
# Not produced by this script.
# 500 um scale bar. Kept here as a manifest note so the panel inventory in this
# file stays complete.
# =============================================================================


cat("\nAll new Figure 4 panels saved to:", output_dir, "\n")

# =============================================================================
# FIGURE 2 - Comprehensive regional activation across all 198 analysed regions
# -----------------------------------------------------------------------------
# WHY THIS FIGURE EXISTS
# ----------------------
# Every other regional panel in the paper shows a SELECTED subset of regions:
# the top 30 by effect size (Figure 3a), the 20 reaching FDR in either sex
# (Fig 4a), or four exemplars (Fig 4c). This figure shows all 198 analysed
# regions in every condition, with nothing selected on significance.
#
# WHY TWO COLOUR SCALES
# ---------------------
# Absolute activation varies far more BETWEEN regions (mean log10 density spans
# roughly 2.2 to 4.5 across regions) than BETWEEN conditions within a region
# (condition medians all fall between about 3.85 and 4.07). A single absolute
# scale is therefore dominated by regional identity and renders the treatment
# effects nearly invisible. The left block gives the absolute level per
# condition on a sequential scale; the right block gives the effect size versus
# vehicle on a diverging scale, which is where the treatment signal is legible.
# Both are needed: the absolute block answers "what is the activation level in
# each condition", the effect block answers "what did the drug change".
#
# DRAWN IN BASE GRAPHICS, at the final placement size, with user units equal to
# MILLIMETRES and device pointsize = 8. Text is therefore absolute 8 pt on the
# printed page and the geometry below is literal millimetres. The whole panel is
# one device - no tiling, no PIL assembly, and nothing is ever rescaled
# (rescaling a rendered panel shrinks fonts and dots unevenly; see FIGURE_SPEC
# section 5).
#
# Drawn in BASE GRAPHICS at the final placement size, with user units equal to
# MILLIMETRES and device pointsize = 8, so text is absolute 8 pt on the printed
# page. Wrapped in local() so its many short-named working variables (D, E, Q,
# n, ord, ...) cannot collide with anything else in this script.
# =============================================================================
local({

  # Inputs and output directory for this figure. output_dir is deliberately
  # re-bound inside local() so it shadows the value the earlier figure sections
  # left behind rather than writing Figure 2 into their folder.
  norm_file    <- file.path(OUTPUT_DIR, "normalized_data_combined.csv")
  primary_file <- file.path(OUTPUT_DIR, "Drug_Statistical_Results_Primary.xlsx")
  output_dir   <- file.path(OUTPUT_DIR, "figures", "Figure2")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)



# ---- Journal geometry -------------------------------------------------------
# Psychopharmacology (large-sized journal): 174 mm full width, 234 mm max height,
# all figure text 8-12 pt. See claude/figure_spec_psychopharmacology.md.
FIG_W_MM    <- 174
MARGIN_MM   <- 4.5     # outer left/right margin (FIGURE_SPEC section 4)
# 4 blocks of 50 rows, not 3 of 66. A 3-block figure is 222 mm tall, which
# leaves under 18 mm of a ~240 mm text block for the caption - about four lines,
# where this figure needs roughly 180 words to identify every element as the
# journal requires. Four narrower blocks bring the figure to ~177 mm and leave
# ample room for the caption on the same page.
N_BLOCKS    <- 4       # side-by-side column blocks the 198 regions are split into
ROW_MM      <- 2.85    # per-region row height; 8 pt needs >= ~2.82 mm
CELL_W_MM   <- 3.7     # heat-map cell width
LAB_W_MM    <- 8.5     # acronym label column (longest acronym is 5 characters)
BAND_W_MM   <- 1.2     # division colour band
GAP_LAB_MM  <- 0.8     # band -> first cell
GAP_MID_MM  <- 1.0     # density block -> effect block
GUTTER_MM   <- 4.0     # between column blocks
# Header band: group label, rule, then the rotated condition labels. The labels
# hang UPWARD from ty = HEAD_MM - 1.2, so the band has to be long enough for the
# longest one. "MorDep" is 3.556 em = 10.03 mm at 8 pt Arial, which at the old
# 16.0 put its ascender at 4.77 mm - just above the rule at gy + 1.9 = 4.9 mm,
# so the label collided with it. 17.5 leaves 1.37 mm of clearance. A runtime
# guard below re-checks this from the rendered string widths.
HEAD_MM     <- 17.5    # header band: group label + rotated condition labels
# Footer built from named components rather than one magic number, so the gap
# above the colour bars and the margin below the key can be tuned independently.
GAP_MAP_BAR <- 6.0     # heat-map bottom -> colour-bar title (keeps the bars clear)
BARTITLE_MM <- 4.2     # colour-bar title band (bold title + clearance above the bar, as in Fig 6)
BAR_H_MM    <- 2.41    # colour-bar thickness - measured off the cartography bars
BAR_W_MM    <- 57.2    # colour-bar length - ditto, so Figs 2, 4, 5 carry one bar size
TICK_MM     <- 0.9     # colour-bar tick length (matches the Fig 6 bars)
TICKLAB_MM  <- 2.2     # tick -> value centre (matches the Fig 6 bars)
GAP_BAR_KEY <- 6.2     # colour bars -> division key (clears the tick values)
KEY_MM      <- 10.2    # division key band: title + two rows of full names (divisions.R spec)
BOTTOM_MM   <- 1.2     # trailing margin below the key
BAR_MM      <- GAP_MAP_BAR + BARTITLE_MM + BAR_H_MM + GAP_BAR_KEY + KEY_MM + BOTTOM_MM
PT          <- 8       # device pointsize; every label is 8 pt

DENS_CONDS <- c("Vehicle", "Morphine", "Chronic Morphine", "Ro")
DENS_LABS  <- c("Veh", "Mor", "MorDep", "Ro")
EFF_SHEETS <- c("Veh vs Mor", "Veh vs MorDep", "Veh vs Ro")
EFF_LABS   <- c("Mor", "MorDep", "Ro")
ROW_ORDER  <- "effect"   # "effect" = by acute-morphine d within division; "alpha" = by acronym

# ---- Load per-condition mean log10 density ----------------------------------
if (!file.exists(norm_file))
  stop("cannot find ", norm_file, " - run regional_analysis.R first")
nd <- utils::read.csv(norm_file, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("name", "acronym") %in% names(nd)))

animal_cols <- setdiff(names(nd), c("name", "acronym"))
cond_of_col <- sub(" - .*$", "", animal_cols)
if (!all(DENS_CONDS %in% unique(cond_of_col)))
  stop("condition columns not found in ", basename(norm_file), ": expected ",
       paste(DENS_CONDS, collapse = ", "), "; saw ",
       paste(unique(cond_of_col), collapse = ", "))

dens <- sapply(DENS_CONDS, function(cd) {
  cols <- animal_cols[cond_of_col == cd]
  apply(nd[, cols, drop = FALSE], 1, function(v) mean(as.numeric(v), na.rm = TRUE))
})
rownames(dens) <- nd$acronym

# ---- Load effect sizes and FDR q-values -------------------------------------
eff <- matrix(NA_real_, nrow = nrow(nd), ncol = length(EFF_SHEETS),
              dimnames = list(nd$acronym, EFF_LABS))
qv  <- eff
for (k in seq_along(EFF_SHEETS)) {
  s <- suppressMessages(as.data.frame(read_excel(primary_file, sheet = EFF_SHEETS[k])))
  i <- match(nd$acronym, s$acronym)
  if (anyNA(i))
    stop("regions missing from sheet '", EFF_SHEETS[k], "': ",
         paste(utils::head(nd$acronym[is.na(i)]), collapse = ", "))
  eff[, k] <- as.numeric(s$cohens_d[i])
  qv[,  k] <- as.numeric(s$p_fdr[i])
}

# SELF-CHECK: the density block and the effect block must describe the same
# regions in the same order, or the two halves of a row belong to different
# structures. Stop rather than draw a misaligned figure.
if (!identical(rownames(dens), rownames(eff)))
  stop("density and effect matrices are not aligned by region")

# ---- Row order: division, then effect size (or acronym) within division -----
divn <- struct_of(nd$acronym)
divn[!divn %in% division_order] <- "Other"
div_f <- factor(divn, levels = c(division_order, "Other"))
ord <- if (ROW_ORDER == "effect") {
  order(div_f, -eff[, "Mor"], nd$acronym)
} else {
  order(div_f, nd$acronym)
}

acr <- nd$acronym[ord]; dvv <- as.character(div_f)[ord]
D <- dens[ord, , drop = FALSE]; E <- eff[ord, , drop = FALSE]; Q <- qv[ord, , drop = FALSE]
n <- length(acr)

# ---- Colour scales ----------------------------------------------------------
NCOL_RAMP <- 256
# Both scales are clipped to robust limits rather than the full range. A handful
# of very low-density regions (periventricular hypothalamic nuclei) stretch the
# absolute range so far that ~90% of cells would occupy a narrow slice of the
# ramp; likewise almost every effect is positive, so full-range symmetric limits
# leave the bulk of the effect cells nearly white. Clipping is stated on the
# colour bars, and values beyond the limits saturate rather than being hidden.
CLIP_Q <- 0.98
dens_lim <- unname(stats::quantile(D, c(1 - CLIP_Q, CLIP_Q), na.rm = TRUE))
dens_lim <- c(floor(dens_lim[1] * 10) / 10, ceiling(dens_lim[2] * 10) / 10)
dens_pal <- viridisLite::viridis(NCOL_RAMP)

eff_max <- ceiling(unname(stats::quantile(abs(E), CLIP_Q, na.rm = TRUE)) * 2) / 2
eff_lim <- c(-eff_max, eff_max)
n_clip  <- sum(D < dens_lim[1] | D > dens_lim[2], na.rm = TRUE) +
           sum(abs(E) > eff_max, na.rm = TRUE)
# ColorBrewer RdBu, identical to the `bwr` ramp in cartography_figures.R (WMDz)
# and .DIV_ANCHORS in coactivation_figures.R (co-activation change). The paper
# uses one convention throughout: viridis for unsigned magnitudes, this ramp for
# signed quantities centred on zero.
eff_pal <- grDevices::colorRampPalette(c(
  "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7",
  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(NCOL_RAMP)

map_col <- function(x, lim, pal) {
  i <- round((x - lim[1]) / diff(lim) * (length(pal) - 1)) + 1
  i[i < 1] <- 1; i[i > length(pal)] <- length(pal)
  out <- pal[i]; out[is.na(x)] <- "#ffffff"; out
}
# luminance-aware mark colour, so the FDR asterisk stays legible on any cell
mark_col <- function(hex) {
  rgbv <- grDevices::col2rgb(hex) / 255
  lum  <- 0.2126 * rgbv[1, ] + 0.7152 * rgbv[2, ] + 0.0722 * rgbv[3, ]
  ifelse(lum < 0.55, "white", "black")
}

# ---- Geometry ---------------------------------------------------------------
per_block   <- ceiling(n / N_BLOCKS)
block_w     <- LAB_W_MM + BAND_W_MM + GAP_LAB_MM +
               length(DENS_CONDS) * CELL_W_MM + GAP_MID_MM +
               length(EFF_LABS) * CELL_W_MM
content_w   <- N_BLOCKS * block_w + (N_BLOCKS - 1) * GUTTER_MM
FIG_H_MM    <- HEAD_MM + per_block * ROW_MM + BAR_MM

if (content_w > FIG_W_MM - 2 * MARGIN_MM)
  stop(sprintf("blocks are %.1f mm wide but only %.1f mm is available - reduce CELL_W_MM or N_BLOCKS",
               content_w, FIG_W_MM - 2 * MARGIN_MM))
if (FIG_H_MM > 234)
  stop(sprintf("figure would be %.1f mm tall, over the 234 mm journal maximum - raise N_BLOCKS",
               FIG_H_MM))
x0 <- (FIG_W_MM - content_w) / 2      # centre the block band in the 174 mm width

cat(sprintf("Figure 2: %d regions, %d blocks of <=%d rows, %.1f x %.1f mm\n",
            n, N_BLOCKS, per_block, FIG_W_MM, FIG_H_MM))
cat(sprintf("  density scale %.1f to %.1f log10 cells/mm3 | effect scale %+.1f to %+.1f Cohen's d\n",
            dens_lim[1], dens_lim[2], eff_lim[1], eff_lim[2]))
cat(sprintf("  %d of %d cells saturate at the clipped limits (%.1f%%)\n",
            n_clip, length(D) + length(E), 100 * n_clip / (length(D) + length(E))))

# ---- Draw -------------------------------------------------------------------
draw_panel <- function() {
  op <- par(mai = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  on.exit(par(op), add = TRUE)
  plot.new()
  plot.window(xlim = c(0, FIG_W_MM), ylim = c(FIG_H_MM, 0))   # y grows downward

  for (b in seq_len(N_BLOCKS)) {
    idx <- ((b - 1) * per_block + 1):min(b * per_block, n)
    bx  <- x0 + (b - 1) * (block_w + GUTTER_MM)
    dx0 <- bx + LAB_W_MM + BAND_W_MM + GAP_LAB_MM
    ex0 <- dx0 + length(DENS_CONDS) * CELL_W_MM + GAP_MID_MM

    # column headers: a group label over each sub-block, then rotated condition
    # labels. HEAD_MM must clear the longest rotated label ("MorDep").
    gy <- 3.0
    dens_mid <- dx0 + length(DENS_CONDS) * CELL_W_MM / 2
    eff_mid  <- ex0 + length(EFF_LABS)  * CELL_W_MM / 2
    text(dens_mid, gy, "Density", adj = c(0.5, 0.5), cex = 1, font = 2)
    text(eff_mid,  gy, "Effect",  adj = c(0.5, 0.5), cex = 1, font = 2)
    segments(dx0, gy + 1.9, dx0 + length(DENS_CONDS) * CELL_W_MM, gy + 1.9,
             col = "grey45", lwd = 0.4)
    segments(ex0, gy + 1.9, ex0 + length(EFF_LABS) * CELL_W_MM, gy + 1.9,
             col = "grey45", lwd = 0.4)
    ty <- HEAD_MM - 1.2
    # The rotated labels hang upward from ty. If a longer condition label is
    # ever added, this is where it starts overlapping the rule above it.
    if (b == 1) {
      .clear <- (ty - max(strwidth(c(DENS_LABS, EFF_LABS)))) - (gy + 1.9)
      if (.clear < 0.5)
        warning(sprintf(paste("Figure 2 header: only %.2f mm between the block rule",
                              "and the longest rotated condition label; increase HEAD_MM."),
                        .clear))
    }
    for (j in seq_along(DENS_LABS))
      text(dx0 + (j - 0.5) * CELL_W_MM, ty, DENS_LABS[j], srt = 90, adj = c(0, 0.5), cex = 1)
    for (j in seq_along(EFF_LABS))
      text(ex0 + (j - 0.5) * CELL_W_MM, ty, EFF_LABS[j], srt = 90, adj = c(0, 0.5), cex = 1)

    for (r in seq_along(idx)) {
      i  <- idx[r]
      yt <- HEAD_MM + (r - 1) * ROW_MM
      yb <- yt + ROW_MM
      # acronym, right-aligned against the division band
      text(bx + LAB_W_MM - 0.4, (yt + yb) / 2, acr[i], adj = c(1, 0.5), cex = 1)
      # division band
      rect(bx + LAB_W_MM, yt, bx + LAB_W_MM + BAND_W_MM, yb,
           col = if (dvv[i] %in% names(division_colors)) division_colors[[dvv[i]]] else "#bdbdbd",
           border = NA)
      # absolute density cells
      for (j in seq_along(DENS_CONDS)) {
        cc <- map_col(D[i, j], dens_lim, dens_pal)
        rect(dx0 + (j - 1) * CELL_W_MM, yt, dx0 + j * CELL_W_MM, yb,
             col = cc, border = "white", lwd = 0.2)
      }
      # effect-size cells, asterisk where the contrast survives FDR
      for (j in seq_along(EFF_LABS)) {
        cc <- map_col(E[i, j], eff_lim, eff_pal)
        rect(ex0 + (j - 1) * CELL_W_MM, yt, ex0 + j * CELL_W_MM, yb,
             col = cc, border = "white", lwd = 0.2)
        if (!is.na(Q[i, j]) && Q[i, j] < 0.05)
          text(ex0 + (j - 0.5) * CELL_W_MM, (yt + yb) / 2 + 0.25, "*",
               col = "black", cex = 1, adj = c(0.5, 0.5))
      }
    }
  }

  # ---- colour bars ----------------------------------------------------------
  map_bottom <- HEAD_MM + per_block * ROW_MM
  bar_y  <- map_bottom + GAP_MAP_BAR + BARTITLE_MM
  bar_h  <- BAR_H_MM
  bar_w  <- BAR_W_MM
  draw_bar <- function(bx, pal, lim, title, fmt, n_tick) {
    xs <- seq(bx, bx + bar_w, length.out = length(pal) + 1)
    for (k in seq_along(pal))
      rect(xs[k], bar_y, xs[k + 1], bar_y + bar_h, col = pal[k], border = NA)
    rect(bx, bar_y, bx + bar_w, bar_y + bar_h, border = "grey30", lwd = 0.3)
    # axis line along the bottom edge, as graphics::axis() draws it under the
    # cartography colour bars (Figure 6), so all colour bars share one idiom
    segments(bx, bar_y + bar_h, bx + bar_w, bar_y + bar_h, col = "black", lwd = 1)
    # NOTE: the y axis runs downward (ylim = c(FIG_H_MM, 0)), which inverts the
    # vertical sense of adj. Centre the title and offset it explicitly instead
    # of relying on adj[2], or it anchors into the bar rather than above it.
    # font = 2 is ignored by plotmath, so the titles carry bold() themselves.
    # Centred 2.3 mm above the bar top: clears the subscript/superscript.
    text(bx + bar_w / 2, bar_y - 2.3, title, adj = c(0.5, 0.5), cex = 1)
    # Ticks and values BELOW the bar, matching the colour bars in Figures 4 and
    # 5. These bars previously carried only their two end values, set beside the
    # bar - a third legend idiom in a three-figure set that otherwise shares one.
    # Ticks span the full range, so the end values are still shown.
    tk <- seq(lim[1], lim[2], length.out = n_tick)
    tx <- bx + (tk - lim[1]) / diff(lim) * bar_w
    segments(tx, bar_y + bar_h, tx, bar_y + bar_h + TICK_MM, col = "black", lwd = 1)
    lab <- sprintf(fmt, tk)
    lab[tk == 0] <- "0"          # a signed or decimal zero reads oddly; plain 0
    lab <- sub("^-", "\u2212", lab)   # true minus sign, as on every other axis
    text(tx, bar_y + bar_h + TICK_MM + TICKLAB_MM, lab, adj = c(0.5, 0.5), cex = 1)
  }
  # Lay the two bars and the FDR note out from content_w rather than from fixed
  # offsets: at 57.2 mm each they no longer fit the gaps that suited a 34 mm bar.
  NOTE_W <- 24
  slack  <- content_w - 2 * bar_w - NOTE_W
  gapb   <- max(6, slack / 3)
  bx_1   <- x0 + max(0, slack - 2 * gapb) / 2
  bx_2   <- bx_1 + bar_w + gapb
  draw_bar(bx_1, dens_pal, dens_lim,
           expression(bold("Mean c-Fos"^"+" * " density (log"[10] * " cells/mm"^3 * ")")), "%.1f", 4)
  draw_bar(bx_2, eff_pal, eff_lim,
           expression(bold("Cohen's ") * bolditalic(d) * bold(" vs. vehicle")), "%.1f", 5)
  # text() with adj[2] = 0.5 centres on the glyph bounding box INCLUDING
  # descenders. "FDR q < 0.05" has a descender and the bar end labels ("2.9",
  # "+2.0") do not, so centring both at the same y leaves this string sitting
  # visibly high. Nudge it down by half a descender so the footer reads level.
  DESC_MM <- 0.35
  text(bx_2 + bar_w + 4, bar_y + bar_h / 2 + DESC_MM,
       expression("*" ~ "FDR" ~ italic(q) < 0.05),
       adj = c(0, 0.5), cex = 1)   # level with the bar, a fixed 4 mm after it

  # division key: names the 13 colour bands, which are otherwise unlabelled.
  # Full names under a "CCFv3 division" title, in a grid of rows, drawn to the
  # shared division-key spec in divisions.R so it matches every other key (the
  # same grid as Figure S1's key). Previously one row of abbreviations.
  keys <- division_order
  sw   <- KEY_SW_MM
  ty   <- bar_y + bar_h + GAP_BAR_KEY + 1.4              # title centre
  text(x0, ty, "CCFv3 division", adj = c(0, 0.5), cex = 1, font = 2)
  ncol <- length(KEY_COL_X_MM)                          # fixed grid shared with Fig S1
  colx <- x0 + KEY_COL_X_MM
  for (k in seq_along(keys)) {
    rr <- (k - 1) %/% ncol; cc <- (k - 1) %% ncol + 1
    ky <- ty + KEY_TITLE_MM + rr * KEY_PITCH_MM
    rect(colx[cc], ky - sw / 2, colx[cc] + sw, ky + sw / 2,
         col = division_colors[[keys[k]]], border = NA)
    text(colx[cc] + sw + KEY_GAP_MM, ky, keys[k], adj = c(0, 0.5), cex = 1)
  }
  invisible(NULL)
}

# ---- Devices: 1200 dpi TIFF (LZW) + matching PNG ----------------------------
w_in <- FIG_W_MM / 25.4; h_in <- FIG_H_MM / 25.4
fam  <- FONT_FAMILY
tiff(file.path(output_dir, "Figure2_regional_heatmap.tiff"),
     width = w_in, height = h_in, units = "in", res = 1200,
     compression = "lzw", type = "cairo", family = fam, pointsize = PT)
draw_panel(); dev.off()
png(file.path(output_dir, "Figure2_regional_heatmap.png"),
    width = w_in, height = h_in, units = "in", res = 1200,
    type = "cairo", family = fam, pointsize = PT)
draw_panel(); dev.off()

# ---- Deposited row order, so the figure and the tables cannot drift ---------
utils::write.csv(
  data.frame(row = seq_len(n), block = rep(seq_len(N_BLOCKS), each = per_block)[seq_len(n)],
             acronym = acr, region_name = nd$name[ord], division = dvv,
             # Prefix each block before writing. D, E and Q all carry the same
             # condition names, so with check.names = FALSE the header came out as
             # ...,Ro,Mor,MorDep,Ro,Mor,MorDep,Ro - "Ro" three times. pandas,
             # csv.DictReader and read.csv() all silently keep only the LAST column
             # of a repeated name, so the deposited file read wrong without ever
             # erroring (an acute-morphine q-value looks like a density).
             `colnames<-`(round(D, 4), paste0("dens_", DENS_LABS)),
             `colnames<-`(round(E, 4), paste0("d_",    EFF_LABS)),
             `colnames<-`(round(Q, 5), paste0("q_",    EFF_LABS)),
             check.names = FALSE),
  file.path(output_dir, "Figure2_row_order.csv"), row.names = FALSE)

cat("  wrote Figure2_regional_heatmap.tiff/.png and Figure2_row_order.csv to\n    ",
    output_dir, "\n", sep = "")

})
