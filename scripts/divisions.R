# divisions.R
# =============================================================================
# CCFv3 anatomical divisions: the single definition of the division grouping
# and colour palette used by every figure script.
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
# This file is sourced by the figure scripts rather than run on its own, so that
# the division grouping and palette are identical across all figures.
#
#   division_order  : canonical legend order (13 divisions)
#   division_colors : 13 divisions -> hex   (main figures: 2, 3, 4, Fig 5 cartography)
#   division_abbr   : 13 full names -> short codes (for space-limited panels)
#   super_order     : canonical order (7 super-divisions)
#   super_map       : 13 divisions -> 7 super-divisions
#   super_colors    : 7 super-divisions -> hex  (network hairballs only; each colour
#                     is taken from a member division so the two levels stay coherent)
#   division_lookup : region acronym -> division   (from region_division_lookup.csv)
#   struct_of(x)    : acronyms -> 13-division      (NA -> "Other")
#   super_of(x)     : acronyms -> 7-super-division (NA -> "Other")
# =============================================================================

# Directory holding region_division_lookup.csv. Resolved in this order:
#   1. the CFOS_DATA_DIR environment variable, if set;
#   2. a DATA_DIR variable already defined by the script that sources this file
#      (so a caller working outside the repo layout does not have to set the
#      environment variable as well);
#   3. "data" -- the repo-relative default, for a run from the repository root.
.DIV_DIR <- local({
  d <- Sys.getenv("CFOS_DATA_DIR", "")
  if (!nzchar(d) && exists("DATA_DIR", envir = globalenv(), inherits = FALSE))
    d <- get("DATA_DIR", envir = globalenv())
  if (!nzchar(d)) "data" else d
})

division_order <- c(
  "Isocortex", "Olfactory", "Hippocampus", "Cortical Subplate",
  "Amygdala", "Striatum", "Pallidum/Septum", "Thalamus",
  "Hypothalamus", "Midbrain", "Pons", "Medulla", "Cerebellum")

division_colors <- c(
  "Isocortex"         = "#3949ab", "Olfactory"    = "#4dad4c",
  "Hippocampus"       = "#ff9800", "Cortical Subplate" = "#00838f",
  "Amygdala"          = "#9c27b0", "Striatum"     = "#ff5da2",  # kept distinct from Midbrain
  "Pallidum/Septum"   = "#f64234", "Thalamus"     = "#607d8b",
  "Hypothalamus"      = "#4cd0e9", "Midbrain"     = "#ec1c67",
  "Pons"              = "#7b5548", "Medulla"      = "#ba7517",
  "Cerebellum"        = "#fdd835")

division_abbr <- c(
  "Isocortex"         = "CTX",  "Olfactory"    = "OLF", "Hippocampus" = "HPF",
  "Cortical Subplate" = "CTXsp","Amygdala"     = "AMY", "Striatum"    = "STR",
  "Pallidum/Septum"   = "PAL",  "Thalamus"     = "TH",  "Hypothalamus"= "HY",
  "Midbrain"          = "MB",   "Pons"         = "P",   "Medulla"     = "MY",
  "Cerebellum"        = "CB")

super_order <- c(
  "Cerebral cortex", "Cerebral nuclei", "Thalamus", "Hypothalamus",
  "Midbrain", "Hindbrain", "Cerebellum")

super_map <- c(
  "Isocortex"         = "Cerebral cortex", "Olfactory"         = "Cerebral cortex",
  "Hippocampus"       = "Cerebral cortex", "Cortical Subplate" = "Cerebral cortex",
  "Amygdala"          = "Cerebral nuclei", "Striatum"          = "Cerebral nuclei",
  "Pallidum/Septum"   = "Cerebral nuclei",
  "Thalamus"          = "Thalamus",        "Hypothalamus"      = "Hypothalamus",
  "Midbrain"          = "Midbrain",
  "Pons"              = "Hindbrain",        "Medulla"          = "Hindbrain",
  "Cerebellum"        = "Cerebellum")

# Each super-colour is a representative member division's colour (all drawn from
# the 13 above), so the coarse hairballs read as a grouping of the same palette.
super_colors <- c(
  "Cerebral cortex" = "#4dad4c",  # Olfactory (cortex family)
  "Cerebral nuclei" = "#9c27b0",  # Amygdala  (nuclei family)
  "Thalamus"        = "#607d8b",  # exact
  "Hypothalamus"    = "#4cd0e9",  # exact
  "Midbrain"        = "#ec1c67",  # exact
  "Hindbrain"       = "#7b5548",  # Pons (hindbrain family)
  "Cerebellum"      = "#fdd835")  # exact

# Region -> division lookup (canonical 198-region CSV; single assignment source).
division_lookup <- local({
  f <- file.path(.DIV_DIR, "region_division_lookup.csv")
  if (!file.exists(f))
    stop("divisions.R: cannot find region_division_lookup.csv at '", f, "'.\n",
         "  Run from the repository root, or set the CFOS_DATA_DIR environment\n",
         "  variable (or define DATA_DIR before sourcing this file) to the\n",
         "  folder that holds it.")
  tab <- utils::read.csv(f, stringsAsFactors = FALSE)
  stats::setNames(tab$Structure, tab$Region)
})

struct_of <- function(regions) {
  s <- unname(division_lookup[regions]); s[is.na(s)] <- "Other"; s
}
super_of <- function(regions) {
  d <- struct_of(regions); s <- unname(super_map[d]); s[is.na(s)] <- "Other"; s
}
