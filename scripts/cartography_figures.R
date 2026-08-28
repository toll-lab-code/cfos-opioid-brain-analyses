#!/usr/bin/env Rscript
# cartography_figures.R
# =============================================================================
# Figure 4, Panels A and B: cartographic node metrics mapped onto the atlas.
#   Panel A : participation coefficient (PC) coronal atlas strips + colorbar
#   Panel B : within-module degree z-score (WMDz) coronal atlas strips + colorbar
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
#   Rscript scripts/cartography_figures.R
# Input and output locations are set in the Config section below; see README.
# =============================================================================
#
# Run network_analysis.R first: this script reads the Node_Roles sheet of
# 03_Modules_and_Roles.xlsx, which that step produces.
#
# ADDITIONAL REQUIREMENTS beyond the other scripts:
#   * Scalable Brain Atlas (SBA) coronal SVG data; see the Config section and
#     the README for the expected folder contents and where to obtain it.
#   * System libraries for rsvg and magick:
#       Debian/Ubuntu  librsvg2-dev, libmagick++-dev
#       macOS/Homebrew librsvg, imagemagick
# =============================================================================

# ---- Packages ----------------------------------------------------------------
required_pkgs <- c("jsonlite", "xml2", "viridisLite", "scales", "stringr",
                   "rsvg", "magick", "readxl")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages. Please run:\n  install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}
suppressPackageStartupMessages({
  library(jsonlite); library(xml2); library(viridisLite); library(scales)
  library(stringr);  library(rsvg); library(magick); library(readxl)
})

# ---- Config (shared by every panel) ------------------------------------------
# Paths are relative to the repository root. Run from the repo root, e.g.:
#   Rscript scripts/cartography_figures.R
# Every location is overridable with an environment variable. CFOS_SBA_DIR is
# the one most likely to need it: the Scalable Brain Atlas data is third-party,
# is not redistributed here, and commonly lives outside the repository.
DATA_DIR    <- path.expand(Sys.getenv("CFOS_DATA_DIR",   "data"))
OUTPUT_ROOT <- path.expand(Sys.getenv("CFOS_OUTPUT_DIR", "results"))
SCRIPT_DIR  <- path.expand(Sys.getenv("CFOS_SCRIPT_DIR", "scripts"))
# Third-party atlas folder: rgb2acr.json, acr2parent.json, acr2full.json, coronal_svg/
SBA_DIR     <- path.expand(Sys.getenv("CFOS_SBA_DIR", file.path(DATA_DIR, "sba_data")))
SHEET_FILE  <- file.path(OUTPUT_ROOT, "03_Modules_and_Roles.xlsx")   # network-pipeline output (Node_Roles sheet)
OUTPUT_DIR  <- file.path(OUTPUT_ROOT, "figures", "Figure4")          # strip + colorbar PNGs are written here
SHEET_NAME  <- "Node_Roles"
FONT        <- Sys.getenv("CFOS_FONT", "Arial")

# Graphics device. "cairo" writes plain sRGB with NO embedded colour profile,
# which is what regional_figures.R and coactivation_figures.R already produce
# and what the PIL assembler expects. The macOS quartz device (the default for
# png()/tiff() there) embeds an Apple "Generic RGB" profile AND stores shifted
# raw pixel values -- e.g. #1B9E77 is written as #1D8F64 -- so a division
# colour composited from a quartz panel does not match the same colour taken
# from the other panels, because PIL pastes raw pixels and ignores the profile.
# Set CFOS_DEVICE=quartz only if cairo cannot resolve the font on your machine.
dev_type <- Sys.getenv("CFOS_DEVICE", "cairo")

# Condition order. The atlas writes one file per condition, so this is
# order-independent here; it is kept for consistency with the other figures.
TREATMENTS <- c("Ro 64-6198", "Vehicle", "Acute Morphine", "Morphine-Dependent")

# #############################################################################
# ## SECTION A / B -- Regional atlas maps (PC, WMDz) + colorbars              ##
# #############################################################################
# =====================================================================
# Regional Atlas Visualization  (R port of generate_regional_atlas.py)
# =====================================================================
# Colors Scalable Brain Atlas (Allen Mouse Brain CCFv3, template ABA_v3) coronal
# SVG regions by per-region participation coefficient (PC) or
# within-module degree z-score (WMDz) for each treatment group.
#
# For each SBA region polygon, the metric value is resolved by:
#   1. exact acronym match in the data
#   2. else: descendants of the SBA region in the data -> averaged
#   3. else: nearest ancestor of the SBA region present in the data
#   4. else: leave the polygon gray (no coverage)
#
# Outputs HTML (interactive SVG grid) and PNG (composite raster).
#
# ---------------------------------------------------------------------
# Required packages:
#   install.packages(c("jsonlite", "xml2", "viridisLite",
#                      "scales", "stringr", "rsvg", "magick"))
# System libraries (for rsvg / magick):
#   librsvg2-dev  and  libmagick++-dev  (Debian/Ubuntu) or
#   librsvg, imagemagick (Homebrew on macOS)
# ---------------------------------------------------------------------


# -- Scalable Brain Atlas inputs ---------------------------------------
# The atlas panels need Scalable Brain Atlas coronal data for the Allen Mouse
# Brain Common Coordinate Framework v3 (SBA template ABA_v3), the same atlas
# version the cell counts are registered to. Expected in a "sba_data" folder
# containing rgb2acr.json, acr2parent.json, acr2full.json, and coronal_svg/.
# This is third-party data and is not redistributed with this repository; see
# the README for how to obtain it.
# SBA_DIR is set in the hardcoded paths block at the top of Config.
SVG_DIR <- file.path(SBA_DIR, "coronal_svg")

# Representative coronal slices evenly spread anterior -> posterior.
# (Olfactory-bulb slices 60 & 100 dropped, plus rightmost section 500.)
SELECTED_SLICES <- c(140, 180, 220, 260, 300, 340, 380, 420, 460)

# Fraction of a slice's width that each slice overlaps its left neighbour in
# the horizontal strip (0 = no overlap, 0.12 = tuck ~12% under the previous).
OVERLAP_FRAC <- 0.38

# Color for region (and brain) outlines. "#000000" gives the black anatomical
# borders seen in published atlas figures. Set to NA to make outlines match
# each region's fill (i.e. no visible borders).
OUTLINE_COLOR <- "#000000"

# Rasterization scale relative to each slice's native SVG width. Slices are
# supersampled at this scale, then the finished strip is downscaled to the
# print size below -- so keep this comfortably above the final resolution.
RENDER_SCALE <- 3

# -- Print sizing (for assembly into the NPP figure) -------------------
# Each strip is exported to be PLACED at this physical width and resolution.
# Layout: a 2-column (178 mm) figure with a PC column beside a WMDz column;
# ~3.3 in per strip fits two side by side with room for row labels + a gutter.
# Raise OUTPUT_DPI to 900-1200 if the strips will be enlarged beyond this width.
STRIP_WIDTH_IN  <- 3.35   # content ~82.4 mm = the A/B panel slot (side-by-side)
OUTPUT_DPI      <- 1200
STRIP_MARGIN_IN <- 0.05   # white margin around each strip (room for leaders)

# Horizontal colorbars sit under each A/B block, rendered at their FINAL
# placement size with pointsize = 8 so text is absolute pt (8 pt title, 8 pt
# ticks; NPP >= 8 pt floor) and the bar drops in 1:1 -- no downscaling.
COLORBAR_WIDTH_IN  <- 2.52   # ~64 mm bar under each A/B block
COLORBAR_HEIGHT_IN <- 0.50   # 8 pt title (top) + bar + 8 pt tick labels (bottom); raised from 0.42 to keep the bar thickness once text -> 8 pt

METRICS <- list(
  # PC: bounded 0->1, unsigned -> sequential map (no canonical colour in the
  # lineage, where PC is encoded as node size). Swap "viridis" for "mako",
  # "magma", "rocket", etc. to taste.
  PC   = list(label = "Participation Coefficient",    cmap = "viridis", diverging = FALSE),
  # WMDz: signed z-score -> canonical blue(low)-white(0)-red(high) diverging,
  # centered at 0 (Kimbrough 2020/2021; Ardinger 2024).
  WMDz = list(label = "Within-Module Degree Z-Score", cmap = "bwr",     diverging = TRUE)
)

# -- Data source config ------------------------------------------------
# The PC / WMDz values come from the Node_Roles sheet of the network
# pipeline's output workbook (03_Modules_and_Roles.xlsx).

# Conditions to plot, in order (must match the Treatment column values).

# Which network(s) to render. One figure set is produced per sex here.
# Options present in the sheet: "Combined", "Males", "Females".
SEXES_TO_PLOT <- c("Combined", "Males", "Females")

# Color-scale groups: sexes within the same group share ONE scale per metric
# (so panels are directly comparable), and get one shared legend. Combined sits
# alone (main figure); Males+Females share a scale (supplemental sex comparison).
# Each group name is used in the legend filename. Sexes not in SEXES_TO_PLOT are
# ignored automatically.
SCALE_GROUPS <- list(
  Combined = c("Combined"),
  BySex    = c("Males", "Females")
)

# Which Node_Roles columns feed each metric.
#   Canonical / binary (matches the original labels and the Kimbrough lineage):
#       PC   -> "Participation_coef"
#       WMDz -> "WM_degree_z"
#   Weighted variants (the ones that drive the cartographic G-A roles):
#       PC   -> "Participation_coef_wt"
#       WMDz -> "WM_strength_z"
PC_COLUMN   <- "Participation_coef_wt"   # weighted PC -- matches the roles, the
WMDZ_COLUMN <- "WM_strength_z"           # network graphs, and Table S5

# -- Fixed colours -----------------------------------------------------
NO_DATA_COLOR    <- "#bfbfbf"  # regions with no metric value
BACKGROUND_COLOR <- "#ffffff"


# -- Colormap LUT (matplotlib-equivalent) ------------------------------
# matplotlib's get_cmap(name)(t) uses a 256-entry lookup table; we mirror
# that so colors match closely. t is expected in [0, 1] (vectorized).
.cmap_lut <- function(cmap_name, n = 256) {
  switch(cmap_name,
         viridis = viridisLite::viridis(n),
         magma   = viridisLite::magma(n),
         inferno = viridisLite::inferno(n),
         plasma  = viridisLite::plasma(n),
         cividis = viridisLite::cividis(n),
         mako    = viridisLite::mako(n),
         rocket  = viridisLite::rocket(n),
         # blue -> white -> red diverging (ColorBrewer RdBu, reversed so blue=low,
         # red=high). Use with a 0-centered (diverging) norm for signed metrics.
         bwr     = grDevices::colorRampPalette(c(
           "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7",
           "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(n),
         viridisLite::viridis(n)  # default fallback
  )
}

cmap_color <- function(cmap_name, t) {
  lut <- .cmap_lut(cmap_name)
  n <- length(lut)
  t <- pmin(pmax(t, 0), 1)            # clamp to [0,1]
  idx <- round(t * (n - 1)) + 1L
  substr(lut[idx], 1, 7)             # strip alpha if present -> #RRGGBB
}


# -- SBA atlas loading -------------------------------------------------

load_sba <- function() {
  # Read file contents explicitly, then parse. (Passing a bare path to
  # fromJSON can make it try to parse the path string itself as JSON.)
  read_json_file <- function(path, required = TRUE) {
    if (!file.exists(path)) {
      if (required) stop("Cannot find atlas file: ", path)
      warning("Optional atlas file not found; continuing without it: ", path)
      return(list())
    }
    jsonlite::fromJSON(paste(readLines(path, warn = FALSE), collapse = "\n"),
                       simplifyVector = FALSE)
  }
  
  rgb2acr <- read_json_file(file.path(SBA_DIR, "rgb2acr.json"), required = TRUE)
  names(rgb2acr) <- toupper(names(rgb2acr))
  acr2parent <- read_json_file(file.path(SBA_DIR, "acr2parent.json"), required = FALSE)
  acr2full   <- read_json_file(file.path(SBA_DIR, "acr2full.json"),   required = FALSE)
  
  # Build children index: parent -> character vector of child acronyms
  children <- list()
  for (acr in names(acr2parent)) {
    parent <- acr2parent[[acr]]
    if (is.null(parent) || is.na(parent) || !nzchar(parent)) next
    children[[parent]] <- c(children[[parent]], acr)
  }
  
  list(rgb2acr = rgb2acr, acr2parent = acr2parent,
       acr2full = acr2full, children = children)
}


all_descendants <- function(acr, children, include_self = TRUE) {
  out <- if (include_self) acr else character(0)
  stack <- children[[acr]]
  if (is.null(stack)) stack <- character(0)
  while (length(stack) > 0L) {
    a <- stack[length(stack)]
    stack <- stack[-length(stack)]
    if (a %in% out) next
    out <- c(out, a)
    kids <- children[[a]]
    if (!is.null(kids)) stack <- c(stack, kids)
  }
  out
}


ancestors <- function(acr, acr2parent) {
  out <- character(0)
  cur <- acr2parent[[acr]]
  while (!is.null(cur) && !is.na(cur) && nzchar(cur)) {
    out <- c(out, cur)
    cur <- acr2parent[[cur]]
  }
  out
}


# -- Metric value resolution -------------------------------------------

# my_values: named list/vector  acronym -> numeric (the data we have).
# Returns a closure: sba_acronym -> numeric (NA if no coverage).
# Resolution order: exact, descendants-mean, nearest-ancestor.
build_resolver <- function(my_values, acr2parent, children) {
  is_good <- vapply(my_values, function(v) !is.null(v) && !is.na(v), logical(1))
  data_keys <- names(my_values)[is_good]
  
  function(sba_acr) {
    if (sba_acr %in% data_keys) {
      return(as.numeric(my_values[[sba_acr]]))
    }
    desc <- all_descendants(sba_acr, children, include_self = FALSE)
    hit_keys <- intersect(desc, data_keys)
    if (length(hit_keys) > 0L) {
      return(mean(as.numeric(unlist(my_values[hit_keys]))))
    }
    for (anc in ancestors(sba_acr, acr2parent)) {
      if (anc %in% data_keys) {
        return(as.numeric(my_values[[anc]]))
      }
    }
    NA_real_
  }
}


# -- Color mapping -----------------------------------------------------
# In Python this returned a matplotlib Normalize object. Here we return
# a small list(vmin, vmax) plus a normalize() helper.

make_norm <- function(values, diverging) {
  arr <- as.numeric(values)
  arr <- arr[!is.na(arr)]
  if (length(arr) == 0L) return(list(vmin = 0, vmax = 1, ticks = c(0, 1)))
  if (diverging) {
    m <- max(abs(min(arr)), abs(max(arr)))
    p <- pretty(c(-m, m)); step <- if (length(p) > 1) p[2] - p[1] else m
    vmax <- ceiling(m / step) * step               # round OUT so bar ends on a tick
    return(list(vmin = -vmax, vmax = vmax, ticks = seq(-vmax, vmax, by = step)))
  }
  lo <- min(arr); hi <- max(arr)
  p <- pretty(c(lo, hi)); step <- if (length(p) > 1) p[2] - p[1] else (hi - lo)
  vmin <- floor(lo / step) * step                  # round the ends onto tick marks so the
  vmax <- ceiling(hi / step) * step                # coloured bar starts/stops at a label
  list(vmin = vmin, vmax = vmax, ticks = seq(vmin, vmax, by = step))
}

normalize_value <- function(v, norm) {
  if (norm$vmax == norm$vmin) return(rep(0.5, length(v)))
  (v - norm$vmin) / (norm$vmax - norm$vmin)
}

hex_color <- function(value, norm, cmap_name) {
  if (is.null(value) || is.na(value)) return(NO_DATA_COLOR)
  cmap_color(cmap_name, normalize_value(value, norm))
}


# -- SVG recoloring ----------------------------------------------------
# The Python version used a regex tuned to the SBA's exported format
# (one <path> per region, identical fill & stroke hex). Here we parse the
# SVG with xml2 instead: more robust to attribute ordering/whitespace,
# and we preserve the xmlns by querying with local-name() rather than
# stripping namespaces. We recolor only paths that already carry both a
# hex fill and a stroke, matching the original regex's intent.

recolor_svg <- function(svg_text, rgb2acr, resolver, norm, cmap_name) {
  doc <- read_xml(svg_text)
  
  hex6 <- "^#?[0-9A-Fa-f]{6}$"
  paths <- xml_find_all(doc, "//*[local-name()='path']")
  for (p in paths) {
    fill <- xml_attr(p, "fill")
    if (is.na(fill) || !grepl(hex6, fill)) next
    # rgb2acr keys are bare 6-char hex (no '#'); strip it before lookup.
    acr <- rgb2acr[[toupper(sub("^#", "", fill))]]
    if (is.null(acr)) next
    val <- resolver(acr)
    new_hex <- hex_color(val, norm, cmap_name)   # no data -> NO_DATA_COLOR (grey)
    xml_set_attr(p, "fill", new_hex)
    xml_set_attr(p, "stroke", if (is.na(OUTLINE_COLOR)) new_hex else OUTLINE_COLOR)
  }
  
  # Make the canvas background transparent (fill="none") so slices can be
  # overlapped without opaque rectangles occluding their neighbours.
  rects <- xml_find_all(doc, "//*[local-name()='rect']")
  for (r in rects) {
    rf <- xml_attr(r, "fill")
    if (!is.na(rf) && tolower(rf) == "#000000") {
      xml_set_attr(r, "fill", "none")
    }
  }
  
  as.character(doc)
}


# -- Rasterization -----------------------------------------------------
# cairosvg.svg2png(scale=) has no direct rsvg equivalent, so we read the
# SVG's intrinsic size and scale the requested width.

slice_to_image <- function(svg_text, scale = 1.0) {
  bytes <- charToRaw(svg_text)
  nw <- svg_native_width(svg_text)
  png_raw <- if (!is.na(nw) && nw > 0) {
    rsvg::rsvg_png(bytes, width = max(1L, round(nw * scale)))
  } else {
    rsvg::rsvg_png(bytes)  # fall back to the SVG's native size
  }
  magick::image_read(png_raw)
}

# Read the SVG's intrinsic pixel width from its width= attribute or, failing
# that, the third value of its viewBox. Replaces rsvg::rsvg_dimensions(),
# which isn't exported in all rsvg versions.
svg_native_width <- function(svg_text) {
  doc <- xml2::read_xml(svg_text)
  w <- xml2::xml_attr(doc, "width")
  if (!is.na(w)) {
    wn <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", w)))
    if (!is.na(wn) && wn > 0) return(wn)
  }
  vb <- xml2::xml_attr(doc, "viewBox")
  if (!is.na(vb)) {
    parts <- suppressWarnings(as.numeric(strsplit(trimws(vb), "[ ,]+")[[1]]))
    if (length(parts) == 4 && !is.na(parts[3]) && parts[3] > 0) return(parts[3])
  }
  NA_real_
}


# -- Composition: PNG --------------------------------------------------
# Replaces matplotlib subplot composition with magick montage.

colorbar_image <- function(norm, cmap_name, label) {
  tf <- tempfile(fileext = ".png")
  # Render at the FINAL placement size with pointsize = 8, so text is absolute:
  # cex 1 -> 8 pt (bold title) and cex.axis 1 -> 8 pt (ticks). NPP >= 8 pt floor;
  # placed 1:1 at assembly.
  grDevices::png(tf,
                 width  = round(COLORBAR_WIDTH_IN  * OUTPUT_DPI),
                 height = round(COLORBAR_HEIGHT_IN * OUTPUT_DPI),
                 res = OUTPUT_DPI, pointsize = 8, type = dev_type)
  # Safety net: if an error occurs mid-plot, still close the device.
  ok <- FALSE
  on.exit(if (!ok && grDevices::dev.cur() > 1L) grDevices::dev.off(), add = TRUE)
  
  # Horizontal bar: bold 8 pt metric title on TOP (side 3), 8 pt tick labels on
  # the BOTTOM (side 1). Sits under each A/B block in the assembly.
  op <- graphics::par(mar = c(1.6, 1.0, 1.4, 1.0))   # L/R 0.8 -> 1.0: 8 pt end labels (0.0 / 0.8, -3 / 3) need more room
  n <- 256L
  cols <- cmap_color(cmap_name, (0:(n - 1)) / (n - 1))
  graphics::image(
    x = seq(norm$vmin, norm$vmax, length.out = n), y = 1,
    z = matrix(seq_len(n), nrow = n, ncol = 1),
    col = cols, axes = FALSE, xlab = "", ylab = ""
  )
  graphics::axis(1, at = norm$ticks, cex.axis = 1, mgp = c(3, 0.4, 0), tcl = -0.3)
  graphics::mtext(label, side = 3, line = 0.3, font = 2, cex = 1)
  graphics::par(op)
  
  grDevices::dev.off()   # close & flush to disk BEFORE reading it back
  ok <- TRUE
  magick::image_read(tf)
}

compose_png <- function(slice_imgs, png_path,
                        overlap_frac = OVERLAP_FRAC, bg = "white") {
  # Trim each slice's transparent margin so neighbours tuck together by their
  # actual brain content rather than their bounding boxes.
  trimmed <- lapply(slice_imgs, magick::image_trim)
  info <- lapply(trimmed, magick::image_info)
  ws <- vapply(info, function(x) x$width,  numeric(1))
  hs <- vapply(info, function(x) x$height, numeric(1))
  H  <- max(hs)
  
  # Horizontal overlap, sized from the mean slice width.
  ov <- round(overlap_frac * mean(ws))
  
  # Left edge x for each slice (each starts ov px before its predecessor ends).
  xs <- numeric(length(trimmed))
  cur <- 0
  for (i in seq_along(trimmed)) {
    xs[i] <- cur
    cur <- cur + ws[i] - ov
  }
  total_w <- cur + ov
  
  # Composite onto a solid background canvas. Slices keep transparent margins,
  # so overlaps reveal the neighbour beneath rather than an opaque box; later
  # (posterior) slices land on top.
  strip <- magick::image_blank(total_w, H, color = bg)
  for (i in seq_along(trimmed)) {
    yoff <- round((H - hs[i]) / 2)
    strip <- magick::image_composite(
      strip, trimmed[[i]], offset = sprintf("+%d+%d", round(xs[i]), yoff)
    )
  }
  
  # Downscale the supersampled strip to the target print width, add the margin,
  # and tag the PNG resolution so it imports at STRIP_WIDTH_IN @ OUTPUT_DPI.
  margin_px <- round(STRIP_MARGIN_IN * OUTPUT_DPI)
  inner_px  <- max(1L, round(STRIP_WIDTH_IN * OUTPUT_DPI) - 2L * margin_px)
  strip <- magick::image_resize(strip, as.character(inner_px))
  if (margin_px > 0) {
    strip <- magick::image_border(strip, bg, sprintf("%dx%d", margin_px, margin_px))
  }
  magick::image_write(strip, path = png_path, format = "png",
                      density = as.character(OUTPUT_DPI))
}

# Write a standalone horizontal colorbar (legend) for one metric's shared scale,
# placed 1:1 under its A/B strip block at assembly.
write_legend <- function(norm, cmap_name, label, png_path) {
  cb <- colorbar_image(norm, cmap_name, label)
  magick::image_write(cb, path = png_path, format = "png",
                      density = as.character(OUTPUT_DPI))
}


# -- Data loader (reads the Node_Roles sheet) --------------------------
# Builds the nested structure the main loop expects:
#   study(= Sex) -> group(= Treatment) -> list(
#       participation_coefficient = named list  region -> PC,
#       wmdz                      = named list  module_id -> (region -> WMDz)
#   )
# WMDz is kept nested by Network_module to mirror the original Python
# contract; flatten_wmdz() collapses it to region -> value downstream.

load_all_studies <- function() {
  if (!file.exists(SHEET_FILE)) {
    stop("Cannot find data workbook: ", SHEET_FILE)
  }
  roles <- as.data.frame(read_excel(SHEET_FILE, sheet = SHEET_NAME),
                         stringsAsFactors = FALSE)
  
  needed <- c("Treatment", "Sex", "Region", "Network_module",
              PC_COLUMN, WMDZ_COLUMN)
  missing <- setdiff(needed, names(roles))
  if (length(missing) > 0) {
    stop("Sheet '", SHEET_NAME, "' is missing column(s): ",
         paste(missing, collapse = ", "))
  }
  
  roles$Region <- trimws(as.character(roles$Region))
  roles <- roles[!is.na(roles$Region) & nzchar(roles$Region), ]
  
  studies <- list()
  
  for (sx in SEXES_TO_PLOT) {
    sub_sex <- roles[roles$Sex == sx, , drop = FALSE]
    if (nrow(sub_sex) == 0L) {
      message("  Warning: no rows for Sex == '", sx, "'")
      next
    }
    groups <- list()
    for (trt in TREATMENTS) {
      sub <- sub_sex[sub_sex$Treatment == trt, , drop = FALSE]
      if (nrow(sub) == 0L) {
        message("  Warning: no rows for Treatment == '", trt, "' (Sex ", sx, ")")
        next
      }
      
      # participation_coefficient: region -> value
      pc <- as.list(suppressWarnings(as.numeric(sub[[PC_COLUMN]])))
      names(pc) <- sub$Region
      
      # wmdz: module_id -> (region -> value)
      wmdz <- list()
      wm_vals <- suppressWarnings(as.numeric(sub[[WMDZ_COLUMN]]))
      for (i in seq_len(nrow(sub))) {
        mod <- as.character(sub$Network_module[i])
        if (is.null(wmdz[[mod]])) wmdz[[mod]] <- list()
        wmdz[[mod]][[sub$Region[i]]] <- wm_vals[i]
      }
      
      groups[[trt]] <- list(participation_coefficient = pc, wmdz = wmdz)
    }
    if (length(groups) > 0L) studies[[sx]] <- groups
  }
  
  if (length(studies) == 0L) {
    stop("No data assembled -- check SEXES_TO_PLOT / TREATMENTS against the sheet.")
  }
  studies
}

safe_filename <- function(s) {
  s <- gsub("[^A-Za-z0-9._-]+", "_", s)
  gsub("^_+|_+$", "", s)
}


# -- Main --------------------------------------------------------------

# Each region belongs to one cluster; merge to {region: value}.
flatten_wmdz <- function(wmdz_by_cluster) {
  out <- list()
  for (cid in names(wmdz_by_cluster)) {
    region_map <- wmdz_by_cluster[[cid]]
    for (region in names(region_map)) {
      out[[region]] <- region_map[[region]]
    }
  }
  out
}

main <- function() {
  sba <- load_sba()
  rgb2acr    <- sba$rgb2acr
  acr2parent <- sba$acr2parent
  acr2full   <- sba$acr2full
  children   <- sba$children
  
  # Cache SVGs for selected slices once.
  slice_svgs <- list()
  for (idx in SELECTED_SLICES) {
    path <- file.path(SVG_DIR, sprintf("Annotation2014_141_%04d.svg", idx))
    slice_svgs[[as.character(idx)]] <-
      paste(readLines(path, warn = FALSE), collapse = "\n")
  }
  
  all_data <- load_all_studies()
  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
  
  # SBA acronyms that actually appear in the selected slices. This depends
  # only on the SVGs, so compute it once and reuse for every group/metric.
  acrs_in_view <- character(0)
  for (svg in slice_svgs) {
    fills <- str_match_all(svg, 'fill="#([0-9A-Fa-f]{6})"')[[1]]
    if (nrow(fills) > 0L) {
      for (f in fills[, 2]) {
        acr <- rgb2acr[[toupper(f)]]
        if (!is.null(acr)) acrs_in_view <- c(acrs_in_view, acr)
      }
    }
  }
  acrs_in_view <- unique(acrs_in_view)
  
  # Determine the single-vs-multi-sex suffix policy: when only one sex is being
  # plotted, filenames stay clean (no sex tag); otherwise every strip is tagged.
  sexes_present <- names(all_data)
  tag_sex <- length(sexes_present) > 1L
  
  for (grp_label in names(SCALE_GROUPS)) {
    sexes <- intersect(SCALE_GROUPS[[grp_label]], sexes_present)
    if (length(sexes) == 0L) next
    
    for (metric_key in names(METRICS)) {
      meta <- METRICS[[metric_key]]
      
      # Build a resolver for every (sex x condition) in this scale group, and
      # pool ALL their values so the colour scale is shared across the group.
      resolvers   <- list()        # keyed "sex||condition"
      pooled_vals <- numeric(0)
      for (sx in sexes) {
        for (group_name in names(all_data[[sx]])) {
          gd <- all_data[[sx]][[group_name]]
          my_values <- if (metric_key == "PC") gd$participation_coefficient
          else flatten_wmdz(gd$wmdz)
          resolver <- build_resolver(my_values, acr2parent, children)
          resolvers[[paste(sx, group_name, sep = "||")]] <- resolver
          for (acr in acrs_in_view) {
            v <- resolver(acr)
            if (!is.na(v)) pooled_vals <- c(pooled_vals, v)
          }
        }
      }
      
      # ONE norm for the whole scale group, applied to every panel in it.
      norm <- make_norm(pooled_vals, meta$diverging)
      cat(sprintf("[%s / %s] shared scale across {%s}: [%.3f, %.3f]\n",
                  grp_label, metric_key, paste(sexes, collapse = ", "),
                  norm$vmin, norm$vmax))
      
      # One shared legend per scale group + metric.
      legend_base <- sprintf("regional_%s_%s_legend", metric_key, grp_label)
      write_legend(norm, meta$cmap, meta$label,
                   file.path(OUTPUT_DIR, paste0(legend_base, ".png")))
      cat(sprintf("  wrote %s.png  (shared legend)\n", legend_base))
      
      # Render each (sex, condition) strip with the group's shared scale.
      for (sx in sexes) {
        for (group_name in names(all_data[[sx]])) {
          resolver <- resolvers[[paste(sx, group_name, sep = "||")]]
          
          slice_imgs <- list()
          for (idx in SELECTED_SLICES) {
            svg_new <- recolor_svg(
              slice_svgs[[as.character(idx)]], rgb2acr, resolver, norm, meta$cmap
            )
            slice_imgs[[length(slice_imgs) + 1L]] <- slice_to_image(svg_new, scale = RENDER_SCALE)
          }
          n_colored <- sum(vapply(acrs_in_view,
                                  function(a) !is.na(resolver(a)), logical(1)))
          
          base <- if (tag_sex)
            sprintf("regional_%s_%s_%s", metric_key, safe_filename(group_name),
                    safe_filename(sx))
          else
            sprintf("regional_%s_%s", metric_key, safe_filename(group_name))
          compose_png(slice_imgs, file.path(OUTPUT_DIR, paste0(base, ".png")))
          
          cat(sprintf("  wrote %s.png  (n=%d colored regions)\n",
                      base, n_colored))
        }
      }
    }
  }
  
  cat("Done.\n")
}

# #############################################################################
# ## SECTION C -- Supplemental Figure S4: cartographic scatter (PC vs WMDz)   ##
# #############################################################################
# =============================================================================
# Supplemental Figure S4 (the joint PC-WMDz cartographic-role scatter; formerly
# numbered S3 before the NPP figure reorganization).
#
# Guimera-Amaral cartographic scatter: within-module degree z-score (WMDz, y)
# vs participation coefficient (PC, x), one panel per condition, Combined
# (pooled-sex) network. Node colour = CCFv3 division (canonical divisions.R).
# Connector hubs (Node_role == "R6 connector hub") are labelled in an open band
# above the cloud, tied to their node by a leader line so nothing overlaps.
#
# Reads the same Node_Roles sheet as Panels A/B (SHEET_FILE / PC_COLUMN /
# WMDZ_COLUMN above). Independent of the SBA atlas -- needs only readxl + base
# graphics -- so it renders even when the atlas inputs are absent.
#
# NPP artwork spec: rendered at final 178 mm width, placed 1:1, pointsize = 8 so
# every text element is 8 pt (>= NPP floor); panel/axis titles carry the
# hierarchy through BOLD, not larger type. 1200 dpi TIFF (LZW) + PNG.
# =============================================================================

# -- Section C config ---------------------------------------------------------
if (!exists("FONT")) FONT <- "Arial"

SC_TREATMENTS <- c("Vehicle", "Acute Morphine", "Morphine-Dependent", "Ro 64-6198")
# Sexes to render (one 2x2 figure each). "Combined" is the S4 figure; add
# "Males","Females" for sex-stratified versions (same code path).
SCATTER_SEXES <- c("Combined")

SC_XLIM      <- c(-0.02, 0.82)          # PC axis
SC_YLIM      <- c(-2.85, 3.90)          # WMDz axis (extra headroom for the label band)
SC_PC_GUIDES <- c(0.05, 0.30, 0.625, 0.75)   # Guimera-Amaral P cut points (dashed)
SC_WMDZ_HUB  <- 1.5                     # hub / non-hub boundary (solid line)
SC_XTICKS    <- c(0.0, 0.2, 0.4, 0.6, 0.8)
SC_YTICKS    <- c(-2, -1, 0, 1, 2, 3)

# Label-band placement: keep sorted order, enforce a minimum centre-to-centre
# gap, stack onto two rows when a panel has > 4 connector hubs.
SC_LX_LO <- 0.05; SC_LX_HI <- 0.78; SC_MIN_GAP <- 0.135
SC_ROW_Y1 <- 3.25                       # single-row band height (<= 4 labels)
SC_ROW_Y2 <- c(3.50, 2.72)              # two-row band heights (> 4 labels)

SC_W_IN <- 178 / 25.4                   # 178 mm double-column width (place 1:1)
SC_H_IN <- 150 / 25.4                   # ~150 mm tall (2x2 + side legend)

# -- Canonical division palette ----------------------------------------------
# Reuse divisions.R (single source of truth) if reachable; otherwise fall back
# to an inline copy of division_order / division_colors so this section renders
# standalone. Either way the colours match Figs 2-5 (Striatum = #ff5da2, etc.).
if (!exists("division_colors") || !exists("division_order")) {
  .div_path <- file.path(if (exists("SCRIPT_DIR")) SCRIPT_DIR else "scripts",
                         "divisions.R")
  if (file.exists(.div_path)) try(source(.div_path), silent = TRUE)
}
if (!exists("division_colors") || !exists("division_order")) {
  division_order <- c(
    "Isocortex", "Olfactory", "Hippocampus", "Cortical Subplate",
    "Amygdala", "Striatum", "Pallidum/Septum", "Thalamus",
    "Hypothalamus", "Midbrain", "Pons", "Medulla", "Cerebellum")
  division_colors <- c(
    "Isocortex" = "#3949ab", "Olfactory" = "#4dad4c",
    "Hippocampus" = "#ff9800", "Cortical Subplate" = "#00838f",
    "Amygdala" = "#9c27b0", "Striatum" = "#ff5da2",
    "Pallidum/Septum" = "#f64234", "Thalamus" = "#607d8b",
    "Hypothalamus" = "#4cd0e9", "Midbrain" = "#ec1c67",
    "Pons" = "#7b5548", "Medulla" = "#ba7517", "Cerebellum" = "#fdd835")
}

# -- Helpers ------------------------------------------------------------------
# White-halo text (8-way offset white copies, then the label on top) so region
# codes stay legible over dots and guide lines.
sc_halo_text <- function(x, y, labels, cex = 1, halo = "white", col = "black") {
  d <- strwidth("o", cex = cex) * 0.18
  for (a in seq(0, 2 * pi, length.out = 9)[-1])
    text(x + cos(a) * d, y + sin(a) * d, labels, cex = cex, col = halo,
         adj = c(0.5, 0.5))
  text(x, y, labels, cex = cex, col = col, adj = c(0.5, 0.5))
}

# 1-D label de-clutter on a sorted vector: push neighbours apart to SC_MIN_GAP,
# then keep the whole row inside [SC_LX_LO, SC_LX_HI] by rolling back from
# whichever end overflows.
sc_declutter <- function(v) {
  n <- length(v)
  if (n > 1) for (i in 2:n) if (v[i] < v[i - 1] + SC_MIN_GAP) v[i] <- v[i - 1] + SC_MIN_GAP
  if (n && v[n] > SC_LX_HI) {
    v[n] <- SC_LX_HI
    if (n > 1) for (i in (n - 1):1) v[i] <- min(v[i], v[i + 1] - SC_MIN_GAP)
  }
  if (n && v[1] < SC_LX_LO) {
    v[1] <- SC_LX_LO
    if (n > 1) for (i in 2:n) v[i] <- max(v[i], v[i - 1] + SC_MIN_GAP)
  }
  v
}

# One cartographic scatter panel (one condition).
sc_draw_panel <- function(sub, panel_title) {
  plot(NA, xlim = SC_XLIM, ylim = SC_YLIM, xaxs = "i", yaxs = "i",
       axes = FALSE, xlab = "", ylab = "")
  for (gx in SC_PC_GUIDES) abline(v = gx, lty = 2, lwd = 0.6, col = "#9a9a9a")
  abline(h = SC_WMDZ_HUB, lwd = 1.0, col = "black")
  # nodes, drawn in canonical division order (fixes legend/z-order)
  for (dv in division_order) {
    m <- sub$Structure == dv & !is.na(sub$PCv) & !is.na(sub$WMv)
    if (any(m))
      points(sub$PCv[m], sub$WMv[m], pch = 21, cex = 0.8,
             bg = division_colors[[dv]], col = "white", lwd = 0.3)
  }
  # L-shaped axes only (no top/right spine)
  axis(1, at = SC_XTICKS, cex.axis = 1, lwd = 0.6, mgp = c(3, 0.35, 0), tcl = -0.3)
  axis(2, at = SC_YTICKS, cex.axis = 1, las = 1, lwd = 0.6, mgp = c(3, 0.55, 0), tcl = -0.3)
  title(main = panel_title, font.main = 2, cex.main = 1, line = 0.5)
  # connector-hub (R6) labels, parked in the band above the cloud
  hubs <- sub[!is.na(sub$Role) & sub$Role == "R6 connector hub", , drop = FALSE]
  hubs <- hubs[order(hubs$PCv), , drop = FALSE]
  n <- nrow(hubs)
  if (n > 0) {
    n_rows <- if (n <= 4) 1L else 2L
    row_y  <- if (n_rows == 1L) SC_ROW_Y1 else SC_ROW_Y2
    for (r in seq_len(n_rows)) {
      idx <- which(((seq_len(n) - 1L) %% n_rows) == (r - 1L))
      if (!length(idx)) next
      g  <- hubs[idx, , drop = FALSE]
      lx <- sc_declutter(pmin(pmax(g$PCv, SC_LX_LO), SC_LX_HI))
      ly <- row_y[r]
      for (k in seq_along(idx)) {
        segments(lx[k], ly, g$PCv[k], g$WMv[k], col = "grey40", lwd = 0.6)
        points(g$PCv[k], g$WMv[k], pch = 1, cex = 1.05, col = "grey25", lwd = 0.6)
        sc_halo_text(lx[k], ly, g$Region[k], cex = 1)
      }
    }
  }
}

# Full 2x2 figure + side legend + shared axis titles (deterministic NDC layout).
sc_draw_figure <- function(dat) {
  PL <- 0.085; PR <- 0.80; PB <- 0.10; PT <- 0.94   # panel-block extent (NDC)
  gx <- 0.045; gy <- 0.075                          # inter-panel gaps
  cw <- (PR - PL - gx) / 2; rh <- (PT - PB - gy) / 2
  pf <- list(c(PL, PL + cw, PB + rh + gy, PT),       # top-left
             c(PR - cw, PR, PB + rh + gy, PT),       # top-right
             c(PL, PL + cw, PB, PB + rh),            # bottom-left
             c(PR - cw, PR, PB, PB + rh))            # bottom-right
  par(family = FONT)
  for (i in seq_along(SC_TREATMENTS)) {
    par(fig = pf[[i]], new = (i > 1), mar = c(2.2, 2.7, 1.8, 0.5))
    sc_draw_panel(dat[dat$Treatment == SC_TREATMENTS[i], , drop = FALSE],
                  SC_TREATMENTS[i])
  }
  # legend (canonical order, divisions present in the data)
  par(fig = c(0.805, 1.0, PB, PT), new = TRUE, mar = c(0, 0, 0, 0))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
       axes = FALSE, xlab = "", ylab = "")
  present <- division_order[division_order %in% unique(dat$Structure)]
  text(0.02, 0.985, "Structure", adj = c(0, 1), font = 2, cex = 1)
  y0 <- 0.90; dy <- 0.058
  for (j in seq_along(present)) {
    yy <- y0 - (j - 1) * dy
    points(0.07, yy, pch = 21, bg = division_colors[[present[j]]],
           col = "white", cex = 1.3, lwd = 0.5)
    text(0.14, yy, present[j], adj = c(0, 0.5), cex = 1)
  }
  # shared axis titles (single centred labels over the panel block only)
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
       axes = FALSE, xlab = "", ylab = "")
  text((PL + PR) / 2, 0.028, "Participation Coefficient (PC)",
       adj = c(0.5, 0.5), font = 2, cex = 1)
  text(0.022, (PB + PT) / 2, "Within-Module Degree Z-score (WMDz)",
       srt = 90, adj = c(0.5, 0.5), font = 2, cex = 1)
}

# -- Runner -------------------------------------------------------------------
run_cartography_scatter <- function() {
  roles <- as.data.frame(read_excel(SHEET_FILE, sheet = SHEET_NAME),
                         stringsAsFactors = FALSE)
  roles$Region    <- trimws(as.character(roles$Region))
  roles$Structure <- as.character(roles$Structure)
  roles$Role      <- as.character(roles$Node_role)
  roles$PCv       <- suppressWarnings(as.numeric(roles[[PC_COLUMN]]))
  roles$WMv       <- suppressWarnings(as.numeric(roles[[WMDZ_COLUMN]]))
  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
  
  for (sx in SCATTER_SEXES) {
    d <- roles[roles$Sex == sx, , drop = FALSE]
    if (nrow(d) == 0L) { message("  [S4] no rows for Sex == '", sx, "'"); next }
    stub <- file.path(OUTPUT_DIR, sprintf("FigureS4_scatter_%s", safe_filename(sx)))
    tiff(paste0(stub, ".tiff"), width = SC_W_IN, height = SC_H_IN, units = "in",
         res = OUTPUT_DPI, compression = "lzw", family = FONT, pointsize = 8,
         type = dev_type)
    sc_draw_figure(d); dev.off()
    png(paste0(stub, ".png"), width = SC_W_IN, height = SC_H_IN, units = "in",
        res = OUTPUT_DPI, family = FONT, pointsize = 8, type = dev_type)
    sc_draw_figure(d); dev.off()
    cat(sprintf("  wrote %s.{tiff,png}  (Fig. S4 cartographic scatter, Sex = %s)\n",
                basename(stub), sx))
  }
}

# =============================================================================
# RUN
# =============================================================================
# Panels A/B (regional atlas maps) -- require the SBA atlas inputs (rsvg/magick
# + coronal SVGs). Wrapped in tryCatch so Supplemental Figure S4 (Section C)
# still renders when those inputs are absent: S4 needs only the Node_Roles
# sheet + base graphics.
tryCatch(main(), error = function(e)
  message("[Panels A/B] atlas render skipped: ", conditionMessage(e)))

# Supplemental Figure S4 -- cartographic scatter (PC vs WMDz).
run_cartography_scatter()