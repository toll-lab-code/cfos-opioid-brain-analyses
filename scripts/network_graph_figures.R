# network_graph_figures.R
# =============================================================================
# Figure 4, Panel C: force-directed co-activation network graphs, plus the
# cartographic role-composition bars.
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
#   Rscript scripts/network_graph_figures.R
# Input and output locations are set in the Config section below; see README.
# =============================================================================
#
# Run network_analysis.R first: this script reads 01_Connectivity.xlsx and
# 03_Modules_and_Roles.xlsx, which that step produces. Graph conventions follow
# Kimbrough 2020 Fig 6, Bloch 2024 Fig 4, and Ardinger 2024.
#
# Network graph encoding (per condition x sex):
#   node size  = participation coefficient (PC, weighted)   smaller = lower PC
#   node fill  = within-module degree z-score (WMDz)         blue<0<red, white=0
#   node ring  = consensus-Louvain module (within-panel)
#   edge width = co-activation strength r (positive edges, density 0.10)
#   layout     = force-directed (Kamada-Kawai; see layout_type below)
#
# Graph edges reproduce threshold_to_density() from network_analysis.R
# verbatim, so the drawn graph is the one the modules/WMDz/PC were computed on.
#
# Outputs (all 600 dpi TIFF + PNG, individual panels for Illustrator assembly):
#   network_<Sex>_<Abbr>.{tiff,png}        12 graphs (4 cond x 3 sex)
#   network_legend.{tiff,png}              shared key
#   roles_<Sex>.{tiff,png}                 3 role-composition bar panels
# =============================================================================

## ------------------------------ PACKAGES -----------------------------------
required_pkgs <- c("readxl", "igraph", "magick")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages. Please run:\n  install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}
library(readxl)
library(igraph)

## ------------------------------- CONFIG ------------------------------------
## Paths are relative to the repository root. Run from the repo root, e.g.:
##   Rscript scripts/network_graph_figures.R
## Override with the CFOS_* environment variables to run against a working
## folder laid out differently.
DATA_DIR    <- path.expand(Sys.getenv("CFOS_DATA_DIR",   "data"))     # region_division_lookup.csv, read by divisions.R
OUTPUT_ROOT <- path.expand(Sys.getenv("CFOS_OUTPUT_DIR", "results"))  # the two input workbooks, and figures out
SCRIPT_DIR  <- path.expand(Sys.getenv("CFOS_SCRIPT_DIR", "scripts"))  # divisions.R

conn_file  <- file.path(OUTPUT_ROOT, "01_Connectivity.xlsx")
roles_file <- file.path(OUTPUT_ROOT, "03_Modules_and_Roles.xlsx")
out_dir    <- file.path(OUTPUT_ROOT, "figures", "Figure4C")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Graphics device. "cairo" writes plain sRGB with NO embedded colour profile,
## which is what regional_figures.R and coactivation_figures.R already produce
## and what the PIL assembler expects. The macOS quartz device (the default for
## png()/tiff() there) embeds an Apple "Generic RGB" profile AND stores shifted
## raw pixel values -- e.g. #1B9E77 is written as #1D8F64 -- so a division
## colour composited from a quartz panel does not match the same colour taken
## from the other panels, because PIL pastes raw pixels and ignores the profile.
## Set CFOS_DEVICE=quartz only if cairo cannot resolve the font on your machine.
dev_type <- Sys.getenv("CFOS_DEVICE", "cairo")

font_fam <- Sys.getenv("CFOS_FONT", "Arial")
dev_res  <- 1200
panel_in <- 69 / 25.4         # square network panels rendered at FINAL size (69 mm); place 1:1
edge_density <- 0.10          # primary_density (matches pipeline). Named so it does
#                               not shadow stats::density.
layout_seed <- 7             # reproducible layout
giant_only <- TRUE           # keep only the largest connected component
# (drops detached islands + the pendant spokes read cleaner)
min_degree <- 2              # draw the k-core: iteratively drop nodes with degree
# < min_degree. 2 removes the loose single-edge pendants
# (tightest look); set to 1 to show every connected node.
# CAPTION when >1: "k-core (degree >= 2) shown; n per panel below"

## hub labelling on the graph: "custom" (explicit list below, +/- hubs),
## "R6" (connector hubs only), "hubs" (all WMDz hubs), or "none".
label_mode <- "R6"    # label the connector hubs only (R6). Every R6 hub also has
##                        WMDz >= 1.5, so the labelled nodes are a strict subset of
##                        the black-ringed set: each labelled node carries a ring.
## Regions named in the manuscript text, used only when label_mode is "custom".
label_regions <- c("AVPV","CA","DMH","HATA","MBO","NOD","PAG","PD","PMd",
                   "PPN","PRE","PVH","RN","SAG","SCs","SNr","VTA")
label_hubs_too <- TRUE   # also label R6 connector hubs (network-defined key nodes)
max_labels <- 10         # cap per panel. Priority: text-named AND hub, then R6
# hubs, then remaining text-named regions, each by PC.
## Always-labelled regions: forced into every panel's label set BEFORE the cap is
## filled, so they appear in every condition even when their PC/hub status would
## otherwise drop them. Total labels are still capped at max_labels (a pinned
## region displaces the lowest-priority auto-pick). Empty c() to disable.
pin_regions <- c("ACB")  # nucleus accumbens, labelled in all conditions
wmdz_hub_cut <- 1.5

## encoding scales (fixed across panels so they are comparable)
zmax_fill  <- 2.5             # WMDz colour clamp (symmetric; white at 0)
pc_size    <- c(base = 2, slope = 8)          # node size = base + slope * PC
edge_r_rng <- c(0.60, 1.00)                   # r range mapped to edge width
edge_w_rng <- c(0.10, 0.90)                   # edge width range (pt) - hairline
uniform_edges <- TRUE                         # TRUE: single hairline. At 10% density
# every kept edge is r>0.8 (Dependent
# >0.93), so width~r barely varies and
# just adds noise. FALSE: width = r.
edge_w_fixed  <- 0.35                          # hairline width when uniform_edges
edge_col   <- "#B9BDC2"                       # light grey edges
edge_alpha <- 0.30                            # edge transparency
frame_w    <- 1.4                             # module-ring line width
lab_cex    <- 1.0                            # region labels = 8 pt (device pointsize = 8, rendered
##                                              at final 69 mm and placed 1:1, so cex 1 == 8 pt).
##                                              Does NOT affect node or legend-dot sizes.
layout_type <- "kk"                          # "kk" (compact, no spokes; default),
#  "drl", or "fr"
## Radial tightening: nodes beyond this radius percentile are pulled inward so the
## informative core fills the panel (the sparse outliers otherwise shrink it after
## rescale). 1 = off; 0.90 keeps the inner 90% fixed and soft-compresses the outer
## 10%. Lower = tighter core.
outlier_pull <- 0.90
## Panel fill: scale the layout so the BULK of nodes (98th pct radius) reach this
## fraction of the panel edge; the few beyond are clamped to the border rather
## than shrinking everything. Higher = less white margin. ~0.95 halves the border.
panel_fill <- 0.99

treatments <- c("Vehicle", "Acute Morphine", "Morphine-Dependent", "Ro 64-6198")
sexes      <- c("Combined", "Males", "Females")
trt_abbr   <- c("Vehicle" = "Veh", "Acute Morphine" = "AcM",
                "Morphine-Dependent" = "MDep", "Ro 64-6198" = "Ro")
## sheet-name pieces in 01_Connectivity.xlsx  ("Corr <ct> <sx>")
corr_trt   <- c("Vehicle" = "Veh", "Acute Morphine" = "AcMor",
                "Morphine-Dependent" = "MorDep", "Ro 64-6198" = "Ro")
corr_sex   <- c("Combined" = "C", "Males" = "M", "Females" = "F")

## module ring palette (used only by color_mode "module_wmdz")
ring_pal <- c("#D62728", "#1F77B4", "#2CA02C", "#9467BD", "#FF7F0E",
              "#8C564B", "#17BECF", "#BCBD22", "#7F7F7F")

## ---- NODE COLOURING ------------------------------------------------------
## "division_super" (default): fill = CCFv3 major division (7 Allen super-groups),
##                             black ring = hub. Same colour = same system set-wide.
## "super_module" : fill = major division, ring = within-panel module
##                  (anatomy AND community structure in one node).
## "division_fill": fill = all 13 CCFv3 divisions, black ring = hub (high-res; busy).
## "group_fill"   : fill = module's dominant major division, black ring = hub.
## "group_ring"   : ring = module's dominant major division, fill = WMDz (Kimbrough).
## "module_wmdz"  : ring = raw within-panel module, fill = WMDz (original).
color_mode <- "division_super" # fill = anatomy (major division), black ring = hub.
# (Module ring dropped: labels were confusing and the
#  ring colours were arbitrary/unmatched across panels.)
module_labels <- FALSE        # off (see above)
show_Q        <- FALSE        # module structure is soft (nothing survived FDR); the
# graphs show anatomy + hubs, not the module story.

## Division grouping and palette come from divisions.R -- the single definition
## shared by every figure script -- so a division is the SAME colour here as in
## Figures 2, 3, 4A/B and S4. (This file previously carried its own palette,
## which did not match them.) The names below are aliases, so the drawing code
## downstream is unchanged.
source(file.path(SCRIPT_DIR, "divisions.R"))
div_order <- division_order    # 13 CCFv3 divisions, canonical order
div_col   <- division_colors   # 13 divisions -> hex (used by color_mode "division_fill")

## 13 CCFv3 divisions -> 7 Allen MAJOR divisions (legible at node size; same
## ontology, one level up). Amygdala -> Cerebral nuclei (Allen splits it between
## cortex and nuclei; one judgment call, note in caption). super_order comes
## from divisions.R unchanged.
super_div <- super_map         # 13 divisions -> 7 super-divisions
super_col <- super_colors      # 7 super-divisions -> hex (used by "division_super")

## Guimera-Amaral role order + colours (cool non-hubs, warm hubs)
role_levels <- c("R1 ultra-peripheral", "R2 peripheral", "R3 non-hub connector",
                 "R4 non-hub kinless",  "R5 provincial hub", "R6 connector hub",
                 "R7 kinless hub")
role_short  <- c("R1", "R2", "R3", "R4", "R5", "R6", "R7")
role_col    <- c("#DEEBF7", "#9ECAE1", "#4292C6", "#08519C",
                 "#FDAE6B", "#E6550D", "#A63603")
names(role_col) <- role_levels

## ----------------------------- LOAD ----------------------------------------
nodes <- as.data.frame(read_excel(roles_file, sheet = "Node_Roles"))
nodes$WM_strength_z         <- as.numeric(nodes$WM_strength_z)
nodes$Participation_coef_wt <- as.numeric(nodes$Participation_coef_wt)
nodes$Network_module        <- as.character(nodes$Network_module)

## module composition -> dominant functional group per (Treatment, Sex, module)
comp <- as.data.frame(read_excel(roles_file, sheet = "Module_Composition"))
comp$Pct            <- as.numeric(comp$Pct)
comp$Network_module <- as.character(comp$Network_module)
comp$Group          <- super_div[comp$Structure]
dom_group <- function(tr, sx) {
  d  <- comp[comp$Treatment == tr & comp$Sex == sx & !is.na(comp$Group), ]
  if (nrow(d) == 0) return(setNames(character(0), character(0)))
  ag <- aggregate(Pct ~ Network_module + Group, data = d, FUN = sum)
  sp <- split(ag, ag$Network_module)
  vapply(sp, function(x) x$Group[which.max(x$Pct)], character(1))
}

read_corr <- function(tr, sx) {
  sheet <- sprintf("Corr %s %s", corr_trt[[tr]], corr_sex[[sx]])
  m  <- as.data.frame(read_excel(conn_file, sheet = sheet))
  rn <- m[[1]]; m[[1]] <- NULL
  M  <- as.matrix(m); rownames(M) <- rn; colnames(M) <- colnames(m)
  storage.mode(M) <- "double"
  M
}

## -------- threshold_to_density(): verbatim from the network pipeline -------
threshold_to_density <- function(cor_mat, density) {
  cm <- cor_mat; cm[!is.finite(cm)] <- 0; diag(cm) <- 0
  pos <- cm; pos[pos < 0] <- 0
  ut         <- upper.tri(pos)
  n_possible <- sum(ut)
  vals       <- pos[ut]
  pos_vals   <- vals[vals > 0]
  n_keep     <- min(ceiling(density * n_possible), length(pos_vals))
  if (n_keep < 1 || length(pos_vals) == 0)
    return(matrix(0, nrow(pos), ncol(pos), dimnames = dimnames(pos)))
  thr    <- sort(pos_vals, decreasing = TRUE)[n_keep]
  thresh <- pos; thresh[thresh < thr] <- 0
  thresh
}

## --------------------------- ENCODING HELPERS ------------------------------
wmdz_fill <- function(z) {
  ramp <- colorRamp(c("#2166AC", "#F7F7F7", "#B2182B"))
  out  <- character(length(z))
  for (i in seq_along(z)) {
    if (is.na(z[i])) { out[i] <- "#B8B8B8"; next }
    zc <- max(min(z[i], zmax_fill), -zmax_fill)
    t  <- (zc + zmax_fill) / (2 * zmax_fill)
    rg <- ramp(t); out[i] <- rgb(rg[1], rg[2], rg[3], maxColorValue = 255)
  }
  out
}
edge_width <- function(r) {
  if (uniform_edges) return(rep(edge_w_fixed, length(r)))
  rc <- pmax(pmin(r, edge_r_rng[2]), edge_r_rng[1])
  edge_w_rng[1] + (rc - edge_r_rng[1]) / diff(edge_r_rng) * diff(edge_w_rng)
}
node_size <- function(pc) pc_size["base"] + pc_size["slope"] * ifelse(is.na(pc), 0, pc)

## ------------------------ BUILD ONE GRAPH ----------------------------------
build_graph <- function(tr, sx) {
  M   <- read_corr(tr, sx)
  adj <- threshold_to_density(M, edge_density)
  g   <- graph_from_adjacency_matrix(adj, mode = "undirected",
                                     weighted = TRUE, diag = FALSE)
  g   <- delete_vertices(g, V(g)[degree(g) == 0])          # drop isolates
  if (giant_only && vcount(g) > 0) {                       # keep largest component
    comp <- components(g)
    g <- induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
  }
  if (min_degree > 1 && vcount(g) > 0) {                   # k-core: drop loose pendants
    g <- induced_subgraph(g, which(coreness(g) >= min_degree))
    if (giant_only && vcount(g) > 0) {                     # k-core may disconnect
      comp <- components(g)
      g <- induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
    }
  }
  nd  <- nodes[nodes$Treatment == tr & nodes$Sex == sx, ]
  mt  <- match(V(g)$name, nd$Region)
  V(g)$PC    <- nd$Participation_coef_wt[mt]
  V(g)$WMDz  <- nd$WM_strength_z[mt]
  V(g)$mod   <- nd$Network_module[mt]
  V(g)$role  <- nd$Node_role[mt]
  V(g)$div   <- nd$Structure[mt]
  V(g)$super <- unname(super_div[nd$Structure[mt]])
  dg <- dom_group(tr, sx)                                   # module -> group
  V(g)$group <- unname(dg[as.character(V(g)$mod)])
  g
}

## white-halo text, centred at (x,y) (8-way offset white, then coloured on top)
halo_text <- function(x, y, labels, cex, halo = "white", col = "black") {
  d <- strwidth("o", cex = cex) * 0.18
  for (a in seq(0, 2 * pi, length.out = 9)[-1])
    text(x + cos(a) * d, y + sin(a) * d, labels, cex = cex, col = halo,
         adj = c(0.5, 0.5))
  text(x, y, labels, cex = cex, col = col, adj = c(0.5, 0.5))
}

## hub labels centred on their node; only labels that would overlap are nudged
## the minimum distance to clear (no leader lines). Lower-priority labels give
## way first, so the most important ones stay dead-centred.
place_net_labels <- function(lay, keep, labs, prio = NULL) {
  if (!any(keep)) return(invisible())
  idx <- which(keep)
  ax <- lay[idx, 1]; ay <- lay[idx, 2]           # node anchors = plotted coords (rescale=FALSE)
  tx <- labs[idx]; n <- length(idx)
  if (is.null(prio)) prio <- seq_len(n) else prio <- prio[idx]
  w  <- strwidth(tx, cex = lab_cex) / 2 + strwidth("n", cex = lab_cex) * 0.35
  h  <- strheight("Mg", cex = lab_cex) * 0.62
  ## start each label pushed off its node, radially outward from the cloud centre,
  ## so it doesn't bury the node and the leader line has a clear direction.
  cx <- mean(lay[, 1]); cy <- mean(lay[, 2])
  dirx <- ax - cx; diry <- ay - cy
  dl <- sqrt(dirx^2 + diry^2); dl[dl < 1e-6] <- 1e-6
  off <- 3.0 * h
  lx <- ax + dirx / dl * off
  ly <- ay + diry / dl * off
  ## separate overlapping label boxes (lower priority moves)
  for (it in 1:600) {
    moved <- FALSE
    if (n > 1) for (i in 1:(n - 1)) for (j in (i + 1):n) {
      dx <- lx[j] - lx[i]; dy <- ly[j] - ly[i]
      ox <- (w[i] + w[j]) - abs(dx); oy <- 2 * h - abs(dy)
      if (ox > 0 && oy > 0) {
        lo <- if (prio[i] <= prio[j]) i else j
        sgn <- if (oy <= ox) c(0, if ((ly[lo] - (ly[i]+ly[j])/2) >= 0) 1 else -1)
        else           c(if ((lx[lo] - (lx[i]+lx[j])/2) >= 0) 1 else -1, 0)
        step <- (if (oy <= ox) oy else ox) + 1e-4
        lx[lo] <- lx[lo] + sgn[1] * step
        ly[lo] <- ly[lo] + sgn[2] * step
        moved <- TRUE
      }
    }
    if (!moved) break
  }
  ## keep every label box inside the panel frame so nothing clips once the cloud
  ## fills the panel (edge-node labels get a short inward leader instead of
  ## shooting past the border). w is per-label half-width; h is the half-height.
  lx <- pmin(pmax(lx, -1 + w), 1 - w)
  ly <- pmin(pmax(ly, -1 + h), 1 - h)
  ## leader line node -> label: white halo under a thin dark line, so it reads over
  ## the grey edges. Text (with its own halo) is drawn last and covers the label end.
  for (k in seq_len(n)) {
    segments(ax[k], ay[k], lx[k], ly[k], col = "white",  lwd = 2.6, lend = 1)
    segments(ax[k], ay[k], lx[k], ly[k], col = "grey20", lwd = 0.8, lend = 1)
    points(ax[k], ay[k], pch = 21, bg = "grey20", col = "white", cex = 0.5, lwd = 0.4)
  }
  for (k in seq_len(n)) halo_text(lx[k], ly[k], tx[k], cex = lab_cex)
}

## ------------------------ DRAW ONE NETWORK PANEL ---------------------------
draw_network <- function(g) {
  hub <- !is.na(V(g)$WMDz) & V(g)$WMDz >= wmdz_hub_cut
  
  ## fill + ring depend on colour mode
  ml <- sort(unique(V(g)$mod[!is.na(V(g)$mod)]))
  rc <- setNames(ring_pal[seq_along(ml)], ml)               # raw-module ring colours
  if (color_mode == "division_super") {        # default: 7 Allen major divisions
    fill  <- ifelse(is.na(V(g)$super), "#C8C8C8", super_col[V(g)$super])
    frame <- ifelse(hub, "black", "white"); fw <- ifelse(hub, 2.0, 0.5)
  } else if (color_mode == "super_module") {   # fill = division, ring = module
    fill  <- ifelse(is.na(V(g)$super), "#C8C8C8", super_col[V(g)$super])
    frame <- ifelse(is.na(V(g)$mod), "#B8B8B8", rc[as.character(V(g)$mod)]); fw <- frame_w
  } else if (color_mode == "division_fill") {  # all 13 divisions
    fill  <- ifelse(is.na(V(g)$div), "#C8C8C8", div_col[V(g)$div])
    frame <- ifelse(hub, "black", "white"); fw <- ifelse(hub, 2.0, 0.5)
  } else if (color_mode == "group_fill") {     # module dominant division
    fill  <- ifelse(is.na(V(g)$group), "#C8C8C8", super_col[V(g)$group])
    frame <- ifelse(hub, "black", "white"); fw <- ifelse(hub, 2.0, 0.5)
  } else if (color_mode == "group_ring") {     # module dominant division as ring
    fill  <- wmdz_fill(V(g)$WMDz)
    frame <- ifelse(is.na(V(g)$group), "#C8C8C8", super_col[V(g)$group]); fw <- frame_w
  } else {                                     # module_wmdz (raw modules)
    fill  <- wmdz_fill(V(g)$WMDz)
    frame <- ifelse(is.na(V(g)$mod), "#B8B8B8", rc[as.character(V(g)$mod)]); fw <- frame_w
  }
  
  ## which nodes get a label
  keep <- rep(FALSE, vcount(g))
  is_r6 <- !is.na(V(g)$role) & V(g)$role == "R6 connector hub"
  pcv  <- ifelse(is.na(V(g)$PC), 0, V(g)$PC)
  label_imp <- pcv                                   # default importance = PC
  if (label_mode == "hubs") keep <- hub
  else if (label_mode == "R6") keep <- is_r6
  else if (label_mode == "custom") {
    is_txt <- V(g)$name %in% label_regions
    ## priority tiers: 3 = text-named & hub, 2 = R6 hub, 1 = text-named, 0 = none
    tier <- ifelse(is_txt & is_r6 & label_hubs_too, 3L,
                   ifelse(is_r6 & label_hubs_too, 2L,
                          ifelse(is_txt, 1L, 0L)))
    label_imp <- tier * 100 + pcv                    # higher = more important
    ## pinned regions are always labelled wherever present, before the cap is
    ## filled (see pin_regions) -- guarantees e.g. ACB in every condition.
    pin  <- which(V(g)$name %in% pin_regions)
    keep[pin] <- TRUE
    cand <- setdiff(which(tier > 0), pin)            # remaining candidates
    cand <- cand[order(-tier[cand], -pcv[cand])]     # tier desc, then PC desc
    keep[head(cand, max(0L, max_labels - length(pin)))] <- TRUE  # fill to the cap
  }
  
  set.seed(layout_seed)
  ## Kamada-Kawai is run UNWEIGHTED on purpose: igraph reads layout weights as
  ## desired edge LENGTHS, so passing r would push the most strongly co-active
  ## pairs furthest apart. The drl / fr alternatives take weights in the usual
  ## attraction sense and are passed them.
  lay <- switch(layout_type,
                kk  = layout_with_kk(g),
                drl = layout_with_drl(g, weights = E(g)$weight),
                layout_with_fr(g, weights = E(g)$weight))
  
  ## tighten sparse outliers so the informative core fills the panel (outlier_pull)
  if (outlier_pull < 1 && nrow(lay) > 3) {
    ctr <- colMeans(lay)
    rel <- sweep(lay, 2, ctr)                         # centre
    rad <- sqrt(rowSums(rel^2))
    cap <- stats::quantile(rad, outlier_pull)
    scl <- ifelse(rad > cap & rad > 0,
                  (cap + (rad - cap) * 0.35) / rad, 1) # soft-compress the excess
    lay <- sweep(rel * scl, 2, ctr, FUN = "+")        # restore centre
  }
  
  ## Fill the panel: rescale EACH axis independently to [-1, 1] so the cloud fills
  ## the square frame (this is exactly what igraph rescale = TRUE does). It is done
  ## here on 'lay' -- not via rescale = TRUE in plot() -- so the label overlay reads
  ## the same coordinates and stays matched. panel_fill (~0.99) leaves a hair of
  ## inset so edge node bodies + haloed labels sit just inside the frame.
  if (nrow(lay) > 3) {
    for (j in 1:2) {
      rj <- range(lay[, j])
      lay[, j] <- if (diff(rj) > 0) ((lay[, j] - rj[1]) / diff(rj) * 2 - 1) * panel_fill else 0
    }
  }
  
  ## z-order: draw hub (black-ring) nodes last so they sit on top of neighbours.
  ## Layout is computed above on the original graph, then everything is reordered
  ## together (non-hubs first, hubs last) so the picture is identical bar z-order.
  zo <- order(hub)
  perm <- integer(vcount(g)); perm[zo] <- seq_along(zo)
  g   <- permute(g, perm)
  lay <- lay[zo, , drop = FALSE]
  fill <- fill[zo]; frame <- frame[zo]
  if (length(fw) > 1) fw <- fw[zo]
  keep <- keep[zo]; label_imp <- label_imp[zo]
  
  par(family = font_fam, mar = c(0.2, 0.2, 0.2, 0.2))
  plot(g, layout = lay,
       vertex.size        = node_size(V(g)$PC),
       vertex.color       = fill,
       vertex.frame.color = frame,
       vertex.frame.width = fw,
       vertex.label       = NA,
       edge.width         = edge_width(E(g)$weight),
       edge.color         = adjustcolor(edge_col, alpha.f = edge_alpha),
       rescale = FALSE, xlim = c(-1, 1), ylim = c(-1, 1))
  
  ## overlay node labels; centred on their node, only colliding ones nudged
  place_net_labels(lay, keep, V(g)$name, prio = label_imp)
  
  ## module labels: number each module M1..Mk at its centroid, coloured to match
  ## the module ring, de-overlapped. Numeric (not anatomical) because these
  ## communities are anatomically mixed, not pure divisions; node fill already
  ## carries anatomy. Labels are arbitrary within a panel, not matched across.
  if (module_labels) {
    nx <- (lay[, 1] - min(lay[, 1])) / (max(lay[, 1]) - min(lay[, 1])) * 2 - 1
    ny <- (lay[, 2] - min(lay[, 2])) / (max(lay[, 2]) - min(lay[, 2])) * 2 - 1
    mm  <- V(g)$mod
    mods <- sort(unique(mm[!is.na(mm)]))
    ## order module numbers by size (largest = M1) for a stable within-panel scheme
    msize <- vapply(mods, function(k) sum(!is.na(mm) & mm == k), numeric(1))
    ord <- order(-msize); mods <- mods[ord]
    mlab <- paste0("M", seq_along(mods))                # M1 (largest) .. Mk
    mcx <- vapply(mods, function(k) mean(nx[!is.na(mm) & mm == k]), numeric(1))
    mcy <- vapply(mods, function(k) mean(ny[!is.na(mm) & mm == k]), numeric(1))
    mcol <- rc[as.character(mods)]                      # module ring colours
    keepm <- mlab != ""
    mcx <- mcx[keepm]; mcy <- mcy[keepm]; mlab <- mlab[keepm]; mcol <- mcol[keepm]
    if (length(mlab)) {
      hw <- strwidth(mlab, cex = lab_cex * 1.25) / 2
      hh <- strheight("Mg", cex = lab_cex * 1.25) * 0.9
      for (it in 1:500) {                               # 2D de-overlap of module tags
        moved <- FALSE; n <- length(mcx)
        if (n > 1) for (i in 1:(n - 1)) for (j in (i + 1):n) {
          dx <- mcx[j] - mcx[i]; dy <- mcy[j] - mcy[i]
          ox <- (hw[i] + hw[j]) - abs(dx); oy <- (2 * hh) - abs(dy)
          if (ox > 0 && oy > 0) {
            if (oy <= ox) { s <- (oy/2 + 1e-4) * if (dy >= 0) 1 else -1
            mcy[i] <- mcy[i]-s; mcy[j] <- mcy[j]+s }
            else          { s <- (ox/2 + 1e-4) * if (dx >= 0) 1 else -1
            mcx[i] <- mcx[i]-s; mcx[j] <- mcx[j]+s }
            moved <- TRUE
          }
        }
        if (!moved) break
      }
      for (k in seq_along(mlab))
        halo_text(mcx[k], mcy[k], mlab[k], cex = lab_cex * 1.25, col = mcol[k])
    }
  }
  
  ## descriptive modularity annotation (Q + module count; no significance claim)
  if (show_Q) {
    mem <- V(g)$mod
    if (any(is.na(mem)))
      mem[is.na(mem)] <- max(mem, na.rm = TRUE) + seq_len(sum(is.na(mem)))
    Qv <- modularity(g, membership = as.integer(factor(mem)),
                     weights = E(g)$weight)
    nmod <- length(unique(V(g)$mod[!is.na(V(g)$mod)]))
    usr <- par("usr")
    text(usr[1] + 0.02 * (usr[2] - usr[1]), usr[4] - 0.02 * (usr[4] - usr[3]),
         sprintf("Q = %.2f  (%d modules)", Qv, nmod),
         adj = c(0, 1), cex = lab_cex * 1.15, col = "grey20")
  }
}

## ------------------------ GENERATE 12 NETWORK PANELS -----------------------
for (sx in sexes) for (tr in treatments) {
  g <- build_graph(tr, sx)
  stub <- file.path(out_dir, sprintf("network_%s_%s", sx, trt_abbr[[tr]]))
  tiff(paste0(stub, ".tiff"), width = panel_in, height = panel_in, units = "in",
       res = dev_res, compression = "lzw", family = font_fam, pointsize = 8,
       type = dev_type); draw_network(g); dev.off()
  png(paste0(stub, ".png"),  width = panel_in, height = panel_in, units = "in",
      res = dev_res, family = font_fam, pointsize = 8,
      type = dev_type); draw_network(g); dev.off()
  cat(sprintf("  %-9s %-18s  nodes=%d  edges=%d\n",
              sx, tr, vcount(g), ecount(g)))
}

## ------------------------ SHARED LEGEND / KEY ------------------------------
draw_legend <- function() {
  par(family = font_fam, mar = c(0.5, 0.5, 0.5, 0.5))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  
  if (color_mode == "division_fill") {
    text(0.02, 0.97, "CCFv3 Division (Node Fill)", adj = 0, font = 2, cex = 0.85)
    nL <- length(div_order); half <- ceiling(nL / 2)
    for (i in seq_len(nL)) {
      col_i <- if (i <= half) 0L else 1L
      row_i <- if (i <= half) i else i - half
      xx <- 0.05 + col_i * 0.50
      yy <- 0.90 - (row_i - 1) * 0.060
      points(xx, yy, pch = 21, cex = 1.9, bg = div_col[div_order[i]], col = "white", lwd = 1)
      text(xx + 0.045, yy, div_order[i], adj = 0, cex = 0.62)
    }
    yh <- 0.90 - half * 0.060
    points(0.05, yh, pch = 21, cex = 1.9, bg = "grey80", col = "black", lwd = 2.2)
    text(0.095, yh, "Hub (WMDz \u2265 1.5)", adj = 0, cex = 0.66)
    ybase <- yh - 0.08
  } else if (color_mode %in% c("division_super", "super_module", "group_fill", "group_ring")) {
    ttl <- switch(color_mode,
                  division_super = "Major Division (Node Fill)",
                  super_module   = "Major Division (Node Fill)",
                  group_fill     = "Module Identity (Node Fill)",
                  group_ring     = "Module Identity (Node Ring)")
    text(0.02, 0.985, sub(" \\(", "\n(", ttl), adj = c(0, 1), font = 2, cex = 8/6)  # 8 pt; wrapped before "(" to fit the 36 mm channel
    ring_swatch <- (color_mode == "group_ring")
    for (i in seq_along(super_order)) {
      yy <- 0.83 - (i - 1) * 0.060
      if (ring_swatch)
        points(0.05, yy, pch = 21, cex = 3.0, bg = "grey90", col = super_col[super_order[i]], lwd = 2.6)
      else
        points(0.05, yy, pch = 21, cex = 3.0, bg = super_col[super_order[i]], col = "white", lwd = 1)
      text(0.10, yy, super_order[i], adj = 0, cex = 8/6)
    }
    yh <- 0.83 - length(super_order) * 0.060 - 0.022   # extra gap so the thick hub ring
    ##                                                        clears the row above (ring size unchanged)
    if (color_mode %in% c("division_super", "group_fill")) {
      points(0.05, yh, pch = 21, cex = 3.0, bg = "grey80", col = "black", lwd = 2.2)
      text(0.10, yh, "Hub (WMDz \u2265 1.5)", adj = 0, cex = 8/6)
      ybase <- yh - 0.07
    } else if (color_mode == "super_module") {
      text(0.04, yh, "ring colour = module (within panel)", adj = 0, cex = 0.64, font = 3)
      ybase <- yh - 0.07
    } else {                                    # group_ring: fill encodes WMDz
      text(0.02, yh - 0.01, "Node fill: WMDz", adj = 0, font = 2, cex = 0.8)
      zs <- seq(-zmax_fill, zmax_fill, length.out = 60); cz <- wmdz_fill(zs)
      for (i in seq_along(zs))
        rect(0.05 + (i - 1)/60*0.5, yh-0.10, 0.05 + i/60*0.5, yh-0.05, col = cz[i], border = NA)
      text(c(0.05,0.30,0.55), yh-0.13, c("-2.5","0","2.5"), cex = 0.6)
      ybase <- yh - 0.18
    }
  } else {
    ## raw-module legend (module_wmdz)
    text(0.02, 0.95, "Module (within panel)", adj = 0, font = 2, cex = 0.85)
    for (i in 1:3) {
      points(0.06 + (i - 1) * 0.07, 0.88, pch = 21, cex = 2.4, bg = "white",
             col = ring_pal[i], lwd = 2)
      text(0.06 + (i - 1) * 0.07, 0.81, c("A","B","C")[i], cex = 0.7)
    }
    text(0.02, 0.68, "Node fill: WMDz", adj = 0, font = 2, cex = 0.85)
    zs <- seq(-zmax_fill, zmax_fill, length.out = 60); cz <- wmdz_fill(zs)
    for (i in seq_along(zs))
      rect(0.05 + (i - 1)/60*0.5, 0.58, 0.05 + i/60*0.5, 0.64, col = cz[i], border = NA)
    text(c(0.05,0.30,0.55), 0.545, c("-2.5","0","2.5"), cex = 0.65)
    ybase <- 0.42
  }
  
  ## PC size (shared)
  text(0.02, ybase, "Participation\nCoefficient (PC)", adj = c(0, 1), font = 2, cex = 8/6)  # wrapped to fit channel
  pcv <- c(0.0, 0.3, 0.6)
  for (i in seq_along(pcv)) {
    points(0.08 + (i - 1) * 0.15, ybase - 0.13, pch = 21, bg = "grey85",
           col = "grey40", cex = node_size(pcv[i]) / 2.55)   # dots sized to the 69 mm render's nodes (PC 0.6 dot == PC-0.6 node)
    text(0.08 + (i - 1) * 0.15, ybase - 0.19, pcv[i], cex = 8/6)  # number tucked under its dot (8 pt)
  }
  ## edge key intentionally omitted: edges are a single uniform hairline, so the
  ## meaning (positive co-activations at 10% density; width not mapped to r) is
  ## stated in the figure caption rather than shown as a swatch.
}
lstub <- file.path(out_dir, "network_legend")
# Render once, trim the surrounding white margin to the content, add a small
# uniform pad, then write both formats so the legend drops in compact.
.leg_tmp <- tempfile(fileext = ".png")
png(.leg_tmp, width = 2.0, height = 2.4, units = "in",
    res = dev_res, family = font_fam, pointsize = 6,
    type = dev_type); draw_legend(); dev.off()
.leg <- magick::image_trim(magick::image_read(.leg_tmp))
.leg <- magick::image_border(.leg, "white", sprintf("%dx%d", round(0.03*dev_res), round(0.03*dev_res)))
magick::image_write(.leg, paste0(lstub, ".png"),  format = "png",
                    density = as.character(dev_res))
magick::image_write(.leg, paste0(lstub, ".tiff"), format = "tiff",
                    compression = "lzw", density = as.character(dev_res))

## ============================================================================
## COMPANION: CARTOGRAPHIC ROLE-COMPOSITION BARS  (one panel per sex)
## Stacked proportion of nodes in each Guimera-Amaral role, across conditions.
## ============================================================================
draw_roles <- function(sx) {
  present <- role_levels[role_levels %in%
                           unique(nodes$Node_role[nodes$Sex == sx & !is.na(nodes$Node_role)])]
  ## proportion matrix: roles (rows) x conditions (cols)
  P <- matrix(0, nrow = length(present), ncol = length(treatments),
              dimnames = list(present, trt_abbr[treatments]))
  for (j in seq_along(treatments)) {
    rr <- nodes$Node_role[nodes$Sex == sx & nodes$Treatment == treatments[j] &
                            !is.na(nodes$Node_role)]
    if (length(rr) > 0) {
      tb <- table(factor(rr, levels = present))
      P[, j] <- as.numeric(tb) / sum(tb)
    }
  }
  par(family = font_fam, mar = c(3.4, 4.2, 1.2, 0.8))
  ## Rendered at final size with pointsize = 8, so cex 1 == 8 pt: every text
  ## element in this panel sits on the same >= 8 pt floor as the other figures.
  bp <- barplot(P, col = role_col[present], border = "white",
                ylim = c(0, 1), las = 1, ylab = "Proportion of nodes",
                cex.names = 1, cex.lab = 1)
  legend("topright", inset = c(-0.0, -0.02), legend = role_short[match(present, role_levels)],
         fill = role_col[present], border = "white", bty = "n",
         cex = 1, ncol = 1, xpd = NA)
}
for (sx in sexes) {
  stub <- file.path(out_dir, sprintf("roles_%s", sx))
  tiff(paste0(stub, ".tiff"), width = 4.4, height = 3.6, units = "in",
       res = dev_res, compression = "lzw", family = font_fam, pointsize = 8,
       type = dev_type); draw_roles(sx); dev.off()
  png(paste0(stub, ".png"),  width = 4.4, height = 3.6, units = "in",
      res = dev_res, family = font_fam, pointsize = 8,
      type = dev_type); draw_roles(sx); dev.off()
}

cat("\nWrote 12 network panels + legend + 3 role-composition panels to:\n  ", out_dir, "\n")