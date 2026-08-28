# network_analysis.R
# =============================================================================
# Whole-brain cFos co-activation network analysis.
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
#   Rscript scripts/network_analysis.R
# Input and output locations are set in the Config section below; see README.
# =============================================================================
#
# Raw cFos counts  ->  within-module degree z-score (WMDz) + participation
# coefficient (PC), plus supporting network metrics. No figures.
# Outputs: topic-grouped Excel workbooks.
#
# PIPELINE (each step literature-anchored):
#   1. Load raw counts (Clustered Males / Clustered Females sheets).
#   2. Clean per region: negatives -> NA; zeros -> half-min-nonzero imputation
#      (standard in cFos/omics); log10 transform. Primary networks are sex-stratified
#      (Males / Females); the pooled Combined network mean-centers each region within
#      sex before pooling (residualizes sex) and is used for resampling-based inference.
#   3. Region x region Pearson correlation across animals (co-activation).
#   4. Density-matched (proportional) thresholding on POSITIVE edges to a COMMON
#      density across conditions, so degree/modularity/roles are comparable
#      (Rubinov & Sporns 2010; van den Heuvel et al. 2017; positive-only:
#      Brynildsen et al. 2020; Xie et al. 2024).
#   5. Network modules via CONSENSUS LOUVAIN on the thresholded graph itself
#      (1000 runs; Lancichinetti & Fortunato 2012) -> stable, graph-native
#      partition for the cartographic metrics.
#   6. WMDz + PC on (graph, consensus modules) per Guimera & Amaral (2005);
#      hub cut = wmdz_hub_cut (1.5; rationale in the Config block); standard P
#      cut points for the 7 node roles. Small modules flagged (size bias:
#      Pedersen et al. 2020).
#   7. Robustness: global metrics across a density sweep + AUC (Garrison et al.
#      2015; Hallquist & Hillary 2019).
#   8. Modularity significance (label-shuffle permutation, density-matched).
#   9. Connectivity differences between conditions (Fisher r-to-z, BH-FDR),
#      with per-pair r difference folded in.
#  10. Cross-treatment module stability (on the consensus network modules).
#  11. Inferential stats on network metrics (Combined): bootstrap CIs per
#      condition + label-permutation tests of between-condition differences,
#      so metric changes are testable, not just descriptive (Vetere 2017).
# =============================================================================

# ---- Packages ----------------------------------------------------------------
required_pkgs <- c("tidyverse", "readxl", "openxlsx", "igraph")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages. Please run:\n  install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}
library(tidyverse)
library(readxl)
library(openxlsx)
library(igraph)

# ================== CONFIG ==================
# Paths are relative to the repository root. Run from the repo root, e.g.:
#   Rscript scripts/network_analysis.R
# Override by editing the two directories below or by setting the environment
# variables CFOS_DATA_DIR / CFOS_OUTPUT_DIR before running.
DATA_DIR    <- path.expand(Sys.getenv("CFOS_DATA_DIR",   "data"))
OUTPUT_DIR  <- path.expand(Sys.getenv("CFOS_OUTPUT_DIR", "results"))
SCRIPT_DIR  <- path.expand(Sys.getenv("CFOS_SCRIPT_DIR", "scripts"))
file_path   <- file.path(DATA_DIR, "Drug_Composite.xlsx")
output_root <- OUTPUT_DIR

treatments <- c("Vehicle", "Acute Morphine", "Morphine-Dependent", "Ro 64-6198")
sex_groups <- c("Males", "Females", "Combined")  # Males/Females = primary (sex-stratified,
# matching Ardinger 2024 / Bloch 2024); Combined = pooled network retained for resampling-based
# inference only, where per-sex n (5-6) is too low for stable bootstrap/permutation.

# Network thresholding / module / role parameters
primary_density <- 0.10                        # common density for primary analysis
sweep_densities <- seq(0.02, 0.20, by = 0.02)  # robustness sweep (2%-20%)
n_consensus     <- 1000                        # consensus Louvain runs (lower for a faster test run)
min_module_size <- 5                           # role-reliability flag threshold
n_perm          <- 5000                        # modularity permutations (lower for a faster test run)
hub_sd_mult     <- 1                           # "high" on a centrality = > mean + this*SD
hub_min_metrics <- 2                           # hub = high on >= this many of 4 centralities
n_boot          <- 1000                        # bootstrap / permutation iterations (lower for a faster test run)
# Cartographic (Guimera-Amaral) hub cut, applied to WEIGHTED WMDz. WM_strength_z here is a
# proper within-module z-score (mean 0, SD ~1), so the canonical "high WMDz" range applies:
# 1.5 is the lower bound of the 1.5-2.5 range (Cohen & D'Esposito, cited by Kimbrough 2020) and
# flags ~4% of nodes (~8 hubs/network) on these data. NOTE: Kimbrough's published 0.8 was on a
# compressed WMDz (range +/-1.5) and is NOT comparable to a unit-SD z; 0.8 here over-calls (~25%
# of nodes). Use 2.0 for a stricter set, 1.28 (~top decile within module) for a broader one.
wmdz_hub_cut    <- 1.5

dir.create(output_root, showWarnings = FALSE, recursive = TRUE)

# Pairwise condition comparisons for Fisher Z connectivity differences
comparisons <- list(
  c("Vehicle",        "Acute Morphine"),
  c("Vehicle",        "Ro 64-6198"),
  c("Vehicle",        "Morphine-Dependent"),
  c("Acute Morphine", "Morphine-Dependent"),
  c("Acute Morphine", "Ro 64-6198")
)

# ================== ANATOMICAL DIVISION LOOKUP ==================
# Region acronym -> CCFv3 division (annotation only). Comes from divisions.R,
# the single definition shared with every figure script, which reads the
# canonical region_division_lookup.csv; struct_of() returns "Other" for any
# acronym absent from that table. This file previously carried its own inline
# map, which was missing 16 of the analysed acronyms (they were typed
# "Unknown") and could drift from the figures.
source(file.path(SCRIPT_DIR, "divisions.R"))

# ================== HELPER: LOAD & PARSE ONE SHEET ==================
load_sheet <- function(file_path, sheet_name) {
  raw <- read_excel(file_path, sheet = sheet_name) |> as.data.frame()
  names(raw) <- trimws(names(raw))
  if (!"acronym" %in% tolower(names(raw)))
    stop("Cannot find 'acronym' column in sheet: ", sheet_name)
  
  acronym_col <- which(tolower(names(raw)) == "acronym")
  regions <- trimws(as.character(raw[[acronym_col]]))
  
  sample_cols  <- setdiff(seq_len(ncol(raw)), c(1, acronym_col))
  sample_names <- names(raw)[sample_cols]
  keep <- !is.na(sample_names) & sample_names != "" &
    sapply(sample_cols, function(j) !all(is.na(raw[[j]])))
  sample_cols  <- sample_cols[keep]
  sample_names <- sample_names[keep]
  
  df <- data.frame(Region = regions, stringsAsFactors = FALSE)
  for (i in seq_along(sample_cols)) {
    df[[sample_names[i]]] <- suppressWarnings(as.numeric(raw[[sample_cols[i]]]))
  }
  
  col_name_map <- list(
    "Vehicle"            = "Vehicle",
    "Acute Morphine"     = "Morphine",
    "Morphine-Dependent" = "Chronic Morphine",
    "Ro 64-6198"         = "Ro"
  )
  
  result <- list()
  for (trt in treatments) {
    excel_prefix <- col_name_map[[trt]]
    if (excel_prefix == "Morphine") {
      pattern <- "^(?!Chronic\\s)Morphine\\s*-"
      cols <- grep(pattern, sample_names, perl = TRUE, ignore.case = TRUE)
    } else {
      pattern <- paste0("^", gsub(" ", "\\\\s+", excel_prefix), "\\s*-")
      cols <- grep(pattern, sample_names, ignore.case = TRUE)
    }
    if (length(cols) == 0) {
      message("  Warning: no columns matched treatment '", trt, "' in sheet '", sheet_name, "'")
      next
    }
    result[[trt]] <- df[, c(1, cols + 1), drop = FALSE]
  }
  result
}

# ================== HELPER: CLEAN & BUILD MATRIX ==================
# negatives -> NA; zeros -> half-min-nonzero imputation; log10. Returns a
# region x animal matrix. Sex pooling correction is applied separately at
# Combined-network construction (within-sex per-region centering), not here.
clean_and_matrix <- function(df, sex_label = NULL) {
  df_clean <- df %>% mutate(across(-Region, as.numeric))
  
  sample_cols <- names(df_clean)[-1]
  all_na_cols <- sample_cols[sapply(sample_cols, function(col) all(is.na(df_clean[[col]])))]
  if (length(all_na_cols) > 0) df_clean <- df_clean %>% dplyr::select(-all_of(all_na_cols))
  
  # Negatives -> NA; zeros -> half of THAT REGION'S minimum non-zero count; log10.
  # The imputation is per region (per feature, i.e. across the row), which is the
  # standard half-minimum rule and what this file's header and the README sheet in
  # 00_Processed_Data.xlsx describe. An earlier version ran column-wise, i.e. per
  # ANIMAL across all 198 regions, so a zero was replaced by half that animal's
  # smallest count anywhere in the brain -- a value with no relation to the region
  # it was substituted into.
  .vals <- as.matrix(df_clean[, -1, drop = FALSE])
  .vals[!is.na(.vals) & .vals < 0] <- NA
  for (.i in seq_len(nrow(.vals))) {
    .z <- !is.na(.vals[.i, ]) & .vals[.i, ] == 0
    if (any(.z)) {
      .pos <- .vals[.i, ][!is.na(.vals[.i, ]) & .vals[.i, ] > 0]
      .vals[.i, .z] <- if (length(.pos)) min(.pos) / 2 else NA
    }
  }
  .out <- as.data.frame(log10(.vals), stringsAsFactors = FALSE)
  names(.out) <- colnames(.vals)
  df_clean <- cbind(as.data.frame(df_clean)["Region"], .out)
  # NOTE: no per-animal z-score. Sexes are pooled (Combined) by within-sex per-region
  # mean-centering at construction time, which removes the sex mean-difference confound
  # that pooling introduces while preserving the cross-animal covariance the network uses.
  # (sex_label retained for signature compatibility; unused.)
  
  region_mat <- df_clean %>%
    group_by(Region) %>%
    summarise(across(everything(), ~ mean(., na.rm = TRUE)), .groups = "drop") %>%
    column_to_rownames("Region") %>%
    as.matrix()
  
  row_sd <- apply(region_mat, 1, function(r) sd(r, na.rm = TRUE))
  bad_rows <- is.na(row_sd) | row_sd == 0 | !apply(region_mat, 1, function(r) any(is.finite(r)))
  if (any(bad_rows)) {
    message("  Dropping ", sum(bad_rows), " region(s) with no variance/data: ",
            paste(rownames(region_mat)[bad_rows], collapse = ", "))
    region_mat <- region_mat[!bad_rows, , drop = FALSE]
  }
  
  region_mat[!is.finite(region_mat)] <- NA
  for (i in seq_len(nrow(region_mat))) {
    if (any(is.na(region_mat[i, ])))
      region_mat[i, is.na(region_mat[i, ])] <- mean(region_mat[i, ], na.rm = TRUE)
  }
  if (nrow(region_mat) < 2) stop("Fewer than 2 valid regions remain after cleaning.")
  for (i in seq_len(nrow(region_mat))) {
    bad <- !is.finite(region_mat[i, ])
    if (any(bad)) {
      fill <- mean(region_mat[i, !bad], na.rm = TRUE)
      region_mat[i, bad] <- if (is.finite(fill)) fill else 0
    }
  }
  region_mat
}

# ================== HELPER: DENSITY-MATCHED THRESHOLD ==================
# Keep the strongest POSITIVE correlations until `density` x [N*(N-1)/2] edges
# survive. Returns thresholded positive adjacency, weight cutoff, achieved density.
threshold_to_density <- function(cor_mat, density) {
  cm <- cor_mat; cm[!is.finite(cm)] <- 0; diag(cm) <- 0
  pos <- cm; pos[pos < 0] <- 0
  ut         <- upper.tri(pos)
  n_possible <- sum(ut)
  vals       <- pos[ut]
  pos_vals   <- vals[vals > 0]
  n_keep     <- min(ceiling(density * n_possible), length(pos_vals))
  if (n_keep < 1 || length(pos_vals) == 0) {
    thresh <- matrix(0, nrow(pos), ncol(pos), dimnames = dimnames(pos))
    return(list(adj = thresh, threshold = NA_real_, achieved_density = 0))
  }
  thr    <- sort(pos_vals, decreasing = TRUE)[n_keep]
  thresh <- pos; thresh[thresh < thr] <- 0
  list(adj = thresh, threshold = thr,
       achieved_density = sum(thresh[ut] > 0) / n_possible)
}

# ================== HELPER: CONSENSUS LOUVAIN NETWORK MODULES ==================
# Communities detected on the thresholded graph itself. Consensus over n_runs
# weighted Louvain runs (co-classification matrix thresholded at tau, re-clustered
# to convergence; Lancichinetti & Fortunato 2012) addresses Louvain stochasticity
# and the modularity resolution limit. Returns named integer membership.
consensus_louvain <- function(g, n_runs = 1000, tau = 0.5, max_iter = 50, seed = 42) {
  vnames <- V(g)$name
  n <- vcount(g)
  if (n < 2) return(setNames(rep(1L, n), vnames))
  set.seed(seed)
  work <- g
  last_membership <- rep(1L, n)
  for (iter in seq_len(max_iter)) {
    M <- matrix(0L, nrow = n_runs, ncol = n)
    for (r in seq_len(n_runs)) {
      cl <- cluster_louvain(work, weights = E(work)$weight)
      M[r, ] <- as.integer(membership(cl))
    }
    last_membership <- M[1, ]
    Co <- matrix(0, n, n)
    for (r in seq_len(n_runs)) Co <- Co + outer(M[r, ], M[r, ], FUN = "==")
    Co <- Co / n_runs
    offdiag <- Co[upper.tri(Co)]
    if (length(offdiag) == 0 || all(offdiag <= 1e-9 | offdiag >= 1 - 1e-9)) break
    D <- Co; D[D < tau] <- 0; diag(D) <- 0
    work <- graph_from_adjacency_matrix(D, mode = "undirected", weighted = TRUE, diag = FALSE)
    V(work)$name <- vnames
    if (ecount(work) == 0) break
  }
  setNames(as.integer(last_membership), vnames)
}

# ================== HELPER: LINEAGE HIERARCHICAL-CLUSTERING MODULE COUNTS ==================
# Reproduces the standard cFos module-detection method (Kimbrough 2020; Kimbrough eNeuro 2021;
# Ardinger 2024): Euclidean distance between regions' interregional-correlation profiles,
# complete-linkage hierarchical clustering, dendrogram cut at a fraction of the max tree height.
# Returns the module count at each cut height, mirroring the tree-cut robustness panels those
# papers report (Kimbrough Fig 5D; Ardinger Fig 4D). Used only to document, on these data, that
# the lineage method collapses to very few modules on strongly-coactivated networks - the
# justification for adopting consensus Louvain. Computed on the FULL (unthresholded) correlation
# matrix, exactly as in the source papers.
hclust_module_count <- function(cor_mat, cut_fractions = c(0.3, 0.4, 0.5, 0.6, 0.7)) {
  nm <- paste0("HC_cut_", cut_fractions)
  cm <- cor_mat; cm[!is.finite(cm)] <- 0; diag(cm) <- 1
  if (is.null(dim(cm)) || nrow(cm) < 3)
    return(setNames(rep(NA_integer_, length(cut_fractions)), nm))
  d  <- dist(cm, method = "euclidean")   # distance between regions' correlation profiles
  hc <- hclust(d, method = "complete")   # complete linkage, as in Kimbrough 2020 / Ardinger 2024
  hmax <- max(hc$height)
  if (!is.finite(hmax) || hmax <= 0)
    return(setNames(rep(NA_integer_, length(cut_fractions)), nm))
  setNames(vapply(cut_fractions, function(f) length(unique(cutree(hc, h = f * hmax))),
                  integer(1)), nm)
}

# ================== HELPER: GUIMERA-AMARAL NODE ROLES (WMDz + PC) ==================
# WMDz (within-module degree z-score) and PC (participation coefficient) on the
# (graph, consensus-module) pair. The WEIGHTED variants drive Node_role (matching the
# Kimbrough/Ardinger weighted-edge cartography); binary variants are reported alongside.
# z is NA for modules < hard_floor (size 1 undefined, size 2 degenerate); Role_reliable
# flags home modules < min_module_size (small-module instability + PC size bias,
# Pedersen et al. 2020). Hub/non-hub cut = z_cut (weighted WMDz; default wmdz_hub_cut = 1.5).
# Computed on the subgraph of assigned nodes.
compute_node_roles <- function(g, membership, min_module_size = 5, hard_floor = 3,
                               z_cut = wmdz_hub_cut) {
  if (is.null(membership) || vcount(g) == 0) return(NULL)
  vnames <- V(g)$name
  mem    <- membership[vnames]; mem[is.na(mem)] <- 0L
  keep   <- which(mem > 0)
  if (length(keep) < 2) return(NULL)
  g_sub     <- induced_subgraph(g, keep)
  sub_names <- V(g_sub)$name
  sub_mem   <- mem[keep]; names(sub_mem) <- sub_names
  
  A <- as.matrix(igraph::as_adjacency_matrix(g_sub, attr = NULL))
  W <- as.matrix(igraph::as_adjacency_matrix(g_sub, attr = "weight"))
  diag(A) <- 0; diag(W) <- 0
  
  mods      <- sort(unique(sub_mem))
  mod_sizes <- table(sub_mem)
  n_sub     <- length(sub_names)
  k_tot <- rowSums(A); s_tot <- rowSums(W)
  
  z_bin <- rep(NA_real_, n_sub); z_wt <- rep(NA_real_, n_sub)
  P_bin <- rep(NA_real_, n_sub); P_wt <- rep(NA_real_, n_sub)
  names(z_bin) <- names(z_wt) <- names(P_bin) <- names(P_wt) <- sub_names
  
  for (mm in mods) {
    idx   <- which(sub_mem == mm)
    msize <- length(idx)
    kin_bin <- rowSums(A[idx, idx, drop = FALSE])
    kin_wt  <- rowSums(W[idx, idx, drop = FALSE])
    if (msize >= hard_floor) {
      sd_bin <- sd(kin_bin); sd_wt <- sd(kin_wt)
      z_bin[idx] <- if (is.finite(sd_bin) && sd_bin > 0) (kin_bin - mean(kin_bin)) / sd_bin else 0
      z_wt[idx]  <- if (is.finite(sd_wt)  && sd_wt  > 0) (kin_wt  - mean(kin_wt))  / sd_wt  else 0
    }
  }
  for (i in seq_len(n_sub)) {
    if (k_tot[i] > 0) {
      frac_b <- sapply(mods, function(mm) sum(A[i, sub_mem == mm]) / k_tot[i])
      P_bin[i] <- 1 - sum(frac_b^2)
    }
    if (s_tot[i] > 0) {
      frac_w <- sapply(mods, function(mm) sum(W[i, sub_mem == mm]) / s_tot[i])
      P_wt[i] <- 1 - sum(frac_w^2)
    }
  }
  
  home_size <- as.integer(mod_sizes[as.character(sub_mem)])
  reliable  <- home_size >= min_module_size
  
  classify <- function(z, P, z_cut) {
    if (is.na(z) || is.na(P)) return(NA_character_)
    if (z < z_cut) {
      if (P <= 0.05) "R1 ultra-peripheral"
      else if (P <= 0.62) "R2 peripheral"
      else if (P <= 0.80) "R3 non-hub connector"
      else "R4 non-hub kinless"
    } else {
      if (P <= 0.30) "R5 provincial hub"
      else if (P <= 0.75) "R6 connector hub"
      else "R7 kinless hub"
    }
  }
  # Roles driven by the WEIGHTED cartography (WM_strength_z, Participation_coef_wt),
  # matching Kimbrough/Ardinger's weighted-edge implementation; binary variants are
  # still reported for reference. Hub/non-hub boundary = wmdz_hub_cut (Fos-scaled).
  roles <- mapply(classify, z_wt, P_wt, MoreArgs = list(z_cut = z_cut))
  roles[!reliable] <- NA_character_
  
  data.frame(
    Region                = sub_names,
    Network_module        = as.integer(sub_mem),
    Home_module_size      = home_size,
    WM_degree_z           = round(z_bin, 4),
    Participation_coef    = round(P_bin, 4),
    WM_strength_z         = round(z_wt, 4),
    Participation_coef_wt = round(P_wt, 4),
    Node_role             = roles,
    Role_reliable         = reliable,
    stringsAsFactors      = FALSE
  )
}

# ================== HELPER: MULTI-CENTRALITY HUB SCORING ==================
# Sporns et al. (2007) / van den Heuvel & Sporns (2013) consensus approach:
# a node is "high" on a centrality if it exceeds mean + sd_mult*SD within the
# network; Hub_score counts how many measures it is high on; Hub flags nodes
# high on >= min_metrics. Centrality_composite is the mean within-network z
# across the measures. Defaults use the four weighted centralities.
add_hub_scores <- function(node_df,
                           metrics = c("Strength", "Weighted_betweenness",
                                       "Weighted_closeness", "Eigenvector"),
                           sd_mult = hub_sd_mult, min_metrics = hub_min_metrics) {
  zmat <- sapply(metrics, function(mc) {
    x <- node_df[[mc]]; mu <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) rep(0, length(x)) else (x - mu) / s
  })
  if (is.null(dim(zmat))) zmat <- matrix(zmat, nrow = nrow(node_df))
  high <- zmat > sd_mult
  node_df$Hub_score            <- rowSums(high, na.rm = TRUE)
  node_df$Centrality_composite <- round(rowMeans(zmat, na.rm = TRUE), 4)
  node_df$Hub                  <- node_df$Hub_score >= min_metrics
  node_df
}

# ================== WEIGHTED NETWORK METRICS (primary, density-matched) ==================
run_network_metrics <- function(region_mat, label,
                                min_module_size = 5, n_consensus = 1000,
                                target_density = primary_density) {
  cor_mat <- cor(t(region_mat)); cor_mat[!is.finite(cor_mat)] <- 0; diag(cor_mat) <- 0
  
  ut_full    <- upper.tri(cor_mat)
  mean_abs_r <- mean(abs(cor_mat[ut_full]), na.rm = TRUE)
  pos_full   <- cor_mat[ut_full]; pos_full <- pos_full[pos_full > 0]
  mean_pos_r <- if (length(pos_full) > 0) mean(pos_full, na.rm = TRUE) else NA_real_
  
  thr_res          <- threshold_to_density(cor_mat, target_density)
  thresh_mat       <- thr_res$adj
  edge_threshold   <- thr_res$threshold
  achieved_density <- thr_res$achieved_density
  
  g <- graph_from_adjacency_matrix(thresh_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
  isolates <- which(degree(g) == 0)
  n_isolated <- length(isolates)
  if (n_isolated > 0) g <- delete_vertices(g, isolates)
  n <- vcount(g); m <- ecount(g)
  
  # Connectivity: after isolate removal a network can still fragment into components.
  # Path-based and eigenvector centralities are only well-defined within a connected
  # component, so we track component structure and compute eigenvector on the giant
  # component (standard in BCT-style analysis), leaving off-giant nodes as NA rather
  # than the artificial ~0 that eigen_centrality() assigns them on disconnected graphs.
  comp        <- igraph::components(g)
  n_components <- comp$no
  giant_frac   <- if (n > 0) max(comp$csize) / n else NA_real_
  in_giant     <- comp$membership == which.max(comp$csize)
  eig_vec      <- setNames(rep(NA_real_, n), V(g)$name)
  if (sum(in_giant) >= 2) {
    g_giant <- induced_subgraph(g, which(in_giant))
    eig_vec[V(g_giant)$name] <- as.numeric(eigen_centrality(g_giant, weights = E(g_giant)$weight)$vector)
  }
  
  max_edges <- n * (n - 1) / 2
  density <- if (max_edges > 0) m / max_edges else 0
  mean_edge_weight <- if (m > 0) mean(E(g)$weight, na.rm = TRUE) else 0
  
  trans_w <- mean(transitivity(g, type = "weighted"), na.rm = TRUE)
  inv_mat_t <- thresh_mat[V(g)$name, V(g)$name]
  inv_mat_t[inv_mat_t > 0] <- 1 / inv_mat_t[inv_mat_t > 0]; diag(inv_mat_t) <- 0
  g_inv <- graph_from_adjacency_matrix(inv_mat_t, mode = "undirected", weighted = TRUE, diag = FALSE)
  avg_path_w <- tryCatch(mean_distance(g_inv, directed = FALSE, weights = E(g_inv)$weight),
                         error = function(e) NA_real_)
  mean_degree <- mean(degree(g), na.rm = TRUE)
  mean_node_strength <- mean(strength(g), na.rm = TRUE)
  
  trans_w <- as.numeric(trans_w)[[1]]; avg_path_w <- as.numeric(avg_path_w)[[1]]
  n_null_graphs <- 100
  if (n > 5 && m > n) {
    rand_C <- numeric(n_null_graphs); rand_L <- numeric(n_null_graphs)
    set.seed(42)
    for (b in seq_len(n_null_graphs)) {
      g_rand <- tryCatch(rewire(g, with = each_edge(prob = 1, loops = FALSE, multiple = FALSE)),
                         error = function(e) NULL)
      if (is.null(g_rand)) { rand_C[b] <- NA; rand_L[b] <- NA; next }
      rand_C[b] <- mean(transitivity(g_rand, type = "weighted"), na.rm = TRUE)
      inv_w_r <- 1 / E(g_rand)$weight; inv_w_r[!is.finite(inv_w_r)] <- 0
      g_rand_inv <- g_rand; E(g_rand_inv)$weight <- inv_w_r
      rand_L[b] <- tryCatch(mean_distance(g_rand_inv, directed = FALSE, weights = E(g_rand_inv)$weight),
                            error = function(e) NA_real_)
    }
    C_rand <- mean(rand_C, na.rm = TRUE); L_rand <- mean(rand_L, na.rm = TRUE)
  } else {
    k_bar  <- as.numeric(mean_degree)[[1]]
    C_rand <- if (isTRUE(n > 1))     as.numeric(k_bar / n)[[1]]           else NA_real_
    L_rand <- if (isTRUE(k_bar > 1)) as.numeric(log(n) / log(k_bar))[[1]] else NA_real_
  }
  sw_sigma <- if (isTRUE(!is.na(avg_path_w)) && isTRUE(!is.na(L_rand)) && isTRUE(L_rand > 0) &&
                  isTRUE(!is.na(trans_w)) && isTRUE(!is.na(C_rand)) && isTRUE(C_rand > 0))
    (trans_w / C_rand) / (avg_path_w / L_rand) else NA_real_
  
  global_df <- data.frame(
    Label = label, N_regions = n, N_edges = m, N_isolated = n_isolated,
    N_components = n_components, Giant_frac = round(giant_frac, 4),
    Target_density = round(target_density, 4), Achieved_density = round(achieved_density, 4),
    Edge_threshold = round(edge_threshold, 4), Network_density = round(density, 4),
    Mean_abs_correlation = round(mean_abs_r, 4),
    Mean_positive_correlation = round(mean_pos_r, 4),
    Mean_edge_weight = round(mean_edge_weight, 4), Mean_degree = round(mean_degree, 2),
    Mean_node_strength = round(mean_node_strength, 4),
    Weighted_transitivity = round(trans_w, 4), Weighted_avg_path = round(avg_path_w, 4),
    SmallWorld_sigma = round(sw_sigma, 4), stringsAsFactors = FALSE
  )
  
  node_df <- data.frame(
    Region = V(g)$name,
    Degree = degree(g), Strength = strength(g),
    Weighted_betweenness = betweenness(g_inv, normalized = TRUE, weights = E(g_inv)$weight),
    Weighted_closeness   = closeness(g_inv, normalized = TRUE, weights = E(g_inv)$weight),
    Eigenvector          = as.numeric(eig_vec[V(g)$name]),
    Local_weighted_transitivity = transitivity(g, type = "weighted"),
    stringsAsFactors = FALSE
  ) %>% arrange(desc(Strength))
  
  net_membership <- consensus_louvain(g, n_runs = n_consensus)
  module_df <- data.frame(Region = names(net_membership),
                          Network_module = as.integer(net_membership),
                          stringsAsFactors = FALSE) %>% arrange(Network_module, Region)
  module_df$Structure <- struct_of(module_df$Region)
  cat("  Modules (consensus Louvain):", length(unique(net_membership)),
      "over", length(net_membership), "regions\n")
  
  roles_df <- compute_node_roles(g, net_membership, min_module_size = min_module_size)
  if (!is.null(roles_df)) {
    node_df <- node_df %>% left_join(roles_df, by = "Region")
    cat("  Roles -", sum(!is.na(node_df$Node_role)), "classified,",
        sum(grepl("hub$", node_df$Node_role), na.rm = TRUE), "G-A hub(s),",
        sum(node_df$Role_reliable == FALSE, na.rm = TRUE), "flagged\n")
  }
  
  # -- Hubs via convergence across centrality measures -----------------------
  # Literature-standard multi-centrality hub identification (Sporns et al. 2007;
  # van den Heuvel & Sporns 2013): a region is flagged "high" on a measure if it
  # exceeds mean + 1 SD within this network; Hub_score counts how many of the
  # four weighted centralities (strength, betweenness, closeness, eigenvector)
  # it is high on; Hub = high on >= hub_min_metrics (default 2). Centrality_
  # composite is the mean within-network z across the four measures (for ranking
  # and cross-condition comparison). Computed at matched density so hub status
  # is comparable across conditions. (The Guimera-Amaral z>=2.5 cartographic hub
  # is reported separately in Node_role; the two answer different questions.)
  node_df <- add_hub_scores(node_df)
  cat("  Centrality hubs:", sum(node_df$Hub, na.rm = TRUE),
      "(>=", hub_min_metrics, "of 4 centralities above mean+1SD)\n")
  
  node_df$Structure <- struct_of(node_df$Region)
  
  list(global = global_df, nodes = node_df, modules = module_df)
}

# ================== DENSITY-SWEEP ROBUSTNESS ==================
run_density_sweep <- function(region_mat, label, densities = sweep_densities) {
  set.seed(42)            # cluster_louvain() below is stochastic
  cor_mat <- cor(t(region_mat)); cor_mat[!is.finite(cor_mat)] <- 0; diag(cor_mat) <- 0
  rows <- list()
  for (d in densities) {
    tr <- threshold_to_density(cor_mat, d)
    g  <- graph_from_adjacency_matrix(tr$adj, mode = "undirected", weighted = TRUE, diag = FALSE)
    iso <- which(degree(g) == 0); if (length(iso) > 0) g <- delete_vertices(g, iso)
    n <- vcount(g); m <- ecount(g)
    if (n < 3 || m < 1) {
      rows[[length(rows) + 1]] <- data.frame(
        Label = label, Target_density = d, Achieved_density = round(tr$achieved_density, 4),
        N_regions = n, N_edges = m, Mean_degree = NA_real_, Mean_strength = NA_real_,
        Weighted_transitivity = NA_real_, Global_efficiency = NA_real_, Modularity_Q = NA_real_)
      next
    }
    W <- as.matrix(igraph::as_adjacency_matrix(g, attr = "weight"))
    invW <- W; invW[invW > 0] <- 1 / invW[invW > 0]
    g_inv <- graph_from_adjacency_matrix(invW, mode = "undirected", weighted = TRUE, diag = FALSE)
    dmat <- distances(g_inv, weights = E(g_inv)$weight)
    offd <- dmat[upper.tri(dmat)]
    # Unreachable pairs (different components) contribute 0, per the standard
    # definition; averaging only over reachable pairs would score a FRAGMENTED
    # network higher than a connected one of the same density.
    eff  <- mean(ifelse(is.finite(offd) & offd > 0, 1 / offd, 0))
    trans_w <- mean(transitivity(g, type = "weighted"), na.rm = TRUE)
    comm <- cluster_louvain(g, weights = E(g)$weight)
    Qd   <- modularity(g, membership(comm), weights = E(g)$weight)
    rows[[length(rows) + 1]] <- data.frame(
      Label = label, Target_density = d, Achieved_density = round(tr$achieved_density, 4),
      N_regions = n, N_edges = m,
      Mean_degree = round(mean(degree(g)), 3), Mean_strength = round(mean(strength(g)), 4),
      Weighted_transitivity = round(trans_w, 4), Global_efficiency = round(eff, 4),
      Modularity_Q = round(Qd, 4))
  }
  sweep_df <- do.call(rbind, rows)
  auc <- function(x, y) {
    ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
    if (length(x) < 2) return(NA_real_)
    o <- order(x); x <- x[o]; y <- y[o]; rng <- max(x) - min(x)
    if (rng <= 0) return(mean(y))
    sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2) / rng
  }
  xs <- sweep_df$Achieved_density
  auc_df <- data.frame(
    Label = label,
    AUC_mean_degree           = round(auc(xs, sweep_df$Mean_degree), 4),
    AUC_mean_strength         = round(auc(xs, sweep_df$Mean_strength), 4),
    AUC_weighted_transitivity = round(auc(xs, sweep_df$Weighted_transitivity), 4),
    AUC_global_efficiency     = round(auc(xs, sweep_df$Global_efficiency), 4),
    AUC_modularity_Q          = round(auc(xs, sweep_df$Modularity_Q), 4),
    stringsAsFactors = FALSE
  )
  list(sweep = sweep_df, auc = auc_df)
}

# ================== STATISTICAL INFERENCE ON NETWORK METRICS ==================
# Each condition yields a single correlation network, so a metric is one point
# estimate with no error bars. To make inferential claims about metric CHANGES
# between conditions, we resample over animals (the unit of replication):
#   (1) bootstrap CIs per condition  -> metric +/- uncertainty;
#   (2) label-permutation test        -> p-value for between-condition differences.
# This is the standard resampling approach for cFos/correlation-network metric
# comparison (e.g. Vetere et al. 2017; Kimbrough et al. 2020). NOTE: with small
# n, bootstrap-with-replacement can duplicate animals (mildly inflating r) and
# CIs will be wide, which appropriately reflects the small sample size.

# Core: build the density-matched network from a region x animal matrix and
# return the comparable global metrics as a named vector.
compute_global_metric_vector <- function(region_mat, target_density = primary_density) {
  cm <- cor(t(region_mat)); cm[!is.finite(cm)] <- 0; diag(cm) <- 0
  g <- graph_from_adjacency_matrix(threshold_to_density(cm, target_density)$adj,
                                   mode = "undirected", weighted = TRUE, diag = FALSE)
  iso <- which(degree(g) == 0); if (length(iso) > 0) g <- delete_vertices(g, iso)
  if (ecount(g) < 1) return(c(Mean_strength = NA_real_, Weighted_transitivity = NA_real_,
                              Global_efficiency = NA_real_, Modularity_Q = NA_real_))
  W <- as.matrix(igraph::as_adjacency_matrix(g, attr = "weight"))
  invW <- W; invW[invW > 0] <- 1 / invW[invW > 0]
  g_inv <- graph_from_adjacency_matrix(invW, mode = "undirected", weighted = TRUE, diag = FALSE)
  dmat <- distances(g_inv, weights = E(g_inv)$weight); offd <- dmat[upper.tri(dmat)]
  eff  <- mean(ifelse(is.finite(offd) & offd > 0, 1 / offd, 0))   # unreachable pairs = 0
  comm <- cluster_louvain(g, weights = E(g)$weight)
  c(Mean_strength         = mean(strength(g)),
    Weighted_transitivity = mean(transitivity(g, type = "weighted"), na.rm = TRUE),
    Global_efficiency     = eff,
    Modularity_Q          = modularity(g, membership(comm), weights = E(g)$weight))
}

# (1) Bootstrap CIs per condition (resample animals with replacement)
bootstrap_network_metrics <- function(region_mat, label, n_boot = 1000,
                                      target_density = primary_density, seed = 42) {
  # Seed BEFORE the observed statistic: compute_global_metric_vector() calls
  # cluster_louvain(), which is stochastic, so the point estimate would otherwise
  # depend on inherited RNG state and move when n_consensus / n_perm are changed.
  set.seed(seed)
  obs   <- compute_global_metric_vector(region_mat, target_density)
  nanim <- ncol(region_mat)
  boot <- matrix(NA_real_, n_boot, length(obs), dimnames = list(NULL, names(obs)))
  for (b in seq_len(n_boot)) {
    idx <- sample(nanim, nanim, replace = TRUE)
    boot[b, ] <- compute_global_metric_vector(region_mat[, idx, drop = FALSE], target_density)
  }
  data.frame(
    Label = label, Metric = names(obs), Observed = round(obs, 4),
    Boot_mean = round(colMeans(boot, na.rm = TRUE), 4),
    Boot_SD   = round(apply(boot, 2, sd, na.rm = TRUE), 4),
    CI_lower  = round(apply(boot, 2, quantile, 0.025, na.rm = TRUE), 4),
    CI_upper  = round(apply(boot, 2, quantile, 0.975, na.rm = TRUE), 4),
    n_boot = n_boot, stringsAsFactors = FALSE, row.names = NULL
  )
}

# (2) Permutation test of metric difference between two conditions.
# Pools animals, shuffles condition labels (preserving group sizes), recomputes
# the metric difference under the null; two-sided p = P(|null diff| >= |obs diff|).
permutation_metric_diff <- function(rm_A, rm_B, label_A, label_B, n_perm = 1000,
                                    target_density = primary_density, seed = 42) {
  shared <- intersect(rownames(rm_A), rownames(rm_B))
  if (length(shared) < 10) return(NULL)
  A <- rm_A[shared, , drop = FALSE]; B <- rm_B[shared, , drop = FALSE]
  nA <- ncol(A); nB <- ncol(B); pooled <- cbind(A, B)
  set.seed(seed)          # seed before the observed difference (Louvain is stochastic)
  obs <- compute_global_metric_vector(A, target_density) -
    compute_global_metric_vector(B, target_density)
  null <- matrix(NA_real_, n_perm, length(obs), dimnames = list(NULL, names(obs)))
  for (p in seq_len(n_perm)) {
    perm <- sample(nA + nB)
    ga <- pooled[, perm[1:nA], drop = FALSE]
    gb <- pooled[, perm[(nA + 1):(nA + nB)], drop = FALSE]
    null[p, ] <- compute_global_metric_vector(ga, target_density) -
      compute_global_metric_vector(gb, target_density)
  }
  # Add-one permutation p (see run_permutation_modularity): (1 + k) / (1 + n_valid),
  # so a p of exactly 0 is never reported.
  pvals <- sapply(seq_along(obs), function(j) {
    nd <- null[, j]; nd <- nd[is.finite(nd)]
    if (length(nd) == 0 || !is.finite(obs[j])) return(NA_real_)
    (1 + sum(abs(nd) >= abs(obs[j]))) / (1 + length(nd))
  })
  data.frame(
    Comparison = paste0(label_A, " vs ", label_B), Metric = names(obs),
    Diff_observed = round(obs, 4), p_perm = round(pvals, 4),
    n_perm = n_perm, stringsAsFactors = FALSE, row.names = NULL
  )
}

# ================== MODULARITY SIGNIFICANCE (degree-preserving null) ==================
# Tests whether observed weighted modularity exceeds what is expected GIVEN THE
# NETWORK'S DEGREE SEQUENCE. The null is a configuration-model ensemble: the
# observed density-matched graph is randomly rewired preserving each node's
# degree (Maslov-Sneppen edge swaps; Rubinov & Sporns 2010), which is the
# standard modularity-significance null. Note that this differs from the null
# used for the small-world index above, which uses each-edge rewiring and does
# not preserve degree. This replaces the earlier animal-label-shuffle null, which
# compared against near-random graphs whose Louvain modularity is inflated and
# produced uninterpretable results. Single-run Louvain is used for observed and
# every null (same procedure); the role partition separately uses consensus
# Louvain. Reports observed Q, null mean/SD, a modularity z-score, and p.
run_permutation_modularity <- function(region_mat, label, n_perm = 1000,
                                       target_density = primary_density) {
  cat("  Modularity vs degree-preserving null:", label, "\n")
  cm <- cor(t(region_mat)); cm[!is.finite(cm)] <- 0; diag(cm) <- 0
  a  <- threshold_to_density(cm, target_density)$adj
  g_obs <- graph_from_adjacency_matrix(a, mode = "undirected", weighted = TRUE, diag = FALSE)
  iso <- which(degree(g_obs) == 0); if (length(iso) > 0) g_obs <- delete_vertices(g_obs, iso)
  
  Q_louvain <- function(g) {
    if (ecount(g) < 1) return(NA_real_)
    cl <- cluster_louvain(g, weights = E(g)$weight)
    modularity(g, membership(cl), weights = E(g)$weight)
  }
  set.seed(42)            # seed before the observed Q (Louvain is stochastic)
  Q_obs <- Q_louvain(g_obs)
  
  Q_null <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    g_rand <- tryCatch(
      rewire(g_obs, with = keeping_degseq(loops = FALSE, niter = ecount(g_obs) * 10)),
      error = function(e) NULL)
    Q_null[p] <- if (is.null(g_rand)) NA_real_ else Q_louvain(g_rand)
  }
  mu <- mean(Q_null, na.rm = TRUE); sdn <- sd(Q_null, na.rm = TRUE)
  # Add-one (Davison & Hinkley) permutation p: (1 + #{null >= obs}) / (1 + n_valid).
  # A permutation test cannot produce p = 0 -- with n_perm resamples the smallest
  # attainable value is 1/(n_perm + 1), and reporting a literal 0 in a deposited
  # workbook is not defensible. p_upper_bound records that floor so the value can
  # be quoted as "p < x" rather than as an exact zero.
  n_valid <- sum(is.finite(Q_null))
  p_perm  <- (1 + sum(Q_null >= Q_obs, na.rm = TRUE)) / (1 + n_valid)
  data.frame(
    Label = label, Observed_Q = round(Q_obs, 4),
    Null_Mean_Q = round(mu, 4), Null_SD_Q = round(sdn, 4),
    Modularity_z = round(if (is.finite(sdn) && sdn > 0) (Q_obs - mu) / sdn else NA_real_, 3),
    p_value_permutation = signif(p_perm, 4),
    p_upper_bound = signif(1 / (1 + n_valid), 3),
    n_permutations = n_perm, n_valid_nulls = n_valid,
    null_model = "degree_preserving",
    stringsAsFactors = FALSE
  )
}

# ================== FISHER Z CONNECTIVITY DIFFERENCES ==================
# Tests whether each region-pair correlation differs between two conditions
# (Fisher r-to-z, two-sided, BH-FDR across all pairs). Per-pair r difference
# (Delta_r = r_A - r_B) folded in. Returns FDR-significant pairs only.
run_fisher_z_test <- function(cor_A, n_A, cor_B, n_B, label_A, label_B, sex_label) {
  # Combined matrices are mean-centered within sex before pooling, which removes one
  # further degree of freedom per correlation, so the Fisher variance uses n - 4
  # there and the usual n - 3 for the single-sex networks.
  df_off <- if (identical(sex_label, "Combined")) 4 else 3
  shared <- intersect(rownames(cor_A), rownames(cor_B))
  if (length(shared) < 3) return(NULL)
  cA <- cor_A[shared, shared]; cB <- cor_B[shared, shared]
  pairs <- which(upper.tri(matrix(0, length(shared), length(shared))), arr.ind = TRUE)
  
  z_vals <- numeric(nrow(pairs)); p_vals <- numeric(nrow(pairs))
  rA_vals <- numeric(nrow(pairs)); rB_vals <- numeric(nrow(pairs))
  for (k in seq_len(nrow(pairs))) {
    i <- pairs[k, 1]; j <- pairs[k, 2]
    r1 <- cA[i, j]; r2 <- cB[i, j]
    rA_vals[k] <- r1; rB_vals[k] <- r2
    r1 <- max(min(r1, 0.9999), -0.9999); r2 <- max(min(r2, 0.9999), -0.9999)
    se <- sqrt(1 / (n_A - df_off) + 1 / (n_B - df_off))
    z_stat <- (atanh(r1) - atanh(r2)) / se
    z_vals[k] <- z_stat; p_vals[k] <- 2 * pnorm(-abs(z_stat))
  }
  p_adj <- p.adjust(p_vals, method = "BH")
  sig <- which(p_adj < 0.05)
  if (length(sig) == 0) return(NULL)
  data.frame(
    Comparison  = paste0(label_A, " vs ", label_B),
    Sex         = sex_label,
    Region_A    = shared[pairs[sig, 1]],
    Region_B    = shared[pairs[sig, 2]],
    r_A         = round(rA_vals[sig], 4),
    r_B         = round(rB_vals[sig], 4),
    Delta_r     = round(rA_vals[sig] - rB_vals[sig], 4),
    Z_statistic = round(z_vals[sig], 4),
    p_value     = signif(p_vals[sig], 4),
    p_adj_BH    = signif(p_adj[sig], 4),
    stringsAsFactors = FALSE
  ) %>% arrange(p_adj_BH)
}

# ================== HELPERS: SHEET NAMING ==================
trt_abbr <- c("Vehicle"="Veh", "Acute Morphine"="AcMor",
              "Morphine-Dependent"="MorDep", "Ro 64-6198"="Ro")
sex_abbr <- c("Males"="M", "Females"="F", "Combined"="C")
safe_sheet <- function(s) substr(gsub("[\\[\\]:*?/\\\\]", "_", s), 1, 31)

# =============================================================================
# MAIN
# =============================================================================
cat("Loading data from:", file_path, "\n\n")
males_by_trt   <- load_sheet(file_path, "Clustered Males")
females_by_trt <- load_sheet(file_path, "Clustered Females")

# Per-condition sample sizes (for Fisher Z)
n_samples <- list()
for (trt in treatments) {
  n_m <- if (trt %in% names(males_by_trt))   ncol(males_by_trt[[trt]])   - 1 else 0
  n_f <- if (trt %in% names(females_by_trt)) ncol(females_by_trt[[trt]]) - 1 else 0
  n_samples[[trt]] <- list(Males = n_m, Females = n_f, Combined = n_m + n_f)
}

# ---- Build cleaned matrices + correlation matrices per treatment x sex ----
region_mats  <- list()
cor_matrices <- list()
for (trt in treatments) {
  region_mats[[trt]] <- list(); cor_matrices[[trt]] <- list()
  has_m <- trt %in% names(males_by_trt); has_f <- trt %in% names(females_by_trt)
  
  if (has_m) {
    rm_ <- clean_and_matrix(males_by_trt[[trt]])
    region_mats[[trt]][["Males"]]  <- rm_
    cm <- cor(t(rm_)); cm[!is.finite(cm)] <- 0; cor_matrices[[trt]][["Males"]] <- cm
  }
  if (has_f) {
    rf_ <- clean_and_matrix(females_by_trt[[trt]])
    region_mats[[trt]][["Females"]] <- rf_
    cm <- cor(t(rf_)); cm[!is.finite(cm)] <- 0; cor_matrices[[trt]][["Females"]] <- cm
  }
  if (has_m && has_f) {
    # Pool sexes for the Combined network. Reuse the SAME per-sex matrices used as the
    # primary analysis, then mean-center each region WITHIN sex before pooling. This
    # residualizes sex out per region (removes the between-sex mean shift that would
    # otherwise create spurious cross-region correlations when M and F are pooled),
    # while leaving the within-sex covariance structure the network relies on intact.
    m_mat <- region_mats[[trt]][["Males"]]
    f_mat <- region_mats[[trt]][["Females"]]
    m_mat <- m_mat - rowMeans(m_mat, na.rm = TRUE)
    f_mat <- f_mat - rowMeans(f_mat, na.rm = TRUE)
    shared <- intersect(rownames(m_mat), rownames(f_mat))
    if (length(shared) >= 2) {
      combined <- cbind(m_mat[shared, , drop = FALSE], f_mat[shared, , drop = FALSE])
      rownames(combined) <- shared
      region_mats[[trt]][["Combined"]] <- combined
      cm <- cor(t(combined)); cm[!is.finite(cm)] <- 0; cor_matrices[[trt]][["Combined"]] <- cm
    }
  }
}

# ---- Network metrics, roles, modules, sweep, permutation per network ----
global_list <- list(); nodes_list <- list(); modules_list <- list()
sweep_list  <- list(); auc_list   <- list(); perm_list <- list()
net_modules <- list()  # for cross-treatment stability

for (trt in treatments) {
  net_modules[[trt]] <- list()
  for (sex in sex_groups) {
    mat <- region_mats[[trt]][[sex]]
    if (is.null(mat)) next
    label <- paste(trt, sex)
    cat("Processing:", label, "\n")
    
    res <- run_network_metrics(mat, label, min_module_size = min_module_size,
                               n_consensus = n_consensus)
    res$global$Treatment <- trt;  res$global$Sex <- sex
    res$nodes$Treatment  <- trt;  res$nodes$Sex  <- sex
    res$modules$Treatment <- trt; res$modules$Sex <- sex
    global_list[[label]]  <- res$global
    nodes_list[[label]]   <- res$nodes
    modules_list[[label]] <- res$modules
    net_modules[[trt]][[sex]] <- res$modules
    
    sw <- run_density_sweep(mat, label)
    sw$sweep$Treatment <- trt; sw$sweep$Sex <- sex
    sw$auc$Treatment   <- trt; sw$auc$Sex   <- sex
    sweep_list[[label]] <- sw$sweep; auc_list[[label]] <- sw$auc
    
    pm <- run_permutation_modularity(mat, label, n_perm = n_perm)
    pm$Treatment <- trt; pm$Sex <- sex
    perm_list[[label]] <- pm
    cat("\n")
  }
}

# ---- Module-detection method comparison (lineage hierarchical clustering vs Louvain) ----
# Documents why consensus Louvain is used: reports the module count the standard cFos
# hierarchical-clustering method (Kimbrough 2020; Ardinger 2024) returns across dendrogram
# cut heights, next to the consensus-Louvain count, for every network.
cat("Module-detection comparison (hierarchical clustering vs consensus Louvain)\n")
modmethod_list <- list()
hc_cut_fracs <- c(0.3, 0.4, 0.5, 0.6, 0.7)
for (trt in treatments) for (sex in sex_groups) {
  cm <- cor_matrices[[trt]][[sex]]; if (is.null(cm)) next
  md <- modules_list[[paste(trt, sex)]]
  louvain_k <- if (!is.null(md)) length(unique(md$Network_module)) else NA_integer_
  hc_counts <- hclust_module_count(cm, hc_cut_fracs)
  modmethod_list[[paste(trt, sex)]] <- data.frame(
    c(list(Treatment = trt, Sex = sex, Louvain_modules = louvain_k),
      as.list(hc_counts)),
    check.names = FALSE, stringsAsFactors = FALSE)
}

# ---- Fisher Z connectivity differences across condition pairs ----
fisher_list <- list()
for (comp in comparisons) {
  for (sex in sex_groups) {
    cA <- cor_matrices[[comp[1]]][[sex]]; cB <- cor_matrices[[comp[2]]][[sex]]
    nA <- n_samples[[comp[1]]][[sex]];    nB <- n_samples[[comp[2]]][[sex]]
    if (!is.null(cA) && !is.null(cB) && nA >= 4 && nB >= 4) {
      fz <- run_fisher_z_test(cA, nA, cB, nB, comp[1], comp[2], sex)
      if (!is.null(fz)) fisher_list[[paste(comp[1], comp[2], sex)]] <- fz
    }
  }
}

# ---- Inferential stats on network metrics (PRIMARY: Combined networks) ----
# Bootstrap CIs per condition + label-permutation tests of between-condition
# differences. These turn the per-condition point estimates into testable claims.
cat("\n========== BOOTSTRAP + PERMUTATION ON NETWORK METRICS (Combined) ==========\n")
boot_list <- list()
for (trt in treatments) {
  mat <- region_mats[[trt]][["Combined"]]
  if (is.null(mat)) next
  cat("  Bootstrapping:", trt, "Combined\n")
  bdf <- bootstrap_network_metrics(mat, paste(trt, "Combined"), n_boot = n_boot)
  bdf$Treatment <- trt; bdf$Sex <- "Combined"
  boot_list[[trt]] <- bdf
}
permdiff_list <- list()
for (comp in comparisons) {
  mA <- region_mats[[comp[1]]][["Combined"]]; mB <- region_mats[[comp[2]]][["Combined"]]
  if (is.null(mA) || is.null(mB)) next
  cat("  Permutation diff:", comp[1], "vs", comp[2], "Combined\n")
  pm_diff <- permutation_metric_diff(mA, mB, comp[1], comp[2], n_perm = n_boot)
  if (!is.null(pm_diff)) permdiff_list[[paste(comp[1], comp[2])]] <- pm_diff
}
metric_comparisons <- NULL
if (length(permdiff_list) > 0) {
  metric_comparisons <- bind_rows(permdiff_list)
  # Mean_strength at matched density mostly re-expresses coactivation STRENGTH (it tracks the
  # weight of the retained top-density edges), so it is tested in its own family, separate from
  # the graph-TOPOLOGY metrics (segregation/integration/community structure). BH-FDR is applied
  # WITHIN each family. This mirrors the cFos lineage, which analyses correlation strength
  # (mean R) separately from graph-theoretic topology (Kimbrough 2020; Bloch 2024) rather than
  # pooling a strength readout into the topology multiple-comparison correction.
  metric_comparisons <- metric_comparisons %>%
    mutate(Metric_family = ifelse(Metric == "Mean_strength", "Strength", "Topology")) %>%
    group_by(Metric_family) %>%
    mutate(p_adj_BH = round(p.adjust(p_perm, method = "BH"), 4)) %>%
    ungroup() %>%
    dplyr::select(Comparison, Metric, Metric_family, Diff_observed, p_perm, p_adj_BH, n_perm)
  cat("  Significant TOPOLOGY differences (p_adj_BH < 0.05):",
      sum(metric_comparisons$p_adj_BH < 0.05 & metric_comparisons$Metric_family == "Topology", na.rm = TRUE),
      "of", sum(metric_comparisons$Metric_family == "Topology"),
      "| Strength:",
      sum(metric_comparisons$p_adj_BH < 0.05 & metric_comparisons$Metric_family == "Strength", na.rm = TRUE),
      "of", sum(metric_comparisons$Metric_family == "Strength"), "\n")
}

# ---- Cross-treatment module stability (on consensus network modules) ----
# Per region, proportion of conditions in which it lands in its MODAL module.
#
# Louvain community numbers are arbitrary WITHIN each network: "module 1" in
# Vehicle and "module 1" in Ro 64-6198 are different sets of regions (on these
# data they overlap by only ~25-30%). Counting how often a region keeps the same
# NUMBER across conditions therefore measures nothing. Each condition's partition
# is first relabelled onto a reference condition's by greedy maximum overlap:
# repeatedly take the (this-condition, reference) module pair sharing the most
# regions, bind them one-to-one and remove both. Modules left unmatched keep a
# distinct label of their own so they can never collide with a reference label.
match_partition <- function(x, ref) {
  lv <- sort(unique(x[!is.na(x)]))
  if (!length(lv)) return(x)
  map <- setNames(rep(NA_integer_, length(lv)), as.character(lv))
  ok  <- !is.na(x) & !is.na(ref)
  if (any(ok)) {
    tab <- table(factor(x[ok], levels = lv), ref[ok])
    while (nrow(tab) > 0 && ncol(tab) > 0 && max(tab) > 0) {
      k <- which(tab == max(tab), arr.ind = TRUE)[1, ]
      map[rownames(tab)[k[1]]] <- as.integer(colnames(tab)[k[2]])
      tab <- tab[-k[1], -k[2], drop = FALSE]
    }
  }
  n_un <- sum(is.na(map))
  if (n_un > 0) map[is.na(map)] <- 1000L + seq_len(n_un)
  unname(map[as.character(x)])
}

stability_list <- list()
for (sex in sex_groups) {
  region_set <- character(0)
  for (trt in treatments) {
    md <- net_modules[[trt]][[sex]]
    if (!is.null(md)) region_set <- union(region_set, md$Region)
  }
  if (length(region_set) == 0) next
  assign_mat <- matrix(NA_integer_, length(region_set), length(treatments),
                       dimnames = list(region_set, treatments))
  for (trt in treatments) {
    md <- net_modules[[trt]][[sex]]
    if (is.null(md)) next
    assign_mat[md$Region, trt] <- md$Network_module
  }
  # Relabel every condition onto the reference condition before comparing.
  has_mod <- vapply(treatments, function(t) !is.null(net_modules[[t]][[sex]]), logical(1))
  ref_trt <- treatments[has_mod][1]
  for (trt in setdiff(treatments[has_mod], ref_trt))
    assign_mat[, trt] <- match_partition(assign_mat[, trt], assign_mat[, ref_trt])

  stab <- apply(assign_mat, 1, function(x) {
    x <- x[!is.na(x)]; if (length(x) == 0) return(NA_real_)
    max(table(x)) / length(x)
  })
  n_assigned <- apply(assign_mat, 1, function(x) sum(!is.na(x)))
  stability_list[[sex]] <- data.frame(
    Sex = sex, Region = names(stab), Structure = struct_of(names(stab)),
    Modal_module_proportion = round(stab, 3), N_treatments_assigned = n_assigned,
    Matched_to = ref_trt,
    stringsAsFactors = FALSE
  ) %>% arrange(desc(Modal_module_proportion), Region)
}

# ---- Hub changes across conditions (pooled Combined network) ----
# Per region, multi-centrality hub status in each condition; identifies regions
# whose hub status shifts relative to Vehicle (baseline). Run on the pooled Combined
# network as a power-favorable companion; the primary sex-stratified hub flags live
# per-sex in Node_Roles.
hub_changes <- NULL
comb_labels <- names(nodes_list)[sapply(nodes_list, function(d) d$Sex[1] == "Combined")]
if (length(comb_labels) > 0) {
  comb_nodes <- bind_rows(nodes_list[comb_labels])
  trt_present <- intersect(treatments, unique(comb_nodes$Treatment))
  score_wide <- comb_nodes %>% dplyr::select(Treatment, Region, Hub_score) %>%
    pivot_wider(names_from = Treatment, values_from = Hub_score)
  hub_wide <- comb_nodes %>% dplyr::select(Treatment, Region, Hub) %>%
    pivot_wider(names_from = Treatment, values_from = Hub)
  hub_mat <- as.matrix(hub_wide[, trt_present, drop = FALSE]); hub_mat[is.na(hub_mat)] <- FALSE
  storage.mode(hub_mat) <- "logical"
  
  hc <- score_wide %>%
    left_join(comb_nodes %>% distinct(Region, Structure), by = "Region")
  hc$Hub_in_n_conditions <- rowSums(hub_mat)
  if ("Vehicle" %in% trt_present) {
    other <- setdiff(trt_present, "Vehicle")
    base_hub <- hub_mat[, "Vehicle"]
    hc$Change_vs_Vehicle <- vapply(seq_len(nrow(hub_mat)), function(i) {
      bh <- base_hub[i]; oh <- hub_mat[i, other]
      gained <- other[oh & !bh]; lost <- if (bh) other[!oh] else character(0)
      parts <- c(if (length(gained)) paste0("gained: ", paste(gained, collapse = "/")),
                 if (length(lost))   paste0("lost: ",   paste(lost,   collapse = "/")))
      if (length(parts) == 0) (if (bh) "stable hub" else "never hub") else paste(parts, collapse = "; ")
    }, character(1))
  }
  hub_changes <- hc %>% filter(Hub_in_n_conditions >= 1) %>%
    arrange(desc(Hub_in_n_conditions), Region) %>%
    dplyr::select(Region, Structure, all_of(trt_present),
                  Hub_in_n_conditions, any_of("Change_vs_Vehicle"))
  cat("\nHub-change analysis (Combined):", nrow(hub_changes),
      "regions are a hub in >=1 condition\n")
}

# =============================================================================
# WRITE TOPIC-GROUPED EXCEL WORKBOOKS
# =============================================================================
hs <- createStyle(textDecoration = "bold", border = "bottom")
add_sheet <- function(wb, name, df, rowNames = FALSE) {
  nm <- safe_sheet(name)
  addWorksheet(wb, nm)
  writeData(wb, nm, df, rowNames = rowNames)
  addStyle(wb, nm, hs, rows = 1, cols = 1:(ncol(df) + as.integer(rowNames)), gridExpand = TRUE)
  setColWidths(wb, nm, cols = 1:(ncol(df) + as.integer(rowNames)), widths = "auto")
}
readme_sheet <- function(wb, fields, descs) {
  addWorksheet(wb, "README")
  writeData(wb, "README", data.frame(Field = fields, Description = descs))
  addStyle(wb, "README", hs, rows = 1, cols = 1:2)
  setColWidths(wb, "README", cols = 1:2, widths = c(34, 110))
}

# ---- Workbook 0: Processed input matrices (raw counts -> log10 region x animal) ----
wb0 <- createWorkbook()
readme_sheet(wb0,
             c("Purpose", "Processing", "Sheets", "Combined networks"),
             c("Cleaned region x animal matrices that feed the correlation networks (the raw-counts -> input step).",
               "Per region: negatives -> NA; zeros -> half of the minimum non-zero value; log10 transform.",
               "One sheet per condition x sex; rows = regions, columns = animals.",
               "Combined sheets mean-center each region within sex before pooling males + females (residualizes the sex mean-difference)."))
for (trt in treatments) for (sex in sex_groups) {
  mat <- region_mats[[trt]][[sex]]; if (is.null(mat)) next
  add_sheet(wb0, paste0("Data ", trt_abbr[[trt]], " ", sex_abbr[[sex]]),
            as.data.frame(mat), rowNames = TRUE)
}
saveWorkbook(wb0, file.path(output_root, "00_Processed_Data.xlsx"), overwrite = TRUE)

# ---- Workbook 1: Connectivity (correlation matrices + Fisher Z differences) ----
wb1 <- createWorkbook()
readme_sheet(wb1,
             c("Correlation matrices", "Method", "FisherZ_SigPairs", "Delta_r", "Significance"),
             c("Region x region Pearson correlation across animals (co-activation), one sheet per condition x sex.",
               "Pearson r on log10 (mean-centered within sex for Combined) region profiles; this is the un-thresholded basis for all network construction.",
               "Region pairs whose correlation differs between two conditions (Fisher r-to-z, two-sided).",
               "Per-pair correlation difference r_A - r_B (difference matrices folded in here).",
               "BH-FDR across all region pairs within each comparison; only FDR-significant pairs (p_adj_BH < 0.05) listed."))
for (trt in treatments) for (sex in sex_groups) {
  cm <- cor_matrices[[trt]][[sex]]; if (is.null(cm)) next
  add_sheet(wb1, paste0("Corr ", trt_abbr[[trt]], " ", sex_abbr[[sex]]),
            as.data.frame(round(cm, 4)), rowNames = TRUE)
}
if (length(fisher_list) > 0) add_sheet(wb1, "FisherZ_SigPairs", bind_rows(fisher_list))
saveWorkbook(wb1, file.path(output_root, "01_Connectivity.xlsx"), overwrite = TRUE)

# ---- Workbook 2: Network metrics (global, density sweep, AUC, permutation) ----
wb2 <- createWorkbook()
readme_sheet(wb2,
             c("Thresholding", "Global_Metrics", "Mean_positive_correlation",
               "Density_Sweep", "Density_Sweep_AUC", "Permutation_Modularity", "SmallWorld_sigma",
               "Bootstrap_CI", "Metric_Comparisons", "Diff_observed / p_perm / p_adj_BH"),
             c("Density-matched (proportional) thresholding on positive edges to a common density (primary_density) across all conditions; positive-only (Brynildsen 2020; Xie 2024); Rubinov & Sporns 2010; van den Heuvel 2017.",
               "One row per condition x sex: whole-network summary on the density-matched graph. Thresholding fixes the EDGE COUNT (density x 198*197/2) identically for every network; because the number of degree-0 isolates then removed differs by condition, Network_density (edges / possible edges among the RETAINED nodes) is slightly higher than Target_density and varies across conditions - report both. N_components / Giant_frac report fragmentation after isolate removal; eigenvector centrality is computed on the giant component (off-giant nodes = NA).",
               "Mean positive Pearson r per condition (overall co-activation strength) - lets the van den Heuvel (2017) proportional-thresholding caveat be assessed.",
               "Global metrics recomputed across densities (2-20%), each condition at matched density: degree, strength, weighted clustering, global efficiency, modularity Q (Garrison 2015; Hallquist & Hillary 2019).",
               "Normalized area under the density-metric curve (= mean metric over the swept range): one threshold-independent value per metric per condition.",
               "Observed weighted modularity vs a DEGREE-PRESERVING (configuration-model) null: the graph is rewired keeping each node's degree (Maslov-Sneppen swaps; Rubinov & Sporns 2010), the standard modularity-significance null. Reports observed Q, null mean/SD, a modularity z-score, and p. Single-run Louvain for observed and all nulls.",
               "Small-world index (C_obs/C_rand)/(L_obs/L_rand) using 100 randomly rewired nulls (igraph each_edge rewiring, which does NOT preserve the degree sequence; Humphries & Gurney 2008); >1 = small-world. Note this differs from the degree-preserving null used for modularity significance.",
               "PRIMARY inferential layer (Combined networks). Bootstrap CIs: animals resampled with replacement, network rebuilt at matched density, metric recomputed (n_boot iterations). Gives metric +/- 95% CI per condition (Vetere 2017; Kimbrough 2020).",
               "Label-permutation tests of between-condition metric differences (Combined): animals pooled, condition labels shuffled, difference recomputed under the null; two-sided p. THIS is the test of whether a network metric changes significantly between conditions.",
               "Diff_observed = metric(A) - metric(B); p_perm = permutation p-value; Metric_family splits Mean_strength (Strength channel, which re-expresses coactivation strength at matched density) from the graph-topology metrics (weighted transitivity, global efficiency, modularity Q). p_adj_BH = BH-FDR applied WITHIN each family, mirroring the lineage's separation of correlation strength (mean R) from graph topology (Kimbrough 2020; Bloch 2024)."))
if (length(global_list) > 0) add_sheet(wb2, "Global_Metrics",
                                       bind_rows(global_list) %>% dplyr::select(Treatment, Sex, everything()))
if (length(sweep_list) > 0)  add_sheet(wb2, "Density_Sweep",
                                       bind_rows(sweep_list) %>% dplyr::select(Treatment, Sex, everything()))
if (length(auc_list) > 0)    add_sheet(wb2, "Density_Sweep_AUC",
                                       bind_rows(auc_list) %>% dplyr::select(Treatment, Sex, everything()))
if (length(perm_list) > 0)   add_sheet(wb2, "Permutation_Modularity",
                                       bind_rows(perm_list) %>% dplyr::select(Treatment, Sex, everything()))
if (length(boot_list) > 0)   add_sheet(wb2, "Bootstrap_CI",
                                       bind_rows(boot_list) %>% dplyr::select(Treatment, Sex, Metric, Observed,
                                                                              Boot_mean, Boot_SD, CI_lower, CI_upper, n_boot))
if (!is.null(metric_comparisons)) add_sheet(wb2, "Metric_Comparisons", metric_comparisons)
saveWorkbook(wb2, file.path(output_root, "02_Network_Metrics.xlsx"), overwrite = TRUE)

# ---- Workbook 3: Modules, node roles (WMDz + PC), and hubs ----
wb3 <- createWorkbook()
readme_sheet(wb3,
             c("Network modules", "Node_Roles", "WM_degree_z", "Participation_coef",
               "WM_strength_z / Participation_coef_wt", "Node_role (G-A cartographic)",
               "Centrality hubs", "Hub_score / Hub", "Centrality_composite",
               "Hub_Regions", "Hub_Changes_Combined", "Role_reliable",
               "Module_Composition", "Module_Stability", "Module_Method_Comparison", "Primary networks"),
             c("Consensus Louvain communities detected on the thresholded graph itself (1000 runs; Lancichinetti & Fortunato 2012) - the graph-native partition for the cartographic metrics.",
               "Per-region node-level metrics for every condition x sex: centralities, the two cartographic metrics, the G-A role, and the multi-centrality hub flags.",
               "Within-module degree z-score (binary form, reported for reference). Cartographic roles are driven by the WEIGHTED WMDz (WM_strength_z) instead. NA for modules < 3 regions.",
               "Participation coefficient (binary, canonical): 0 = all edges within own module, ->1 = edges spread across modules. Module-size bias (Pedersen et al. 2020).",
               "Strength-based weighted variants of WMDz and PC.",
               "Guimera-Amaral role (R1-R7) on the WEIGHTED cartography (WM_strength_z, Participation_coef_wt). Hub/non-hub boundary = wmdz_hub_cut = 1.5 (lower bound of the canonical 1.5-2.5 high-WMDz range; Cohen & D'Esposito, cited by Kimbrough 2020). WM_strength_z is a unit-SD within-module z-score here, so the canonical range applies directly. Multi-centrality hubs below are reported as convergent evidence.",
               "Hub identification by convergence across centrality measures (Sporns et al. 2007; van den Heuvel & Sporns 2013): strength, weighted betweenness, weighted closeness, eigenvector. A region is 'high' on a measure if > mean + 1 SD within its network.",
               "Hub_score = number of the four centralities a region is high on (0-4); Hub = TRUE when high on >= hub_min_metrics (default 2). Computed at matched density so hub status is comparable across conditions.",
               "Mean within-network z across the four centralities (continuous hub-ness, for ranking and cross-condition comparison).",
               "All region x condition x sex rows flagged Hub = TRUE, with their centralities and roles.",
               "PRIMARY hub-change analysis (Combined networks): Hub_score per condition, number of conditions each region is a hub, and gained/lost hub status vs Vehicle.",
               "TRUE if home module >= min_module_size; FALSE flags roles to read with caution.",
               "Anatomical (CCF major-structure) composition of each network module per condition x sex.",
               "Per region, proportion of conditions it stays in its modal network module (consensus modules). Louvain module numbers are arbitrary within each network, so each condition's partition is first relabelled onto the reference condition (Matched_to) by greedy maximum overlap; without that step the column would be comparing meaningless labels.",
               "Module count from the standard cFos hierarchical-clustering method (Euclidean distance on correlation profiles, complete linkage; Kimbrough 2020, Ardinger 2024) at dendrogram cut heights 0.3-0.7 of max tree height, shown next to the consensus-Louvain count. Documents that the lineage method collapses to very few modules on these strongly-coactivated networks, justifying consensus Louvain.",
               "Single-sex networks (Males / Females) are the primary, sex-stratified analysis (matching the sex-difference focus of Ardinger 2024 / Bloch 2024). The pooled Combined network (n = 10-12, sex residualized) is retained for resampling-based inference, where per-sex n (5-6) is too low for stable bootstrap/permutation."))
if (length(nodes_list) > 0) {
  node_all <- bind_rows(nodes_list) %>%
    dplyr::select(Treatment, Sex, Region, Structure, Network_module, Home_module_size,
                  WM_degree_z, Participation_coef, WM_strength_z, Participation_coef_wt,
                  Node_role, Role_reliable, Degree, Strength, Weighted_betweenness,
                  Weighted_closeness, Eigenvector, Local_weighted_transitivity,
                  Hub_score, Centrality_composite, Hub)
  add_sheet(wb3, "Node_Roles", node_all)
  hub_regions <- node_all %>% filter(Hub) %>%
    arrange(Sex != "Combined", Treatment, desc(Centrality_composite)) %>%
    dplyr::select(Treatment, Sex, Region, Structure, Hub_score, Centrality_composite,
                  Strength, Weighted_betweenness, Weighted_closeness, Eigenvector, Node_role)
  if (nrow(hub_regions) > 0) add_sheet(wb3, "Hub_Regions", hub_regions)
}
if (!is.null(hub_changes) && nrow(hub_changes) > 0)
  add_sheet(wb3, "Hub_Changes_Combined", hub_changes)
if (length(modules_list) > 0) add_sheet(wb3, "Module_Assignments",
                                        bind_rows(modules_list) %>% dplyr::select(Treatment, Sex, Region, Structure, Network_module))
# Module anatomical composition
comp_list <- list()
for (lab in names(modules_list)) {
  md <- modules_list[[lab]]
  comp_list[[lab]] <- md %>% group_by(Treatment, Sex, Network_module, Structure) %>%
    summarise(N = n(), .groups = "drop") %>% group_by(Treatment, Sex, Network_module) %>%
    mutate(Pct = round(100 * N / sum(N), 1)) %>% arrange(Treatment, Sex, Network_module, desc(N))
}
if (length(comp_list) > 0) add_sheet(wb3, "Module_Composition", bind_rows(comp_list))
if (length(stability_list) > 0) add_sheet(wb3, "Module_Stability", bind_rows(stability_list))
if (length(modmethod_list) > 0) add_sheet(wb3, "Module_Method_Comparison", bind_rows(modmethod_list))
saveWorkbook(wb3, file.path(output_root, "03_Modules_and_Roles.xlsx"), overwrite = TRUE)

writeLines(capture.output(sessionInfo()),
           file.path(output_root, "sessionInfo_network.txt"))

cat("\nAll workbooks written to:", output_root, "\n")
cat("  00_Processed_Data.xlsx | 01_Connectivity.xlsx | 02_Network_Metrics.xlsx | 03_Modules_and_Roles.xlsx\n")
cat("Done.\n")
