# regional_analysis.R
# =============================================================================
# Whole-brain cFos regional activation analysis.
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
#   Rscript scripts/regional_analysis.R
# Input and output locations are set in the Config section below; see README.
# =============================================================================
#
# ANALYSIS PRIORITY:
#   PRIMARY   - One-way ANOVA on pooled males + females (condition only, no sex term)
#               Output: Drug_Statistical_Results_Primary.xlsx
#
#   SECONDARY - Sex-stratified one-way ANOVA (males and females separately)
#               Output: Drug_Statistical_Results_BySex.xlsx
#
#   SUPPLEMENTARY - Two-way ANOVA (condition x sex), Negative Binomial, and a
#               global-activation sensitivity analysis
#               Output: Drug_Statistical_Results_TwoWay.xlsx
#                       Drug_Statistical_Results_NegBin.xlsx
#                       Drug_Statistical_Results_GlobalNorm.xlsx
#                       Model_Assumption_Checks.xlsx
#
# DATA STRUCTURE:
#   Drug_Composite.xlsx: "Clustered Males" and "Clustered Females" sheets
#   198 CCFv3 regions, 4 conditions, n=5-6 per group per sex
#
# NORMALIZATION:
#   ANOVA methods: log10(cells/mm3) using CCFv3 atlas volumes (Wang et al. 2020).
#     Density is the field-standard outcome for whole-brain cFos mapping (Kimbrough
#     et al. 2020; Hu et al.; Ishii et al.) and is the PRIMARY analysis here.
#   NegBin method: bilateral INTEGER counts (hemisphere average x 2) with
#     log(bilateral volume) as offset. See bilateral_count() for why the x2 matters.
#   GLOBAL-NORM sensitivity: each region's count is divided by that ANIMAL'S
#     whole-brain total before log10. This is NOT the primary analysis and is not
#     proposed as a replacement for density -- it exists to answer the reviewer
#     question "is the regional pattern selective, or a rescaling of a global
#     difference in total activation?", which matters here because whole-brain
#     totals differ by condition in males (see Figure 2B).
#
# ZERO COUNTS:
#   log10(0) is undefined, so an animal with 0 cells in a region is EXCLUDED from
#   that region's ANOVA statistics (omnibus_F/p, estimate, std_error, t_statistic,
#   contrast_F, cohens_d, p_tukey). The raw-count columns (group_a_mean,
#   group_b_mean, fold_change, direction) are means over ALL animals, zeros
#   included, so they describe the complete group. The two therefore rest on
#   different animals whenever a group contains a zero: n_a / n_b report how many
#   animals entered the test and n_zero_a / n_zero_b how many were dropped, so
#   every such row is identifiable. The negative-binomial supplement models
#   bilateral counts and keeps the zeros, so its n_a / n_b are the full group sizes.
#   In this dataset zeros are rare but sex-asymmetric: 1 value in males, 17 in
#   females across 11 regions (incl. ASO and AVPV, both plotted). Report that.
#   (network_analysis.R takes the other route for the network inputs and imputes
#   zeros at half the region's minimum non-zero value; see its header.)
#
# EFFECT DIRECTION CONVENTION:
#   For every comparison (A vs B), all signed quantities are reported as
#   (B - A) = treatment - reference, so a treatment-induced INCREASE is POSITIVE.
#   This applies to estimate, t_statistic / z_statistic, cohens_d, and (via
#   exp) IRR, and is uniform across the primary, by-sex, two-way, and NegBin
#   tables. It matches the signed Cohen's d convention used in Figure 3C.
#   fold_change is likewise treatment/reference (mean_B / mean_A). Tukey p-values
#   are two-sided and unaffected by direction. The `direction` column reports the
#   raw-count ordering (A > B or A < B).
#
# PAIRWISE COMPARISONS (5 of 6 - Morphine-Dependent vs Ro excluded):
#   Vehicle vs Morphine, Vehicle vs Morphine-Dependent, Vehicle vs Ro,
#   Morphine vs Morphine-Dependent, Morphine vs Ro
#
# Required packages:
#   install.packages(c("readxl", "openxlsx", "dplyr", "MASS", "emmeans"))

# ---- Packages ----------------------------------------------------------------
required_pkgs <- c("readxl", "openxlsx", "dplyr", "MASS", "emmeans")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages. Please run:\n  install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}
library(readxl); library(openxlsx); library(dplyr); library(MASS); library(emmeans)

# ---- Config ------------------------------------------------------------------
# Paths are relative to the repository root. Run from the repo root, e.g.:
#   Rscript scripts/regional_analysis.R
# To use a different layout, override the two directories below, either by
# editing them here or by setting the environment variables before running.
DATA_DIR    <- Sys.getenv("CFOS_DATA_DIR",   "data")
OUTPUT_DIR  <- Sys.getenv("CFOS_OUTPUT_DIR", "results")
RAW_FILE    <- file.path(DATA_DIR, "Drug_Composite.xlsx")
VOLUME_FILE <- file.path(DATA_DIR, "ccfv3_volumes.xlsx")
FDR_ALPHA   <- 0.05

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Load data ---------------------------------------------------------------
cat("Loading data...\n")
df_m <- as.data.frame(read_excel(RAW_FILE, sheet = "Clustered Males"))
df_f <- as.data.frame(read_excel(RAW_FILE, sheet = "Clustered Females"))

get_data_cols <- function(df) {
  cols <- grep("Vehicle|Morphine|Ro", names(df), value = TRUE)
  cols[!grepl("^Unnamed", cols)]
}
data_cols_m <- get_data_cols(df_m)
data_cols_f <- get_data_cols(df_f)

make_groups <- function(cols) {
  list(
    Vehicle            = grep("^Vehicle",  cols, value = TRUE),
    Morphine           = grep("^Morphine", cols, value = TRUE),
    `Morphine-Dependent` = grep("^Chronic",  cols, value = TRUE),
    Ro                 = grep("^Ro",       cols, value = TRUE)
  )
}
groups_m   <- make_groups(data_cols_m)
groups_f   <- make_groups(data_cols_f)
cond_names <- names(groups_m)

cat(sprintf("  Males:   %d regions, %d animals\n", nrow(df_m), length(data_cols_m)))
cat(sprintf("  Females: %d regions, %d animals\n", nrow(df_f), length(data_cols_f)))

# ---- Column-integrity guard --------------------------------------------------
# A single header typo (e.g. "Morphne" instead of "Morphine") makes a real
# animal invisible to get_data_cols()/make_groups() and it is dropped silently,
# corrupting every statistic for that condition. Any header that looks like an
# animal column (contains "Male"/"Female" followed by a number) MUST be both
# captured and assigned to exactly one condition group, otherwise the run aborts with an error.
check_columns <- function(df, data_cols, groups, sex_label) {
  looks_data <- grep("(Male|Female)[[:space:]]*[0-9]", names(df), value = TRUE)
  grouped    <- unlist(groups, use.names = FALSE)
  dropped    <- setdiff(looks_data, data_cols)   # look like data but not captured
  ungrouped  <- setdiff(data_cols,  grouped)     # captured but no condition group
  if (length(dropped) || length(ungrouped)) {
    stop(sprintf(paste0(
      "[%s] column-integrity check FAILED.\n",
      "  Look like data but not captured: %s\n",
      "  Captured but not grouped:        %s\n",
      "  Fix the header typo (e.g. 'Morphne' -> 'Morphine') and re-run."),
      sex_label,
      if (length(dropped))   paste(dropped,   collapse = ", ") else "none",
      if (length(ungrouped)) paste(ungrouped, collapse = ", ") else "none"))
  }
  cat(sprintf("  [%s] integrity OK: %d columns grouped (%s)\n",
              sex_label, length(data_cols),
              paste(sprintf("%s=%d", names(groups), lengths(groups)), collapse = ", ")))
}
check_columns(df_m, data_cols_m, groups_m, "Males")
check_columns(df_f, data_cols_f, groups_f, "Females")

# ---- Row-alignment guard -----------------------------------------------------
# The pooled, two-way and negative-binomial sections read the female sheet by the
# MALE row number (df_f[i, ] alongside df_m[i, ]), so both sheets must list the
# same regions in the same order. A row inserted, deleted or re-sorted in one
# sheet would silently pair one region's male counts with another's female
# counts, so the run aborts rather than producing quietly wrong statistics.
if (!identical(df_m$acronym, df_f$acronym)) {
  only_m <- setdiff(df_m$acronym, df_f$acronym)
  only_f <- setdiff(df_f$acronym, df_m$acronym)
  mism   <- if (length(df_m$acronym) == length(df_f$acronym))
              which(df_m$acronym != df_f$acronym) else integer(0)
  stop(sprintf(paste0(
    "row-alignment check FAILED: 'Clustered Males' and 'Clustered Females' must\n",
    "  list the same regions in the same order.\n",
    "  Males only:           %s\n",
    "  Females only:         %s\n",
    "  First order mismatch: %s\n",
    "  Fix the sheets so both acronym columns are identical and re-run."),
    if (length(only_m)) paste(only_m, collapse = ", ") else "none",
    if (length(only_f)) paste(only_f, collapse = ", ") else "none",
    if (length(mism)) sprintf("row %d (%s vs %s)", mism[1],
                              df_m$acronym[mism[1]], df_f$acronym[mism[1]])
    else "n/a (row counts differ)"))
}
cat(sprintf("  Row alignment OK: %d regions, identical order in both sheets\n",
            nrow(df_m)))

# ---- Volumes (Wang et al. 2020 CCFv3 Table S4) --------------------------------
cat("Loading volumes...\n")
vol_df          <- as.data.frame(read_excel(VOLUME_FILE, skip = 1))
vol_map         <- setNames(as.numeric(vol_df[["Mean Volume (m)"]]),
                            trimws(as.character(vol_df[["abbreviation"]])))
df_m$volume_mm3 <- vol_map[trimws(df_m$acronym)]
df_f$volume_mm3 <- vol_map[trimws(df_f$acronym)]
cat(sprintf("  Matched %d/%d regions\n", sum(!is.na(df_m$volume_mm3)), nrow(df_m)))

# Volume-match guard. An unmatched acronym gives volume = NA, hence all-NA
# densities and a silently empty result row rather than an error. Fail loudly.
.unmatched <- df_m$acronym[is.na(df_m$volume_mm3) | df_m$volume_mm3 <= 0]
if (length(.unmatched)) {
  stop(sprintf(paste0(
    "volume lookup FAILED for %d region(s): %s\n",
    "  Check spelling/whitespace in the 'abbreviation' column of %s,\n",
    "  and confirm 'Mean Volume (m)' is in mm^3."),
    length(.unmatched), paste(.unmatched, collapse = ", "), VOLUME_FILE))
}

# ---- Log10 density normalization ---------------------------------------------
log10_density <- function(df, data_cols) {
  out <- df[, c("name", "acronym")]
  for (col in data_cols) {
    vals <- df[[col]] / df$volume_mm3
    vals[vals == 0] <- NA
    out[[col]] <- log10(vals)
  }
  out
}
norm_m <- log10_density(df_m, data_cols_m)
norm_f <- log10_density(df_f, data_cols_f)

# Build pooled normalized data (males + females combined, same column names preserved)
norm_comb   <- df_m[, c("name", "acronym")]
groups_comb <- list()
df_comb     <- df_m[, c("name", "acronym", "volume_mm3")]
for (g in cond_names) {
  for (col in groups_m[[g]]) { norm_comb[[col]] <- norm_m[[col]]; df_comb[[col]] <- df_m[[col]] }
  for (col in groups_f[[g]]) { norm_comb[[col]] <- norm_f[[col]]; df_comb[[col]] <- df_f[[col]] }
  groups_comb[[g]] <- c(groups_m[[g]], groups_f[[g]])
}
cat("Normalization: log10(cells/mm3)\n")

# ---- Helpers -----------------------------------------------------------------
# 5 biologically relevant pairs - Morphine-Dependent vs Ro excluded
pairs <- list(
  c("Vehicle",  "Morphine"),
  c("Vehicle",  "Morphine-Dependent"),
  c("Vehicle",  "Ro"),
  c("Morphine", "Morphine-Dependent"),
  c("Morphine", "Ro")
)

# emmeans::contrast(method="pairwise") returns rows in combn(levels) order:
# (1,2),(1,3),(1,4),(2,3),(2,4),(3,4) for levels in the given order. We select
# the row for a pair by POSITION, not by matching the printed contrast label.
# Label matching with grepl is unsafe here because "Morphine" is a substring of
# "Morphine-Dependent", so a fallback search silently returns the wrong contrast.
pairwise_row <- function(ga, gb, levs) {
  ia <- match(ga, levs); ib <- match(gb, levs)
  if (is.na(ia) || is.na(ib)) return(NA_integer_)
  lo <- min(ia, ib); hi <- max(ia, ib)
  combos <- utils::combn(length(levs), 2)      # columns are (i, j), i < j
  which(combos[1, ] == lo & combos[2, ] == hi)
}

# Recover the bilateral integer count from a hemisphere-averaged value.
# Values in Drug_Composite.xlsx are (left + right) / 2, so they are integers or
# half-integers; x2 is exactly the bilateral count. Anything else means the input
# is not a hemisphere-averaged count and the negative-binomial model would be
# fitting rounded, non-count data -- so stop rather than round silently.
bilateral_count <- function(x) {
  y <- x * 2
  if (abs(y - round(y)) > 1e-8)
    stop(sprintf(paste0(
      "value %.4f is not a hemisphere-averaged count (x2 = %.4f is not an integer).\n",
      "  The negative-binomial supplement models bilateral integer counts; check\n",
      "  how Drug_Composite.xlsx was generated before proceeding."), x, y))
  as.integer(round(y))
}

sig_stars <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "ns",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.1   ~ "\u2020",
    TRUE      ~ "ns"
  )
}

# Count regions below a significance threshold, read from the p-value column
# directly so the count does not depend on the significance-marker encoding.
n_below <- function(p, alpha) sum(p < alpha, na.rm = TRUE)

add_fdr <- function(res, p_col, fdr_col, star_col) {
  valid <- !is.na(res[[p_col]])
  res[[fdr_col]] <- NA
  res[[fdr_col]][valid] <- p.adjust(res[[p_col]][valid], method = "BH")
  res[[star_col]] <- sig_stars(res[[fdr_col]])
  res
}

write_excel <- function(results_list, filename, sig_col = "sig_tukey") {
  wb  <- createWorkbook()
  hdr <- createStyle(fontName="Arial", fontSize=10, textDecoration="bold",
                     border="Bottom", borderColour="#000000")
  dat <- createStyle(fontName="Arial", fontSize=10)
  sig <- createStyle(fontName="Arial", fontSize=10, fgFill="#FFF2CC")
  for (label in names(results_list)) {
    res   <- results_list[[label]]
    sname <- gsub("Morphine-Dependent", "MorDep", label)
    sname <- gsub("Morphine", "Mor", sname)
    sname <- gsub("Vehicle",  "Veh", sname)
    sname <- gsub("Females",  "Fem", sname)
    sname <- substr(sname, 1, 31)
    addWorksheet(wb, sname)
    writeData(wb, sname, res, startRow=1, startCol=1, rowNames=FALSE)
    addStyle(wb, sname, hdr, rows=1, cols=1:ncol(res), gridExpand=TRUE)
    setRowHeights(wb, sname, rows=1, heights=18)
    for (i in seq_len(nrow(res))) {
      sig_val <- res[[sig_col]][i]
      is_sig  <- isTRUE(!is.na(sig_val) && sig_val != "ns")
      addStyle(wb, sname, if (is_sig) sig else dat,
               rows=i+1, cols=1:ncol(res), gridExpand=TRUE)
    }
    setColWidths(wb, sname, cols=1:ncol(res), widths="auto")
    freezePane(wb, sname, firstActiveRow=2)
  }
  saveWorkbook(wb, filename, overwrite=TRUE)
  cat(sprintf("  Saved: %s\n", basename(filename)))
}

# Shared one-way ANOVA function - used by both primary and by-sex analyses
run_oneway <- function(norm_df, df_raw_sex, groups_sex, label_suffix = "") {
  out <- list()
  for (pair in pairs) {
    ga <- pair[1]; gb <- pair[2]
    label   <- if (nchar(label_suffix)) paste(ga, "vs", gb, "-", label_suffix)
    else paste(ga, "vs", gb)
    cols_a  <- groups_sex[[ga]]; cols_b <- groups_sex[[gb]]
    
    rows_list <- vector("list", nrow(norm_df))
    n_err     <- 0
    for (i in seq_len(nrow(norm_df))) {
      acr   <- norm_df$acronym[i]; rgn <- norm_df$name[i]
      va    <- as.numeric(norm_df[i, cols_a, drop=TRUE]);    va    <- va[!is.na(va)]
      vb    <- as.numeric(norm_df[i, cols_b, drop=TRUE]);    vb    <- vb[!is.na(vb)]
      raw_a <- as.numeric(df_raw_sex[i, cols_a, drop=TRUE]); raw_a <- raw_a[!is.na(raw_a)]
      raw_b <- as.numeric(df_raw_sex[i, cols_b, drop=TRUE]); raw_b <- raw_b[!is.na(raw_b)]
      mean_a_raw <- if (length(raw_a)) mean(raw_a) else NA
      mean_b_raw <- if (length(raw_b)) mean(raw_b) else NA
      
      all_v <- c(); all_l <- c()
      for (g in cond_names) {
        v <- as.numeric(norm_df[i, groups_sex[[g]], drop=TRUE]); v <- v[!is.na(v)]
        all_v <- c(all_v, v); all_l <- c(all_l, rep(g, length(v)))
      }
      
      empty <- data.frame(region_name=rgn, acronym=acr, comparison=paste(ga,"-",gb),
                          omnibus_F=NA, omnibus_p=NA, df_num=NA, df_den=NA,
                          contrast_F=NA,
                          estimate=NA, std_error=NA, ci_low=NA, ci_high=NA,
                          t_statistic=NA, p_unadj=NA, p_tukey=NA,
                          fold_change=NA, cohens_d=NA, hedges_g=NA,
                          group_a_mean=mean_a_raw, group_b_mean=mean_b_raw,
                          n_a=length(va), n_b=length(vb),
                          n_zero_a=sum(raw_a == 0), n_zero_b=sum(raw_b == 0),
                          direction=NA, stringsAsFactors=FALSE)

      if (length(va) < 2 || length(vb) < 2) { rows_list[[i]] <- empty; next }

      tryCatch({
        # Explicit factor levels so Tukey row lookup is deterministic (never
        # dependent on R's default alphabetical ordering).
        fit   <- aov(all_v ~ factor(all_l, levels = cond_names))
        aov_s <- summary(fit)[[1]]
        ms_e  <- aov_s[["Mean Sq"]][2]
        df_n  <- aov_s[["Df"]][1]     # condition df  -> report as F(df_num, df_den)
        df_d  <- aov_s[["Df"]][2]     # residual df; varies with zero exclusion
        tukey <- TukeyHSD(fit)[[1]]
        # p adj is two-sided (direction-independent); take it from whichever
        # row name matches this pair.
        key1  <- paste0(gb,"-",ga); key2 <- paste0(ga,"-",gb)
        p_tukey <- if (key1 %in% rownames(tukey)) tukey[key1, "p adj"] else
          if (key2 %in% rownames(tukey)) tukey[key2, "p adj"] else NA
        
        # All signed quantities: treatment - reference = (B - A), increase-positive.
        n_a <- length(va); n_b <- length(vb)
        estimate <- mean(vb) - mean(va)            # log10 mean difference (B - A)
        se   <- sqrt(ms_e * (1/n_a + 1/n_b))
        t_s  <- if (!is.na(estimate) && se > 0) estimate / se else NA
        fc   <- if (!is.na(mean_a_raw) && mean_a_raw > 0) mean_b_raw / mean_a_raw else NA
        sd_p <- sqrt(((n_a-1)*var(va) + (n_b-1)*var(vb)) / (n_a+n_b-2))
        cd   <- if (sd_p > 0) (mean(vb) - mean(va)) / sd_p else NA   # Cohen's d, B - A
        # Hedges' g: Cohen's d with the small-sample bias correction. At n=5-6 per
        # group (the sex-stratified analysis) d is inflated by ~8%, so report g
        # alongside it rather than leaving the correction for a reviewer to raise.
        hg   <- if (!is.na(cd)) cd * (1 - 3 / (4*(n_a + n_b) - 9)) else NA
        # Uncorrected two-sided p for the contrast, on the ANOVA's residual df.
        # This is the genuinely UNCORRECTED pairwise p. p_tukey below is already
        # adjusted across the 6 within-region pairwise comparisons, so the two must
        # never be described interchangeably in the manuscript.
        p_un <- if (!is.na(t_s) && !is.na(df_d) && df_d > 0)
                  2 * stats::pt(-abs(t_s), df = df_d) else NA
        tcrit <- if (!is.na(df_d) && df_d > 0) stats::qt(0.975, df = df_d) else NA
        ci_l  <- if (!is.na(tcrit) && !is.na(se)) estimate - tcrit * se else NA
        ci_h  <- if (!is.na(tcrit) && !is.na(se)) estimate + tcrit * se else NA
        dir  <- if (!is.na(mean_a_raw) && mean_a_raw > mean_b_raw)
          paste(ga,">",gb) else paste(ga,"<",gb)

        rows_list[[i]] <- data.frame(region_name=rgn, acronym=acr,
                                     comparison=paste(ga,"-",gb),
                                     omnibus_F=aov_s[["F value"]][1], omnibus_p=aov_s[["Pr(>F)"]][1],
                                     df_num=df_n, df_den=df_d,
                                     contrast_F=if(!is.na(t_s)) t_s^2 else NA,
                                     estimate=estimate, std_error=se,
                                     ci_low=ci_l, ci_high=ci_h,
                                     t_statistic=t_s, p_unadj=p_un, p_tukey=p_tukey,
                                     fold_change=fc, cohens_d=cd, hedges_g=hg,
                                     group_a_mean=mean_a_raw, group_b_mean=mean_b_raw,
                                     n_a=n_a, n_b=n_b,
                                     n_zero_a=sum(raw_a == 0), n_zero_b=sum(raw_b == 0),
                                     direction=dir, stringsAsFactors=FALSE)
      }, error=function(e) {
        rows_list[[i]] <<- empty; n_err <<- n_err + 1
        if (n_err <= 5) message(sprintf("    [%s] %s: %s", label, acr, conditionMessage(e)))
      })
    }
    
    res <- do.call(rbind, rows_list)
    valid_t <- !is.na(res$p_tukey);  res$p_fdr <- NA
    res$p_fdr[valid_t]    <- p.adjust(res$p_tukey[valid_t], method="BH")
    valid_a <- !is.na(res$omnibus_p); res$omnibus_fdr <- NA
    res$omnibus_fdr[valid_a] <- p.adjust(res$omnibus_p[valid_a], method="BH")
    res$sig_tukey       <- sig_stars(res$p_tukey)
    res$sig_fdr         <- sig_stars(res$p_fdr)
    res$sig_omnibus_fdr <- sig_stars(res$omnibus_fdr)
    
    out[[label]] <- res
    cat(sprintf("    %-40s  Tukey: %3d (%3d)   FDR: %3d (%3d)%s\n", label,
                n_below(res$p_tukey, 0.05), n_below(res$p_tukey, 0.10),
                n_below(res$p_fdr,   0.05), n_below(res$p_fdr,   0.10),
                if (n_err) sprintf("   [%d region(s) errored]", n_err) else ""))
  }
  out
}

# ==============================================================================
# PRIMARY: One-way ANOVA - males + females pooled, condition only
# ==============================================================================
cat("\n--- PRIMARY: One-way ANOVA (pooled, no sex term) ---\n")
primary_results <- run_oneway(norm_comb, df_comb, groups_comb, label_suffix = "")
write_excel(primary_results,
            file.path(OUTPUT_DIR, "Drug_Statistical_Results_Primary.xlsx"))

# ==============================================================================
# SECONDARY: Sex-stratified one-way ANOVA
# ==============================================================================
cat("\n--- SECONDARY: Sex-stratified one-way ANOVA ---\n")
bysex_results <- c(
  run_oneway(norm_m, df_m, groups_m, label_suffix = "Males"),
  run_oneway(norm_f, df_f, groups_f, label_suffix = "Females")
)
write_excel(bysex_results,
            file.path(OUTPUT_DIR, "Drug_Statistical_Results_BySex.xlsx"))

# ==============================================================================
# SUPPLEMENTARY 1: Two-way ANOVA - condition x sex
# ==============================================================================
cat("\n--- SUPPLEMENTARY: Two-way ANOVA (condition x sex) ---\n")

make_long_twoway <- function(i) {
  rows <- list()
  for (g in cond_names) {
    for (col in groups_m[[g]]) {
      v <- norm_m[i, col]
      if (!is.na(v)) rows[[length(rows)+1]] <-
          data.frame(value=v, condition=g, sex="Male",   stringsAsFactors=FALSE)
    }
    for (col in groups_f[[g]]) {
      v <- norm_f[i, col]
      if (!is.na(v)) rows[[length(rows)+1]] <-
          data.frame(value=v, condition=g, sex="Female", stringsAsFactors=FALSE)
    }
  }
  long <- do.call(rbind, rows)
  long$condition <- factor(long$condition, levels=cond_names)
  long$sex       <- factor(long$sex, levels=c("Male","Female"))
  long
}

twoway_results <- list()
sex_rows       <- vector("list", nrow(norm_m))
int_rows       <- vector("list", nrow(norm_m))

for (i in seq_len(nrow(norm_m))) {
  acr  <- norm_m$acronym[i]; rgn <- norm_m$name[i]
  long <- make_long_twoway(i)
  sex_rows[[i]] <- data.frame(region_name=rgn, acronym=acr, F_sex=NA, p_sex=NA, stringsAsFactors=FALSE)
  int_rows[[i]] <- data.frame(region_name=rgn, acronym=acr, F_interaction=NA, p_interaction=NA, stringsAsFactors=FALSE)
  if (nrow(long) < 6) next
  tryCatch({
    fit    <- aov(value ~ condition * sex, data=long)
    summ   <- summary(fit)[[1]]
    rnames <- trimws(rownames(summ))
    sex_idx <- which(rnames == "sex")
    int_idx <- which(rnames == "condition:sex")
    sex_rows[[i]] <- data.frame(region_name=rgn, acronym=acr,
                                F_sex = if (length(sex_idx)) summ[["F value"]][sex_idx] else NA,
                                p_sex = if (length(sex_idx)) summ[["Pr(>F)"]][sex_idx]  else NA,
                                stringsAsFactors=FALSE)
    int_rows[[i]] <- data.frame(region_name=rgn, acronym=acr,
                                F_interaction = if (length(int_idx)) summ[["F value"]][int_idx] else NA,
                                p_interaction = if (length(int_idx)) summ[["Pr(>F)"]][int_idx]  else NA,
                                stringsAsFactors=FALSE)
  }, error=function(e) {
    sex_rows[[i]] <<- data.frame(region_name=rgn, acronym=acr, F_sex=NA, p_sex=NA, stringsAsFactors=FALSE)
    int_rows[[i]] <<- data.frame(region_name=rgn, acronym=acr, F_interaction=NA, p_interaction=NA, stringsAsFactors=FALSE)
  })
}

sex_effect_df <- add_fdr(do.call(rbind, sex_rows), "p_sex",         "q_sex",         "sig_sex")
int_effect_df <- add_fdr(do.call(rbind, int_rows), "p_interaction", "q_interaction", "sig_interaction")
cat(sprintf("  Sex main effect:               %3d (%3d) regions\n",
            n_below(sex_effect_df$q_sex, 0.05), n_below(sex_effect_df$q_sex, 0.10)))
cat(sprintf("  Condition x Sex interaction:   %3d (%3d) regions\n",
            n_below(int_effect_df$q_interaction, 0.05), n_below(int_effect_df$q_interaction, 0.10)))

for (pair in pairs) {
  ga <- pair[1]; gb <- pair[2]; label <- paste(ga, "vs", gb)
  rows_list <- vector("list", nrow(norm_m))
  n_err     <- 0
  
  for (i in seq_len(nrow(norm_m))) {
    acr  <- norm_m$acronym[i]; rgn <- norm_m$name[i]
    long <- make_long_twoway(i)
    # log10 densities entering the model (zeros became NA upstream)
    va <- c(as.numeric(norm_m[i, groups_m[[ga]], drop=TRUE]),
            as.numeric(norm_f[i, groups_f[[ga]], drop=TRUE]))
    vb <- c(as.numeric(norm_m[i, groups_m[[gb]], drop=TRUE]),
            as.numeric(norm_f[i, groups_f[[gb]], drop=TRUE]))
    va <- va[!is.na(va)]; vb <- vb[!is.na(vb)]
    raw_a <- c(as.numeric(df_m[i, groups_m[[ga]], drop=TRUE]),
               as.numeric(df_f[i, groups_f[[ga]], drop=TRUE]))
    raw_b <- c(as.numeric(df_m[i, groups_m[[gb]], drop=TRUE]),
               as.numeric(df_f[i, groups_f[[gb]], drop=TRUE]))
    raw_a <- raw_a[!is.na(raw_a)]; raw_b <- raw_b[!is.na(raw_b)]
    mean_a <- if (length(raw_a)) mean(raw_a) else NA
    mean_b <- if (length(raw_b)) mean(raw_b) else NA
    
    empty <- data.frame(region_name=rgn, acronym=acr, comparison=paste(ga,"-",gb),
                        omnibus_F=NA, omnibus_p=NA, contrast_F=NA,
                        estimate=NA, std_error=NA, t_statistic=NA, p_tukey=NA,
                        fold_change=NA, cohens_d=NA,
                        group_a_mean=mean_a, group_b_mean=mean_b,
                        n_a=length(va), n_b=length(vb),
                        n_zero_a=sum(raw_a == 0), n_zero_b=sum(raw_b == 0),
                        direction=NA, stringsAsFactors=FALSE)
    
    if (nrow(long) < 6) { rows_list[[i]] <- empty; next }
    
    tryCatch({
      fit    <- aov(value ~ condition * sex, data=long)
      aov_s  <- summary(fit)[[1]]
      rnames <- trimws(rownames(aov_s))
      cidx   <- which(rnames == "condition")
      omnibus_F <- if (length(cidx)) aov_s[["F value"]][cidx] else NA
      omnibus_p <- if (length(cidx)) aov_s[["Pr(>F)"]][cidx]  else NA
      
      em    <- emmeans::emmeans(fit, ~ condition)
      con   <- emmeans::contrast(em, method="pairwise", adjust="tukey")
      cdf   <- as.data.frame(con)
      # Select by position among the levels the MODEL actually fitted. A condition
      # with no observations in this region is dropped from the fit, so emmeans
      # returns fewer than 6 contrast rows; indexing on cond_names would then
      # silently return a different comparison.
      emm_lv <- as.character(as.data.frame(em)$condition)
      ridx   <- pairwise_row(ga, gb, emm_lv)
      crow  <- if (!is.na(ridx) && ridx <= nrow(cdf)) cdf[ridx, , drop=FALSE] else cdf[0, ]
      
      # Direction-safe estimate straight from the marginal means: B - A
      # (treatment - reference), independent of the contrast's internal ordering.
      emm_df <- as.data.frame(em)
      m_a    <- emm_df$emmean[emm_df$condition == ga]
      m_b    <- emm_df$emmean[emm_df$condition == gb]
      estimate <- if (length(m_a) && length(m_b)) m_b - m_a else NA
      se       <- if (nrow(crow)) crow$SE[1]      else NA
      p_tukey  <- if (nrow(crow)) crow$p.value[1] else NA
      t_s      <- if (!is.na(estimate) && !is.na(se) && se > 0) estimate / se else NA
      
      n_a <- length(va); n_b <- length(vb)
      sd_p <- sqrt(((n_a-1)*var(va)+(n_b-1)*var(vb))/(n_a+n_b-2))
      cd   <- if (!is.na(sd_p) && sd_p > 0) (mean(vb)-mean(va))/sd_p else NA  # B - A
      fc   <- if (!is.na(mean_a) && mean_a > 0) mean_b/mean_a else NA
      dir  <- if (!is.na(mean_a) && !is.na(mean_b)) {
        if (mean_a > mean_b) paste(ga,">",gb) else paste(ga,"<",gb)
      } else NA
      
      rows_list[[i]] <- data.frame(region_name=rgn, acronym=acr,
                                   comparison=paste(ga,"-",gb),
                                   omnibus_F=omnibus_F, omnibus_p=omnibus_p,
                                   contrast_F=if(!is.na(t_s)) t_s^2 else NA,
                                   estimate=estimate, std_error=se, t_statistic=t_s, p_tukey=p_tukey,
                                   fold_change=fc, cohens_d=cd,
                                   group_a_mean=mean_a, group_b_mean=mean_b,
                                   n_a=n_a, n_b=n_b,
                                   n_zero_a=sum(raw_a == 0), n_zero_b=sum(raw_b == 0),
                                   direction=dir, stringsAsFactors=FALSE)
    }, error=function(e) { rows_list[[i]] <<- empty; n_err <<- n_err + 1 })
  }
  
  res <- do.call(rbind, rows_list)
  valid_t <- !is.na(res$p_tukey);   res$p_fdr <- NA
  res$p_fdr[valid_t] <- p.adjust(res$p_tukey[valid_t], method="BH")
  valid_a <- !is.na(res$omnibus_p); res$omnibus_fdr <- NA
  res$omnibus_fdr[valid_a] <- p.adjust(res$omnibus_p[valid_a], method="BH")
  res$sig_tukey       <- sig_stars(res$p_tukey)
  res$sig_fdr         <- sig_stars(res$p_fdr)
  res$sig_omnibus_fdr <- sig_stars(res$omnibus_fdr)
  
  twoway_results[[label]] <- res
  cat(sprintf("    %-35s  Tukey: %3d (%3d)   FDR: %3d (%3d)%s\n", label,
              n_below(res$p_tukey, 0.05), n_below(res$p_tukey, 0.10),
              n_below(res$p_fdr,   0.05), n_below(res$p_fdr,   0.10),
              if (n_err) sprintf("   [%d region(s) errored]", n_err) else ""))
}

twoway_results[["Sex Main Effect"]]          <- sex_effect_df
twoway_results[["Condition x Sex Interact"]] <- int_effect_df
write_excel(twoway_results,
            file.path(OUTPUT_DIR, "Drug_Statistical_Results_TwoWay.xlsx"))

# ==============================================================================
# SUPPLEMENTARY 2: Negative binomial - condition x sex, raw counts
# ==============================================================================
cat("\n--- SUPPLEMENTARY: Negative binomial (condition x sex, raw counts) ---\n")

negbin_results <- list()
for (pair in pairs) {
  ga <- pair[1]; gb <- pair[2]; label <- paste(ga, "vs", gb)
  cat(sprintf("  %s\n", label))
  
  rows_list <- vector("list", nrow(df_m))
  n_err     <- 0
  for (i in seq_len(nrow(df_m))) {
    acr <- df_m$acronym[i]; rgn <- df_m$name[i]
    vol <- df_m$volume_mm3[i]
    raw_a <- c(as.numeric(df_m[i, groups_m[[ga]], drop=TRUE]),
               as.numeric(df_f[i, groups_f[[ga]], drop=TRUE]))
    raw_b <- c(as.numeric(df_m[i, groups_m[[gb]], drop=TRUE]),
               as.numeric(df_f[i, groups_f[[gb]], drop=TRUE]))
    raw_a <- raw_a[!is.na(raw_a)]; raw_b <- raw_b[!is.na(raw_b)]
    mean_a <- if (length(raw_a)) mean(raw_a) else NA
    mean_b <- if (length(raw_b)) mean(raw_b) else NA
    
    empty <- data.frame(region_name=rgn, acronym=acr, comparison=paste(ga,"-",gb),
                        negbin_deviance=NA, negbin_df=NA, negbin_p=NA,
                        estimate=NA, std_error=NA, ci_low=NA, ci_high=NA,
                        z_statistic=NA, p_contrast=NA, fold_change=NA, IRR=NA,
                        group_a_mean=mean_a, group_b_mean=mean_b,
                        n_a=length(raw_a), n_b=length(raw_b),
                        n_zero_a=sum(raw_a == 0), n_zero_b=sum(raw_b == 0),
                        direction=NA, stringsAsFactors=FALSE)
    
    if (is.na(vol) || vol <= 0) { rows_list[[i]] <- empty; next }
    
    # BILATERAL SUMS, not hemisphere averages. The composite stores counts
    # AVERAGED ACROSS HEMISPHERES, so ~half of all values end in .5 and are not
    # integers. Feeding those to a count model required rounding, and R's round()
    # is banker's rounding (0.5 -> 0, 2.5 -> 2), which perturbs half the dataset
    # and turns the smallest real observations into structural zeros. Multiplying
    # by 2 recovers the underlying bilateral count exactly (an integer), and the
    # matching offset uses the bilateral volume. This is what makes the
    # "negative binomial GLM on counts" claim in the Methods literally true.
    rows <- list()
    for (g in cond_names) {
      for (col in groups_m[[g]]) {
        cnt <- df_m[i, col]
        if (!is.na(cnt) && cnt >= 0)
          rows[[length(rows)+1]] <- data.frame(count=bilateral_count(cnt),
                                               condition=g, sex="Male",   log_vol=log(2*vol), stringsAsFactors=FALSE)
      }
      for (col in groups_f[[g]]) {
        cnt <- df_f[i, col]
        if (!is.na(cnt) && cnt >= 0)
          rows[[length(rows)+1]] <- data.frame(count=bilateral_count(cnt),
                                               condition=g, sex="Female", log_vol=log(2*vol), stringsAsFactors=FALSE)
      }
    }
    if (length(rows) < 8) { rows_list[[i]] <- empty; next }
    long <- do.call(rbind, rows)
    long$condition <- factor(long$condition, levels=cond_names)
    long$sex       <- factor(long$sex, levels=c("Male","Female"))
    
    tryCatch({
      fit      <- MASS::glm.nb(count ~ condition * sex + offset(log_vol), data=long)
      # The null drops the condition terms but KEEPS sex, so the LRT below tests the
      # CONDITION effect (main + condition:sex) and is comparable to omnibus_p in the
      # ANOVA workbooks. A "~ 1" null would test the whole model (condition + sex +
      # interaction) against nothing, which is a different hypothesis.
      fit_null <- MASS::glm.nb(count ~ sex             + offset(log_vol), data=long)
      lrt      <- anova(fit_null, fit)
      negbin_p  <- lrt[["Pr(Chi)"]][2]
      negbin_dev <- lrt[["LR stat."]][2]
      # df of the LRT, taken from the models rather than the (space-padded)
      # anova column name, so the chi-square test is reportable as X2(df).
      negbin_df <- fit_null$df.residual - fit$df.residual
      
      em   <- emmeans::emmeans(fit, ~ condition)
      con  <- emmeans::contrast(em, method="pairwise", adjust="tukey")
      cdf  <- as.data.frame(con)
      # Position among the levels the MODEL actually fitted -- see the two-way note.
      emm_lv <- as.character(as.data.frame(em)$condition)
      ridx <- pairwise_row(ga, gb, emm_lv)
      crow <- if (!is.na(ridx) && ridx <= nrow(cdf)) cdf[ridx, , drop=FALSE] else cdf[0, ]
      if (nrow(crow) == 0) stop("no matching contrast found")
      
      # Direction-safe log-rate difference from marginal means: B - A
      # (treatment - reference), so IRR = exp(est_log) is a treatment/reference
      # rate ratio and reads the same direction as fold_change.
      emm_df  <- as.data.frame(em)
      l_a     <- emm_df$emmean[emm_df$condition == ga]
      l_b     <- emm_df$emmean[emm_df$condition == gb]
      est_log <- if (length(l_a) && length(l_b)) l_b - l_a else NA
      se_log  <- crow$SE[1]
      z_stat  <- if (!is.na(est_log) && !is.na(se_log) && se_log > 0) est_log/se_log else NA
      p_con   <- crow$p.value[1]
      IRR     <- exp(est_log)
      fc  <- if (!is.na(mean_a) && mean_a > 0) mean_b/mean_a else NA
      dir <- if (!is.na(mean_a) && !is.na(mean_b)) {
        if (mean_a > mean_b) paste(ga,">",gb) else paste(ga,"<",gb)
      } else NA
      
      rows_list[[i]] <- data.frame(region_name=rgn, acronym=acr,
                                   comparison=paste(ga,"-",gb),
                                   negbin_deviance=negbin_dev, negbin_df=negbin_df, negbin_p=negbin_p,
                                   estimate=est_log, std_error=se_log,
                                   ci_low=est_log - 1.959964*se_log, ci_high=est_log + 1.959964*se_log,
                                   z_statistic=z_stat, p_contrast=p_con,
                                   fold_change=fc, IRR=IRR, group_a_mean=mean_a, group_b_mean=mean_b,
                                   n_a=length(raw_a), n_b=length(raw_b),
                                   n_zero_a=sum(raw_a == 0), n_zero_b=sum(raw_b == 0),
                                   direction=dir, stringsAsFactors=FALSE)
    }, error=function(e) {
      rows_list[[i]] <<- empty; n_err <<- n_err + 1
      if (n_err <= 5) message(sprintf("    [%s] %s: %s", label, acr, conditionMessage(e)))
    })
  }
  
  res <- do.call(rbind, rows_list)
  valid_c <- !is.na(res$p_contrast); res$p_fdr <- NA
  res$p_fdr[valid_c] <- p.adjust(res$p_contrast[valid_c], method="BH")
  valid_m <- !is.na(res$negbin_p);   res$negbin_fdr <- NA
  res$negbin_fdr[valid_m] <- p.adjust(res$negbin_p[valid_m], method="BH")
  res$sig_tukey      <- sig_stars(res$p_contrast)
  res$sig_fdr        <- sig_stars(res$p_fdr)
  res$sig_negbin_fdr <- sig_stars(res$negbin_fdr)
  
  negbin_results[[label]] <- res
  cat(sprintf("    sig contrast: %3d (%3d)   FDR: %3d (%3d)%s\n",
              n_below(res$p_contrast, 0.05), n_below(res$p_contrast, 0.10),
              n_below(res$p_fdr,      0.05), n_below(res$p_fdr,      0.10),
              if (n_err) sprintf("   [%d region(s) errored]", n_err) else ""))
}
write_excel(negbin_results,
            file.path(OUTPUT_DIR, "Drug_Statistical_Results_NegBin.xlsx"))

# ==============================================================================
# SUPPLEMENTARY 3: Global-activation sensitivity analysis
# ------------------------------------------------------------------------------
# Density (cells/mm3) remains the primary outcome -- it is what this literature
# uses. But whole-brain total activation itself differs by condition (in males,
# one-way ANOVA p = 0.03), and because density is an ABSOLUTE measure, a global
# difference propagates into every region. A reviewer will ask whether the
# regional findings are selective or a rescaling of that global difference.
#
# This section answers it directly: each animal's regional count is expressed as
# a fraction of that animal's whole-brain total before log10, which removes the
# per-animal global level. Regions that remain significant are activated
# DISPROPORTIONATELY; regions that drop out moved with the whole brain.
#
# Interpret with care in both directions: proportions are compositional (they
# must sum to 1), so this analysis is conservative and cannot itself establish
# absence of a regional effect. It is a robustness check, not a replacement.
# ==============================================================================
cat("\n--- SUPPLEMENTARY: Global-activation sensitivity (proportion of whole-brain total) ---\n")

# Per-animal whole-brain total, over the same 198 analysis regions.
animal_totals <- sapply(unlist(groups_comb, use.names = FALSE),
                        function(cc) sum(df_comb[[cc]], na.rm = TRUE))

norm_glob <- df_comb[, c("name", "acronym")]
for (cc in unlist(groups_comb, use.names = FALSE)) {
  frac <- df_comb[[cc]] / animal_totals[[cc]]
  frac[frac == 0] <- NA
  norm_glob[[cc]] <- log10(frac)
}

globnorm_results <- run_oneway(norm_glob, df_comb, groups_comb, label_suffix = "")
write_excel(globnorm_results,
            file.path(OUTPUT_DIR, "Drug_Statistical_Results_GlobalNorm.xlsx"))

# Side-by-side concordance with the primary density analysis, for the Supplement.
concord <- do.call(rbind, lapply(names(primary_results), function(lab) {
  a <- primary_results[[lab]]; b <- globnorm_results[[lab]]
  data.frame(comparison   = lab,
             n_fdr_density = sum(a$p_fdr < FDR_ALPHA, na.rm = TRUE),
             n_fdr_globnorm= sum(b$p_fdr < FDR_ALPHA, na.rm = TRUE),
             n_fdr_both    = sum(a$p_fdr < FDR_ALPHA & b$p_fdr < FDR_ALPHA, na.rm = TRUE),
             pct_d_positive_density  = round(100*mean(a$cohens_d > 0, na.rm = TRUE), 1),
             pct_d_positive_globnorm = round(100*mean(b$cohens_d > 0, na.rm = TRUE), 1),
             mean_abs_d_density  = round(mean(abs(a$cohens_d), na.rm = TRUE), 3),
             mean_abs_d_globnorm = round(mean(abs(b$cohens_d), na.rm = TRUE), 3),
             sig_tukey = "ns",   # no per-row highlighting in this summary sheet
             stringsAsFactors = FALSE)
}))
print(concord[, setdiff(names(concord), "sig_tukey")], row.names = FALSE)
write_excel(list("Density vs GlobalNorm" = concord),
            file.path(OUTPUT_DIR, "Drug_GlobalNorm_Concordance.xlsx"))

# ==============================================================================
# Model assumption checks (for the Supplementary Methods)
# ------------------------------------------------------------------------------
# Per-region Shapiro-Wilk on the pooled within-group residuals and Levene-style
# (Brown-Forsythe) test of variance homogeneity, so the Supplement can state what
# fraction of regions depart from the ANOVA assumptions rather than leaving a
# reviewer to wonder.
# ==============================================================================
cat("\n--- Model assumption checks ---\n")
assum <- vector("list", nrow(norm_comb))
for (i in seq_len(nrow(norm_comb))) {
  gv <- lapply(cond_names, function(g) {
    v <- as.numeric(norm_comb[i, groups_comb[[g]], drop = TRUE]); v[!is.na(v)]
  })
  names(gv) <- cond_names
  row <- data.frame(region_name = norm_comb$name[i], acronym = norm_comb$acronym[i],
                    shapiro_p = NA_real_, bf_p = NA_real_, stringsAsFactors = FALSE)
  if (all(lengths(gv) >= 3)) {
    resid <- unlist(lapply(gv, function(v) v - mean(v)), use.names = FALSE)
    row$shapiro_p <- tryCatch(stats::shapiro.test(resid)$p.value, error = function(e) NA_real_)
    # Brown-Forsythe: one-way ANOVA on absolute deviations from the group median
    ad  <- unlist(lapply(gv, function(v) abs(v - stats::median(v))), use.names = FALSE)
    grp <- factor(rep(cond_names, lengths(gv)), levels = cond_names)
    row$bf_p <- tryCatch(summary(aov(ad ~ grp))[[1]][["Pr(>F)"]][1], error = function(e) NA_real_)
  }
  assum[[i]] <- row
}
assum_df <- do.call(rbind, assum)
assum_df$sig_shapiro <- sig_stars(assum_df$shapiro_p)
assum_df$sig_bf      <- sig_stars(assum_df$bf_p)
cat(sprintf("  Shapiro-Wilk p<0.05 (non-normal residuals): %d/%d (%.1f%%)\n",
            n_below(assum_df$shapiro_p, 0.05), nrow(assum_df),
            100*mean(assum_df$shapiro_p < 0.05, na.rm = TRUE)))
cat(sprintf("  Brown-Forsythe p<0.05 (unequal variance):   %d/%d (%.1f%%)\n",
            n_below(assum_df$bf_p, 0.05), nrow(assum_df),
            100*mean(assum_df$bf_p < 0.05, na.rm = TRUE)))
write_excel(list("Assumption checks" = assum_df),
            file.path(OUTPUT_DIR, "Model_Assumption_Checks.xlsx"),
            sig_col = "sig_shapiro")

# ---- Save normalized data ----------------------------------------------------
write.csv(norm_comb, file.path(OUTPUT_DIR, "normalized_data_combined.csv"), row.names=FALSE)
write.csv(norm_m,    file.path(OUTPUT_DIR, "normalized_data_males.csv"),    row.names=FALSE)
write.csv(norm_f,    file.path(OUTPUT_DIR, "normalized_data_females.csv"),  row.names=FALSE)

# ---- Summary -----------------------------------------------------------------
cat("\n", strrep("=",65), "\n", sep="")
cat("SUMMARY\n")
cat(strrep("=",65), "\n", sep="")
cat("Normalization: log10(cells/mm3) | ANOVA on pooled males + females\n")
cat("Counts are regions significant at 0.05, with the trend-inclusive count\n")
cat("(0.10) in parentheses. Reported results use the 0.05 counts.\n\n")

cat("PRIMARY - One-way ANOVA (pooled):\n")
for (label in names(primary_results)) {
  res <- primary_results[[label]]
  cat(sprintf("  %-35s  Tukey: %3d (%3d)   FDR: %3d (%3d)\n", label,
              n_below(res$p_tukey, 0.05), n_below(res$p_tukey, 0.10),
              n_below(res$p_fdr,   0.05), n_below(res$p_fdr,   0.10)))
}
cat("\nSECONDARY - Sex-stratified ANOVA:\n")
for (label in names(bysex_results)) {
  res <- bysex_results[[label]]
  cat(sprintf("  %-40s  Tukey: %3d (%3d)   FDR: %3d (%3d)\n", label,
              n_below(res$p_tukey, 0.05), n_below(res$p_tukey, 0.10),
              n_below(res$p_fdr,   0.05), n_below(res$p_fdr,   0.10)))
}
cat("\nSUPPLEMENTARY - Two-way ANOVA (condition x sex):\n")
cat(sprintf("  Sex main effect:                    %3d (%3d) regions\n",
            n_below(sex_effect_df$q_sex, 0.05),
            n_below(sex_effect_df$q_sex, 0.10)))
cat(sprintf("  Condition x Sex interaction:        %3d (%3d) regions\n",
            n_below(int_effect_df$q_interaction, 0.05),
            n_below(int_effect_df$q_interaction, 0.10)))
for (label in setdiff(names(twoway_results),
                      c("Sex Main Effect","Condition x Sex Interact"))) {
  res <- twoway_results[[label]]
  cat(sprintf("  %-35s  Tukey: %3d (%3d)   FDR: %3d (%3d)\n", label,
              n_below(res$p_tukey, 0.05), n_below(res$p_tukey, 0.10),
              n_below(res$p_fdr,   0.05), n_below(res$p_fdr,   0.10)))
}
cat("\nOutputs:\n")
cat("  Drug_Statistical_Results_Primary.xlsx  [PRIMARY]\n")
cat("  Drug_Statistical_Results_BySex.xlsx    [SECONDARY]\n")
cat("  Drug_Statistical_Results_TwoWay.xlsx   [SUPPLEMENTARY]\n")
cat("  Drug_Statistical_Results_NegBin.xlsx   [SUPPLEMENTARY]\n")
cat("  Drug_Statistical_Results_GlobalNorm.xlsx [SUPPLEMENTARY - sensitivity]\n")
cat("  Drug_GlobalNorm_Concordance.xlsx       [SUPPLEMENTARY - sensitivity]\n")
cat("  Model_Assumption_Checks.xlsx           [SUPPLEMENTARY]\n")
cat("  normalized_data_combined/males/females.csv\n")
# ---- Provenance --------------------------------------------------------------
# Written on every run so the deposited results always carry the exact package
# versions that produced them (referenced by the header and the Data Availability
# statement).
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "sessionInfo.txt"))
cat("  Wrote sessionInfo.txt\n")

cat("\nDone.\n")
