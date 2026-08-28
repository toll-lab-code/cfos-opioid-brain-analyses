# preprocess_neuroinfo.R
# =============================================================================
# NeuroInfo export cleanup: hemisphere averaging and region filtering.
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
#   Rscript scripts/preprocess_neuroinfo.R
# Input and output locations are set in the Config section below; see README.
# =============================================================================
#
# This is the first stage of the pipeline. It takes per-animal NeuroInfo region
# tables exported as CSV, averages the left and right hemisphere rows for each
# region, restricts the result to the analysis region set, and writes one
# "<animal>_counts.csv" per animal. Those per-animal files are then compiled
# into Drug_Composite.xlsx, which is the input to regional_analysis.R and
# network_analysis.R.
#
# Every file under RAW_DIR matching FILE_PATTERN is processed, searched
# recursively, so all conditions and both sexes are handled in a single run.
# No particular input file-naming scheme is assumed; see the Config section.
# =============================================================================

# ---- Packages ----------------------------------------------------------------
required_pkgs <- c("dplyr", "jsonlite", "purrr", "tibble")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages. Please run:\n  install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}
library(dplyr); library(jsonlite); library(tools); library(purrr); library(tibble)

# ---- Config ------------------------------------------------------------------
# Paths are relative to the repository root. Run from the repo root, e.g.:
#   Rscript scripts/preprocess_neuroinfo.R
# Override by editing the directories below or by setting the environment
# variables CFOS_DATA_DIR / CFOS_OUTPUT_DIR before running.
DATA_DIR   <- Sys.getenv("CFOS_DATA_DIR",   "data")
OUTPUT_DIR <- Sys.getenv("CFOS_OUTPUT_DIR", "results")

# Raw exports, searched recursively so any subfolder structure (per condition,
# per sex, per cohort) is picked up in a single run. The default pattern matches
# every .csv under RAW_DIR, so no particular file-naming scheme is assumed.
# Narrow it if the raw folder also holds CSVs that should not be processed, e.g.
#   FILE_PATTERN <- " - Raw\\.csv$"
RAW_DIR      <- file.path(DATA_DIR, "neuroinfo_raw")
CLEANED_DIR  <- file.path(OUTPUT_DIR, "neuroinfo_cleaned")
FILE_PATTERN <- "\\.csv$"

# Output naming. If the input file name ends in any of RAW_SUFFIXES, that suffix
# is stripped; CLEAN_SUFFIX is then appended. So "Veh_F06 - Raw.csv" and
# "Veh_F06.csv" both yield "Veh_F06_counts.csv". Set RAW_SUFFIXES to
# character(0) to append without stripping anything.
RAW_SUFFIXES <- c(" - Raw", "_raw", "_Raw")
CLEAN_SUFFIX <- "_counts"

# Export layout. Defaults match a standard NeuroInfo region-table export.
# Adjust these if the export has a different column arrangement.
HEADER_SKIP <- 1                       # lines above the header row
DROP_COLS   <- 7:9                     # columns dropped before averaging (NULL = drop none)
FLUOR_COL   <- 5                       # column of the averaged table renamed below (NA = skip)
FLUOR_NAME  <- "fluorescence strength"
BLANK_COL   <- 6                       # column of the averaged table cleared (NA = skip)

# Allen CCFv3 structure ontology, used only to fill in region names where the
# export leaves the name field blank. It does not affect which regions are kept.
# Download (structure graph 1), saved with its default file name:
#   http://api.brain-map.org/api/v2/structure_graph_download/1.json
allen_json <- file.path(DATA_DIR, "1.json")

# Regions excluded from analysis after hemisphere averaging. Study-specific:
# edit for a different analysis.
remove_regions <- c("BA", "ARH", "PVp", "VMPO","VLPO", "RCH", "ME", "NTB", "PPY",
                    "SCO", "IO", "RPA", "SCH", "TU", "OV", "VII", "LRN", "SO", "AMB", 
                    "MARN", "PMv", "ACVII", "AI", "LIN", "RM", "RO", "PG", "SOC",
                    "NLOT", "AVP", "TRN", "VMH")

# Analysis region set (CCFv3 acronyms). A region is carried through only if it is
# present in the export AND has both a left and a right hemisphere row, so the
# number of regions in the cleaned output is smaller than the length of this
# vector: regions absent or unilateral in a given export drop out at this step.
# The 198 grey-matter regions analysed in the paper are the regions surviving
# both this whitelist and the remove_regions exclusion above.
brain_regions <- c(
  "FRP","MO","SS","GU","VISC","AUD","VIS","ACA","PL","ILA","ORB","AI","RSP","PTLp","TEa","PERI","ECT",
  "MOB","AOB","AON","TT","DP","PIR","NLOT","COA","PAA","TR",
  "CA","DG","FC","IG",
  "ENT","PAR","POST","PRE","SUB","ProS","HATA","APr",
  "CLA","EP","LA","BLA","BMA","PA",
  "CP",
  "ACB","FS","OT","LSS",
  "LS","SF","SH",
  "AAA","BA","CEA","IA","MEA",
  "GPe","GPi",
  "SI","MA",
  "MS","NDB","TRS",
  "BST","BAC",
  "VENT","SPF","SPA","PP","GENd",
  "LAT","ATN","MED","MTN","ILM","RT","GENv","EPI",
  "SO","ASO","PVH","PVa","PVi","ARH",
  "ADP","AHA","AVP","AVPV","DMH","MEPO","MPO","OV","PD","PS","PSCH","PVp","PVpo","SBPV","SCH","SFO","VMPO","VLPO",
  "AHN","MBO","MPN","PMd","PMv","PVHd","VMH","PH",
  "LHA","LPO","PST","PSTN","PeF","RCH","STN","TU","ZI",
  "ME",
  "SCs","IC","NB","SAG","PBG","MEV","SCO",
  "SNr","VTA","PN","RR","MRN","SCm","PAG","PRT","InCo","CUN","RN","III","MA3","EW","IV","Pa4","VTN","AT","LT","DT","MT","SNl",
  "SNc","PPN","RAmb",
  "NLL","PSV","PB","SOC",
  "B","DTN","LTN","PDTg","PCG","PG","PRNc","PRNv","SG","SSN","SUT","TRN","V","P5","Acs5","PC5","I5",
  "CS","LC","LDT","NI","PRNr","RPO","SLC","SLD",
  "AP","CN","DCN","ECU","NTB","NTS","SPVC","SPVI","SPVO","Pa5","z",
  "VI","ACVI","VII","ACVII","EV","AMB","DMX","ECO","GRN","ICB","IO","IRN","ISN","LIN","LRN","MARN","MDRN","PARN","PAS","PGRN","PHY","PMR","PPY","VNC","x","XII","y","INV",
  "RM","RPA","RO",
  "LING","CENT2","CENT3","CUL4","CUL5","DEC","FOTU","PYR","UVU","NOD",
  "SIM","ANcr1","ANcr2","PRM","COPY","PFL","FL",
  "FN","IP","DN","VeCB"
)

# ---- Allen ontology lookup (loaded once) -------------------------------------

stopifnot(file.exists(allen_json))
j <- jsonlite::fromJSON(allen_json, simplifyVector = FALSE)

# Walks the ontology and returns every acronym/name pair. The Allen download is
# an API response object whose structure tree sits under "msg", so the walker
# recurses into any list-valued element rather than assuming a fixed shape; this
# also makes it tolerant of ontology files saved with a different wrapper.
# Each node is visited once and rows are combined in a single rbind.
collect_nodes <- function(x) {
  if (!is.list(x)) return(NULL)
  rows <- list()
  if (!is.null(x$acronym) && !is.null(x$name)) {
    rows[[length(rows) + 1L]] <- data.frame(acronym = as.character(x$acronym),
                                            name    = as.character(x$name),
                                            stringsAsFactors = FALSE)
  }
  for (el in x) {
    if (is.list(el)) rows[[length(rows) + 1L]] <- collect_nodes(el)
  }
  if (length(rows) == 0L) {
    return(data.frame(acronym = character(), name = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

lookup <- collect_nodes(j) %>%
  filter(!is.na(acronym), acronym != "", !is.na(name), name != "") %>%
  distinct(acronym, .keep_all = TRUE)

# ---- Core cleaning function --------------------------------------------------

process_one_file <- function(in_csv, raw_dir, cleaned_dir, lookup, brain_regions, remove_regions) {
  
  # Output path mirrors the raw file's location under RAW_DIR, so any
  # subfolder structure is preserved in the output.
  base     <- file_path_sans_ext(basename(in_csv))
  raw_norm <- normalizePath(raw_dir, winslash = "/", mustWork = FALSE)
  in_norm  <- normalizePath(in_csv,  winslash = "/", mustWork = FALSE)
  rel_path <- substr(in_norm, nchar(raw_norm) + 2L, nchar(in_norm))
  rel_dir  <- dirname(rel_path)
  out_dir  <- if (rel_dir == ".") cleaned_dir else file.path(cleaned_dir, rel_dir)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  # Strip any configured raw-file suffix, then append CLEAN_SUFFIX. Matching is
  # literal, so suffixes containing regex characters are safe.
  out_base <- base
  for (sfx in RAW_SUFFIXES) {
    if (nzchar(sfx) && endsWith(out_base, sfx)) {
      out_base <- substr(out_base, 1L, nchar(out_base) - nchar(sfx))
      break
    }
  }
  out_csv  <- file.path(out_dir, paste0(out_base, CLEAN_SUFFIX, ".csv"))
  
  message("\nProcessing: ", basename(in_csv))
  message("Output:     ", basename(out_csv))
  
  # --- Load CSV ---
  rawdata <- read.csv(in_csv, skip = HEADER_SKIP, check.names = FALSE)
  
  # --- Make headers unique ---
  names(rawdata) <- trimws(names(rawdata))
  names(rawdata) <- make.unique(names(rawdata), sep = "_dup")
  
  # Drop the export's bookkeeping columns (DROP_COLS; positions 7:9, spreadsheet
  # columns G:I, in a standard NeuroInfo export). Positions beyond the width of
  # this file are ignored.
  if (!is.null(DROP_COLS)) {
    drop_idx <- intersect(DROP_COLS, seq_along(rawdata))
    if (length(drop_idx) > 0) rawdata <- rawdata[, -drop_idx, drop = FALSE]
  }
  
  # --- Drop any 'section' column(s) ---
  if ("section" %in% names(rawdata)) {
    rawdata <- rawdata[, names(rawdata) != "section", drop = FALSE]
  }
  
  # Ensure expected columns exist
  if (!"acronym" %in% names(rawdata)) stop("Column 'acronym' not found in CSV: ", basename(in_csv))
  if (!"name" %in% names(rawdata)) rawdata$name <- NA_character_
  if (!"hemisphere" %in% names(rawdata)) rawdata$hemisphere <- NA_character_
  
  # Clean key fields. Hemisphere labels are lower-cased and mapped onto
  # "left"/"right"; common single-letter and abbreviated variants are accepted
  # so that exports labelled L/R or lh/rh work without modification.
  rawdata$acronym    <- trimws(rawdata$acronym)
  hemi <- tolower(trimws(as.character(rawdata$hemisphere)))
  hemi[hemi %in% c("l", "lh", "left")]  <- "left"
  hemi[hemi %in% c("r", "rh", "right")] <- "right"
  rawdata$hemisphere <- hemi
  
  # Fill missing names using Allen lookup
  fill_names <- lookup$name[match(rawdata$acronym, lookup$acronym)]
  rawdata$name <- ifelse(is.na(rawdata$name) | rawdata$name == "", fill_names, rawdata$name)
  
  # Intersect whitelist with what is present
  valid_regions <- brain_regions[brain_regions %in% rawdata$acronym]
  missing_regions <- setdiff(brain_regions, valid_regions)
  if (length(missing_regions) > 0) {
    message("Regions not found in data: ", paste(missing_regions, collapse = ", "))
  }
  
  # Averaging logic
  skipped_info <- list()
  
  average_region <- function(region_data) {
    left_row  <- region_data %>% filter(hemisphere == "left")
    right_row <- region_data %>% filter(hemisphere == "right")
    
    if (nrow(left_row) == 0 | nrow(right_row) == 0) {
      skipped_info[[region_data$acronym[1]]] <<- paste("Left rows:", nrow(left_row),
                                                       "Right rows:", nrow(right_row))
      return(NULL)
    }
    
    avg_row <- tibble::tibble(
      acronym = region_data$acronym[1],
      name    = region_data$name[1]
    )
    
    exclude <- c("hemisphere","acronym","name")
    candidate_cols <- setdiff(names(region_data), exclude)
    
    for (col in candidate_cols) {
      vals <- suppressWarnings(as.numeric(c(left_row[[col]], right_row[[col]])))
      if (length(vals) == 2 && all(!is.na(vals))) {
        avg_row[[col]] <- mean(vals)
      }
    }
    avg_row
  }
  
  results_list <- lapply(valid_regions, function(region) {
    region_data <- rawdata %>% filter(acronym == region)
    average_region(region_data)
  })
  
  new_df <- dplyr::bind_rows(results_list)
  
  # Label the fluorescence measure and clear the placeholder column that the
  # downstream compilation step expects. Positions are configurable (FLUOR_COL /
  # BLANK_COL); set either to NA to skip that step.
  if (!is.na(FLUOR_COL) && ncol(new_df) >= FLUOR_COL) colnames(new_df)[FLUOR_COL] <- FLUOR_NAME
  if (!is.na(BLANK_COL) && ncol(new_df) >= BLANK_COL) new_df[[BLANK_COL]] <- ""
  
  # Remove specified regions
  new_df <- new_df %>% filter(!acronym %in% remove_regions)
  
  # Save cleaned CSV
  write.csv(new_df, out_csv, row.names = FALSE)
  
  # Return summary info (for batch report). region_signature is the ordered
  # acronym vector, used by the cross-animal consistency guard after the batch.
  tibble::tibble(
    input_file = basename(in_csv),
    output_file = basename(out_csv),
    n_rows_out = nrow(new_df),
    region_signature = paste(new_df$acronym, collapse = "|"),
    skipped_regions = paste(names(skipped_info), collapse = ", ")
  )
}

# ---- Batch run ---------------------------------------------------------------
# Recursive search picks up every raw export under RAW_DIR in one pass, across
# all conditions and both sexes.
dir.create(CLEANED_DIR, showWarnings = FALSE, recursive = TRUE)

files <- list.files(RAW_DIR, pattern = FILE_PATTERN, full.names = TRUE, recursive = TRUE)
if (length(files) == 0) {
  stop("No files matching '", FILE_PATTERN, "' found under: ", RAW_DIR)
}
message("Found ", length(files), " raw file(s) under ", RAW_DIR)

batch_summary <- purrr::map_dfr(
  files,
  ~ tryCatch(
    process_one_file(.x, RAW_DIR, CLEANED_DIR, lookup, brain_regions, remove_regions),
    error = function(e) tibble::tibble(input_file = basename(.x), error = e$message)
  )
)

# Batch report: one row per input file, including any that errored.
summary_csv <- file.path(CLEANED_DIR, "batch_cleaning_summary.csv")
write.csv(batch_summary, summary_csv, row.names = FALSE)

# ---- Cross-animal consistency guard ------------------------------------------
# Every animal must end up with the SAME regions in the SAME order. A region is
# kept only when it has both a left and a right hemisphere row, so an export that
# is missing a region, or has it unilaterally, silently yields one row fewer for
# that animal. The next step (compiling Drug_Composite.xlsx) lays the per-animal
# columns side by side in a fixed grid, so a single animal with a different row
# set would shift every value below it and pair one region's counts with another
# region's -- with no error anywhere downstream. Fail here instead.
if ("error" %in% names(batch_summary) && any(!is.na(batch_summary$error))) {
  failed <- batch_summary[!is.na(batch_summary$error), ]
  stop(sprintf("%d file(s) failed to process:\n  %s",
               nrow(failed),
               paste(sprintf("%s: %s", failed$input_file, failed$error), collapse = "\n  ")))
}
sigs <- unique(batch_summary$region_signature)
if (length(sigs) != 1) {
  modal   <- names(sort(table(batch_summary$region_signature), decreasing = TRUE))[1]
  ref     <- strsplit(modal, "|", fixed = TRUE)[[1]]
  offend  <- batch_summary[batch_summary$region_signature != modal, ]
  detail  <- vapply(seq_len(nrow(offend)), function(i) {
    got  <- strsplit(offend$region_signature[i], "|", fixed = TRUE)[[1]]
    miss <- setdiff(ref, got); extra <- setdiff(got, ref)
    sprintf("  %s (%d regions)\n      missing: %s\n      extra:   %s",
            offend$input_file[i], length(got),
            if (length(miss))  paste(miss,  collapse = ", ") else "none",
            if (length(extra)) paste(extra, collapse = ", ") else "none")
  }, character(1))
  stop(sprintf(paste0(
    "cross-animal consistency check FAILED.\n",
    "  %d of %d animals do not share the modal region set (%d regions).\n",
    "  Compiling these into a fixed-grid composite would misalign regions.\n%s\n",
    "  Check the raw export for the listed animal(s): a region kept for the others\n",
    "  is absent or present in only one hemisphere here."),
    nrow(offend), nrow(batch_summary), length(ref),
    paste(detail, collapse = "\n")))
}
message("\nConsistency OK: all ", nrow(batch_summary), " animals share the same ",
        length(strsplit(sigs, "|", fixed = TRUE)[[1]]), " regions in the same order.")

# ---- Provenance --------------------------------------------------------------
# Written to a preprocess-specific filename so it does not overwrite the
# sessionInfo.txt that regional_analysis.R writes into the same results folder.
# The two stages load different packages, so both records are worth keeping.
writeLines(capture.output(sessionInfo()),
           file.path(CLEANED_DIR, "sessionInfo_preprocess.txt"))

message("\nBatch complete. ", length(files), " file(s) processed.")
message("Cleaned files: ", CLEANED_DIR)
message("Summary saved to: ", summary_csv)
message("Session info:  ", file.path(CLEANED_DIR, "sessionInfo_preprocess.txt"))

