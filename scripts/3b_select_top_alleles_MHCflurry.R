#!/usr/bin/env Rscript
# ============================================================
# JUGNU_IMMUNOPEPTIDOME — Step 3b: Select top MHCflurry allele per peptide
# July 2026 | Gaurav Raichand | The Institute of Cancer Research
#
# Input:  canonical_0Xmers_flank_mhcflurry_YYYYMMDD.csv  (from Step 3a, x4)
# Output: canonical_0Xmer_mhcflurry_top_YYYYMMDD.tsv  (x4)
#
# Logic: for each peptide, keep the allele with the highest
# mhcflurry_presentation_score (cohort-wide best -- no per-sample logic,
# no OptiType, no samples.txt).
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

input_dir    <- Sys.getenv("OUTPUT_DIR")
directory_out <- Sys.getenv("STEP14_OUTPUT_DIR")

if (nchar(input_dir)     == 0) stop("OUTPUT_DIR not set -- source config.sh first")
if (nchar(directory_out) == 0) directory_out <- input_dir
dir.create(directory_out, showWarnings = FALSE, recursive = TRUE)

current_date <- format(Sys.Date(), "%Y%m%d")

# Auto-detect newest MHCflurry output file
latest_file <- function(pattern) {
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("[ERROR] No file matching: ", pattern,
                                " in ", input_dir)
  files[which.max(file.info(files)$mtime)]
}

# Read both normal MHCflurry CSVs and the Step 3a CSV written with an
# 8-column header above 11-column data rows. In the malformed legacy file,
# affinity is field 7 and presentation score is field 10.
read_mhcflurry <- function(path) {
  probe <- readLines(path, n = 2L, warn = FALSE)
  if (length(probe) < 2L) stop("[ERROR] Empty MHCflurry file: ", path)

  header_n <- lengths(strsplit(probe[1L], ",", fixed = TRUE))
  data_n   <- lengths(strsplit(probe[2L], ",", fixed = TRUE))

  if (header_n == 8L && data_n == 11L) {
    cat("[WARN] Detected legacy 8-header/11-data CSV; parsing fields 1,2,7,10\n")
    return(fread(
      path, skip = 1L, header = FALSE,
      select = c(1L, 2L, 7L, 10L),
      col.names = c("peptide", "allele", "mhcflurry_affinity",
                    "mhcflurry_presentation_score"),
      na.strings = c("", "NA")
    ))
  }

  header <- names(fread(path, nrows = 0L))
  affinity <- intersect(
    c("mhcflurry_affinity", "mhcflurry_binding_affinity"), header
  )[1L]
  required <- c("peptide", "allele", "mhcflurry_presentation_score")
  missing <- setdiff(required, header)
  if (length(missing) || is.na(affinity)) {
    stop("[ERROR] Unsupported MHCflurry columns in ", path,
         "; header: ", paste(header, collapse = ","))
  }

  dt <- fread(
    path, select = c("peptide", "allele", affinity,
                     "mhcflurry_presentation_score"),
    na.strings = c("", "NA")
  )
  if (affinity != "mhcflurry_affinity") {
    setnames(dt, affinity, "mhcflurry_affinity")
  }
  dt
}

###########################################################################
# Process one length at a time to control memory (files are large)
###########################################################################

for (len in c("08", "09", "10", "11")) {
  cat(sprintf("\n--- %smer ---\n", len))

  files <- list.files(
    input_dir,
    pattern = sprintf("^canonical_%smers_flank_mhcflurry_[0-9_]+\\.csv$", len),
    full.names = TRUE
  )
  if (length(files) == 0L) {
    cat(sprintf("[WARN] No MHCflurry output for %smer -- skipping\n", len))
    next
  }
  f <- files[which.max(file.info(files)$mtime)]
  cat(sprintf("[INFO] Loading: %s\n", basename(f)))

  dt <- read_mhcflurry(f)
  dt[, mhcflurry_presentation_score := as.numeric(mhcflurry_presentation_score)]
  dt[, mhcflurry_affinity           := as.numeric(mhcflurry_affinity)]
  if (all(is.na(dt$mhcflurry_presentation_score)) ||
      all(is.na(dt$mhcflurry_affinity))) {
    stop("[ERROR] Parsed MHCflurry score columns contain only NA: ", f)
  }
  cat(sprintf("[INFO] %d rows loaded\n", nrow(dt)))

  # Best allele per peptide: highest presentation score
  setorder(dt, peptide, -mhcflurry_presentation_score)
  top <- dt[, .SD[1L], by = peptide]
  cat(sprintf("[INFO] %d unique peptides after best-allele selection\n", nrow(top)))

  out_file <- file.path(directory_out,
    paste0("canonical_", len, "mer_mhcflurry_top_", current_date, ".tsv"))
  fwrite(top, out_file, sep = "\t", col.names = TRUE, quote = FALSE)
  cat(sprintf("[INFO] Written: %s\n", out_file))

  rm(dt, top); gc()
}

cat("\n[DONE] Step 3b complete.\n")
