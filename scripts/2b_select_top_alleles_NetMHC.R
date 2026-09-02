#!/usr/bin/env Rscript
# ============================================================
# JUGNU_IMMUNOPEPTIDOME — Step 2b: Select top NetMHCpan allele per peptide
# July 2026 | Gaurav Raichand | The Institute of Cancer Research
#
# Input:  canonical_0Xmer_netmhcpan_YYYY_MMDD.tsv  (from Step 2a, x4)
# Output: canonical_0Xmer_netmhcpan_selected_alleles_YYYYMMDD.tsv  (x4)
#
# Logic: keep strong binder rows (binder == "SB"), then for each peptide
# keep the allele with the lowest EL rank (best presenter).
# ============================================================

suppressPackageStartupMessages(library(data.table))

input_dir    <- Sys.getenv("OUTPUT_DIR")
directory_out <- Sys.getenv("STEP13_OUTPUT_DIR")

if (nchar(input_dir)     == 0) stop("OUTPUT_DIR not set -- source config.sh first")
if (nchar(directory_out) == 0) directory_out <- input_dir
dir.create(directory_out, showWarnings = FALSE, recursive = TRUE)

# Auto-detect run_date from newest canonical_09mer_netmhcpan_*.tsv
nmp_files <- list.files(input_dir,
                         pattern = "^canonical_09mer_netmhcpan_[0-9]{4}_[0-9]{4}\\.tsv$")
if (length(nmp_files) == 0)
  stop("[ERROR] No canonical_09mer_netmhcpan_YYYY_MMDD.tsv in: ", input_dir,
       "\n  Has Step 2a finished?")

run_date <- sub("^canonical_09mer_netmhcpan_(.+)\\.tsv$", "\\1",
                nmp_files[which.max(file.info(
                  file.path(input_dir, nmp_files))$mtime)])
out_date <- gsub("_", "", run_date)

cat("[INFO] Reading 2a files from:    ", input_dir,     "\n")
cat("[INFO] Writing 2b files to:      ", directory_out, "\n")
cat("[INFO] Detected 2a run_date:     ", run_date,      "\n")
cat("[INFO] Output files dated:       ", out_date,      "\n\n")

###########################################################################
# Helper: load, filter SB, select best allele per peptide
###########################################################################

select_top_allele <- function(length_label) {

  infile  <- file.path(input_dir,
               paste0("canonical_", length_label, "mer_netmhcpan_",
                      run_date, ".tsv"))
  outfile <- file.path(directory_out,
               paste0("canonical_", length_label,
                      "mer_netmhcpan_selected_alleles_", out_date, ".tsv"))

  if (!file.exists(infile)) {
    warning("[WARN] File not found, skipping: ", infile)
    return(invisible(NULL))
  }

  # Stream-filter the complete Step 2a table before R reads it. With 70
  # alleles the full TSV is tens of GB; fread()-then-filter duplicates that
  # table during sorting and can exceed the Slurm memory limit. Step 2a keeps
  # all SB/WB/NB rows, while only SB rows enter the in-memory selection here.
  awk_program <- paste(
    'BEGIN { FS="\\t" }',
    'NR == 1 || $8 == "SB" || $8 == "<= SB"'
  )
  awk_command <- sprintf("awk %s %s", shQuote(awk_program), shQuote(infile))
  sb <- fread(cmd = awk_command, na.strings = c("", "NA"))
  cat(sprintf("[INFO] %smer: %d SB rows streamed from full Step 2a TSV\n",
              length_label, nrow(sb)))

  sb[, netmhcpan_EL_rank  := as.numeric(netmhcpan_EL_rank)]
  sb[, netmhcpan_EL_score := as.numeric(netmhcpan_EL_score)]

  if (nrow(sb) == 0) {
    warning("[WARN] No strong binders for ", length_label, "mer")
    fwrite(sb, outfile, sep = "\t", na = "NA", col.names = TRUE, quote = FALSE)
    return(invisible(sb))
  }

  # Best allele per peptide: lowest EL rank, ties broken by highest EL score
  setorder(sb, peptide, netmhcpan_EL_rank, -netmhcpan_EL_score)
  best <- sb[, .SD[1], by = peptide]
  cat(sprintf("[INFO] %smer: %d unique peptides after best-allele selection\n",
              length_label, nrow(best)))

  fwrite(best, outfile, sep = "\t", na = "NA", col.names = TRUE, quote = FALSE)
  cat(sprintf("[INFO] Written: %s\n\n", outfile))
  invisible(best)
}

for (len in c("08", "09", "10", "11")) {
  select_top_allele(len)
}

# Summary
cat("\n=== 2b Summary ===\n")
for (len in c("08", "09", "10", "11")) {
  f <- file.path(directory_out,
                 paste0("canonical_", len,
                        "mer_netmhcpan_selected_alleles_", out_date, ".tsv"))
  if (file.exists(f)) {
    cat(sprintf("  %smer: %d peptides -> %s\n",
                len, nrow(fread(f, select = 1L)), basename(f)))
  } else {
    cat(sprintf("  %smer: MISSING\n", len))
  }
}
cat("\n[DONE] Step 2b complete.\n")
