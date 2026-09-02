#!/usr/bin/env Rscript
# ============================================================
# JUGNU_IMMUNOPEPTIDOME — Step 1: Gen_nmers.R
# July 2026 | Gaurav Raichand | The Institute of Cancer Research
#
# Purpose: Slide a window across every human protein-coding transcript
#          sequence (from EnsDb) to extract 8-11mers with guaranteed
#          accurate ENST-relative amino acid coordinates.
#
# Why EnsDb instead of UniProt FASTA:
#   Using UniProt sequences and then joining ENST IDs risks coordinate
#   mismatches -- UniProt and Ensembl protein sequences can differ in
#   length due to isoform differences and annotation choices. By pulling
#   sequences directly from EnsDb per ENST, aa_start/aa_end are always
#   correct relative to the ENST protein they came from.
#
# Input  : EnsDb.Hsapiens.v86 (built-in Bioconductor package)
#          No external FASTA needed for tiling.
#          FASTA_CANONICAL is still used by downstream steps (2a, 3a)
#          but not by this script.
#
# Output (per length, all in OUTPUT_DIR):
#   canonical_0Xmers_YYYYMMDD.tsv
#     Columns: enst | aa_start | aa_end | n_mer | ctex_up | ctex_dn
#     aa_start/aa_end: 1-based, inclusive, relative to the ENST protein
#
#   canonical_0Xmer_netmhcpan_input_YYYYMMDD.pep
#     Unique peptide sequences, one per line (for Step 2a NetMHCpan)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ensembldb)
  library(EnsDb.Hsapiens.v86)
})

start_time <- proc.time()

###########################################################################
# 1. Config
###########################################################################

output_dir <- Sys.getenv("OUTPUT_DIR")
if (output_dir == "") stop("[ERROR] OUTPUT_DIR not set -- source config.sh first.")

run_date <- format(Sys.Date(), "%Y%m%d")

cat("[CONFIG] OUTPUT_DIR:      ", output_dir, "\n")
cat("[CONFIG] HLA_ALLELE_FILE: ", Sys.getenv("HLA_ALLELE_FILE"), "\n")
cat("[CONFIG] Run date:        ", run_date, "\n\n")

###########################################################################
# 2. Pull all protein-coding transcript sequences from EnsDb
#
# proteins() returns one row per transcript with its full translated
# protein sequence. We filter to protein_coding biotype only.
# aa_start/aa_end will be positions in THIS sequence -- guaranteed correct.
###########################################################################

cat("[INFO] Loading protein sequences from EnsDb.Hsapiens.v86...\n")

prot_df <- tryCatch({
  as.data.table(
    proteins(
      EnsDb.Hsapiens.v86,
      columns = c("tx_id", "protein_sequence"),
      filter  = TxBiotypeFilter("protein_coding")
    )
  )
}, error = function(e) {
  stop("[ERROR] Failed to load proteins from EnsDb: ", conditionMessage(e))
})

setnames(prot_df, "tx_id", "enst")

# Drop any rows with missing or empty protein sequence
prot_df <- prot_df[!is.na(protein_sequence) & nchar(protein_sequence) > 0]

# Quality filter: remove sequences with stop codons or ambiguous AAs
before <- nrow(prot_df)
prot_df <- prot_df[!grepl("[*BJOUXZ]", protein_sequence)]
cat(sprintf("[INFO] Transcripts loaded: %d | After quality filter: %d (removed %d)\n\n",
            before, nrow(prot_df), before - nrow(prot_df)))

if (nrow(prot_df) == 0) stop("[ERROR] No protein sequences survived quality filter.")

enst_vec <- prot_df$enst
seq_vec  <- prot_df$protein_sequence
n_prot   <- length(seq_vec)
cat(sprintf("[INFO] %d protein-coding transcripts to tile\n\n", n_prot))

###########################################################################
# 3. Sliding window tiling per length
#
# For each transcript and each window start position:
#   enst     = Ensembl transcript ID
#   aa_start = 1-based start position in the ENST protein sequence
#   aa_end   = 1-based end position (inclusive)
#   n_mer    = peptide sequence
#   n_flank  = up to 30 AA N-terminal context (raw, before padding)
#   c_flank  = up to 30 AA C-terminal context (raw, before padding)
#
# Flanks are padded to exactly 30 chars with "-" (MHCflurry convention).
# All positions guaranteed correct relative to the ENST protein.
###########################################################################

for (h in 8:11) {
  lenpad <- sprintf("%02d", h)

  coord_file  <- file.path(output_dir,
                   sprintf("canonical_%smers_%s.tsv",               lenpad, run_date))
  netmhc_file <- file.path(output_dir,
                   sprintf("canonical_%smer_netmhcpan_input_%s.pep", lenpad, run_date))

  # Find most recent existing outputs for this length (any date)
  existing_coord  <- sort(Sys.glob(file.path(output_dir,
                      sprintf("canonical_%smers_*.tsv", lenpad))), decreasing = TRUE)
  existing_netmhc <- sort(Sys.glob(file.path(output_dir,
                      sprintf("canonical_%smer_netmhcpan_input_*.pep", lenpad))), decreasing = TRUE)

  if (length(existing_coord) > 0 && length(existing_netmhc) > 0 &&
      file.size(existing_coord[1])  > 100 &&
      file.size(existing_netmhc[1]) > 100) {
    cat(sprintf("[INFO] %smer: recent outputs found -- skipping\n", lenpad))
    cat(sprintf("[INFO]   coord:  %s\n", basename(existing_coord[1])))
    cat(sprintf("[INFO]   netmhc: %s\n", basename(existing_netmhc[1])))
    next
  }

  cat(sprintf("[INFO] === %smer: tiling %d transcripts ===\n", lenpad, n_prot))

  tiles <- vector("list", n_prot)

  for (i in seq_len(n_prot)) {
    seq <- seq_vec[i]
    len <- nchar(seq)
    if (len < h) next

    starts <- seq_len(len - h + 1L)
    tiles[[i]] <- data.table(
      enst     = enst_vec[i],
      aa_start = starts,
      aa_end   = starts + h - 1L,
      n_mer    = substring(seq, starts, starts + h - 1L),
      n_flank  = substring(seq, pmax(1L, starts - 30L), starts - 1L),
      c_flank  = substring(seq, starts + h, pmin(len, starts + h + 29L))
    )
  }

  cat(sprintf("[INFO] %smer: combining tiles...\n", lenpad))
  dt <- rbindlist(tiles, use.names = FALSE, fill = FALSE)
  rm(tiles); gc()

  # Drop non-standard amino acids
  dt <- dt[!grepl("[^ACDEFGHIKLMNPQRSTVWY]", n_mer)]

  # Pad flanks to exactly 30 characters (MHCflurry convention)
  dt[, ctex_up := str_pad(n_flank, width = 30L, side = "left",  pad = "-")]
  dt[, ctex_dn := str_pad(c_flank, width = 30L, side = "right", pad = "-")]

  cat(sprintf("[INFO] %smer: %d total peptide-position rows across %d transcripts\n",
              lenpad, nrow(dt), uniqueN(dt$enst)))

  # --- Output 1: Coordinate table ---------------------------------------
  fwrite(dt[, .(enst, aa_start, aa_end, n_mer, ctex_up, ctex_dn)],
         coord_file, sep = "\t", col.names = TRUE, quote = FALSE)
  cat(sprintf("[INFO] Coord TSV written: %s\n", basename(coord_file)))

  # --- Output 2: NetMHCpan input ----------------------------------------
  # Unique peptide sequences only -- protein origin irrelevant for scoring
  unique_peps <- unique(dt$n_mer)
  writeLines(unique_peps, netmhc_file)
  cat(sprintf("[INFO] NetMHCpan .pep written: %s (%d unique peptides)\n",
              basename(netmhc_file), length(unique_peps)))

  rm(dt, unique_peps); gc()
  cat(sprintf("[INFO] %smer done.\n\n", lenpad))
}

runtime <- proc.time() - start_time
cat(sprintf("[DONE] Step 1 complete in %.1f min\n", runtime[3] / 60))
cat("[DONE] Columns: enst | aa_start | aa_end | n_mer | ctex_up | ctex_dn\n")
cat("[DONE] Coordinates guaranteed correct -- derived from EnsDb ENST protein sequences\n")
cat("[DONE] Output in:", output_dir, "\n")
