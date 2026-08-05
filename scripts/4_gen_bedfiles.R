#!/usr/bin/env Rscript

# ============================================================
# July 2026 | Gaurav Raichand | The Institute of Cancer Research
#
# JUGNU immunopeptidome - Step 4 
# Cross-tool concordance tiering followed by Tier-1 BED generation.
# ============================================================

suppressPackageStartupMessages(library(data.table))

start_time <- proc.time()
output_dir <- Sys.getenv("OUTPUT_DIR")
if (output_dir == "") stop("[ERROR] OUTPUT_DIR not set -- source config.sh first.")

NMP_SB_RANK <- 0.5
NMP_WB_RANK <- 2.0
MHC_AFF_NM <- 500
run_date <- format(Sys.Date(), "%Y%m%d")
lengths <- c("08", "09", "10", "11")

cat("[CONFIG] OUTPUT_DIR: ", output_dir, "\n", sep = "")
cat("[CONFIG] Tier 1: NMP EL rank <0.5 AND MHCflurry affinity <500 nM\n")
cat("[CONFIG] Tier 2: NMP EL rank >=0.5 and <2.0 AND MHCflurry affinity <500 nM\n")
cat("[CONFIG] Tier 3: one tool binds and the other does not\n")

newest <- function(pattern) {
  files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  if (!length(files)) return(NULL)
  files[which.max(file.info(files)$mtime)]
}

normalize_allele <- function(x) {
  toupper(gsub("\\*", "", trimws(x)))
}

bed_allele <- function(x) {
  gsub("[*:]", "", normalize_allele(x))
}

# Step 3a previously wrote an 8-field header above 11-field data rows.
# In that legacy file affinity is field 7 and presentation score is field 10.
read_mhcflurry <- function(path) {
  probe <- readLines(path, n = 2L, warn = FALSE)
  if (length(probe) < 2L) stop("[ERROR] Empty MHCflurry file: ", path)
  header_n <- lengths(strsplit(probe[1L], ",", fixed = TRUE))
  data_n <- lengths(strsplit(probe[2L], ",", fixed = TRUE))

  if (header_n == 8L && data_n == 11L) {
    cat("[WARN] Legacy 8-header/11-data MHCflurry CSV detected; using fields 1,2,7,10\n")
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
  if (length(setdiff(required, header)) || is.na(affinity)) {
    stop("[ERROR] Unsupported MHCflurry header: ", paste(header, collapse = ","))
  }

  dt <- fread(
    path,
    select = c("peptide", "allele", affinity,
               "mhcflurry_presentation_score"),
    na.strings = c("", "NA")
  )
  if (affinity != "mhcflurry_affinity") {
    setnames(dt, affinity, "mhcflurry_affinity")
  }
  dt
}

is_low_complexity <- function(pep) {
  chars <- strsplit(as.character(pep), "", fixed = TRUE)[[1L]]
  max(tabulate(match(chars, unique(chars)))) / length(chars) > 0.6
}

retained_parts <- vector("list", length(lengths))

for (idx in seq_along(lengths)) {
  len <- lengths[idx]
  cat(sprintf("\n=== Cross-analysis: %smer ===\n", len))

  nmp_file <- newest(sprintf("^canonical_%smer_netmhcpan_[0-9_]+\\.tsv$", len))
  mhc_file <- newest(sprintf("^canonical_%smers_flank_mhcflurry_[0-9_]+\\.csv$", len))
  if (is.null(nmp_file) || is.null(mhc_file)) {
    cat(sprintf("[WARN] Missing one or both raw prediction files for %smer - skipping\n", len))
    next
  }

  cat("[INFO] NetMHCpan: ", basename(nmp_file), "\n", sep = "")
  nmp <- fread(
    nmp_file,
    select = c("allele", "peptide", "netmhcpan_EL_score",
               "netmhcpan_EL_rank", "netmhcpan_BA_score",
               "netmhcpan_BA_rank"),
    na.strings = c("", "NA")
  )
  nmp[, allele := normalize_allele(allele)]
  nmp[, netmhcpan_EL_rank := as.numeric(netmhcpan_EL_rank)]
  nmp[, netmhcpan_EL_score := as.numeric(netmhcpan_EL_score)]
  setorder(nmp, allele, peptide, netmhcpan_EL_rank, -netmhcpan_EL_score)
  nmp <- unique(nmp, by = c("allele", "peptide"))

  # Only NMP SB/WB rows can enter a retained concordance tier by themselves.
  # NMP NB rows are recovered implicitly when a MHC binder has no NMP SB/WB row.
  nmp_relevant <- nmp[!is.na(netmhcpan_EL_rank) & netmhcpan_EL_rank < NMP_WB_RANK]
  cat(sprintf("[INFO] NMP relevant (EL rank <2): %d rows\n", nrow(nmp_relevant)))
  rm(nmp); gc(verbose = FALSE)

  cat("[INFO] MHCflurry: ", basename(mhc_file), "\n", sep = "")
  mhc <- read_mhcflurry(mhc_file)
  mhc[, allele := normalize_allele(allele)]
  mhc[, mhcflurry_affinity := as.numeric(mhcflurry_affinity)]
  mhc[, mhcflurry_presentation_score := as.numeric(mhcflurry_presentation_score)]
  if (all(is.na(mhc$mhcflurry_affinity)) ||
      all(is.na(mhc$mhcflurry_presentation_score))) {
    stop("[ERROR] Parsed MHCflurry score columns contain only NA")
  }
  setorder(mhc, allele, peptide, mhcflurry_affinity,
           -mhcflurry_presentation_score)
  mhc <- unique(mhc, by = c("allele", "peptide"))
  mhc_binder <- mhc[!is.na(mhcflurry_affinity) &
                    mhcflurry_affinity < MHC_AFF_NM]
  cat(sprintf("[INFO] MHC binders (<500 nM): %d rows\n", nrow(mhc_binder)))
  rm(mhc); gc(verbose = FALSE)

  # Full join of all potentially retained calls. Missing tool values mean NB.
  cross <- merge(
    nmp_relevant, mhc_binder,
    by = c("allele", "peptide"), all = TRUE
  )
  rm(nmp_relevant, mhc_binder); gc(verbose = FALSE)

  cross[, nmp_call := fcase(
    !is.na(netmhcpan_EL_rank) & netmhcpan_EL_rank < NMP_SB_RANK, "SB",
    !is.na(netmhcpan_EL_rank) & netmhcpan_EL_rank < NMP_WB_RANK, "WB",
    default = "NB"
  )]
  cross[, mhc_call := fifelse(
    !is.na(mhcflurry_affinity) & mhcflurry_affinity < MHC_AFF_NM,
    "binder", "NB"
  )]
  cross[, concordance_tier := fcase(
    nmp_call == "SB" & mhc_call == "binder", "Tier1_HighConfidence",
    nmp_call == "WB" & mhc_call == "binder", "Tier2_MediumConfidence",
    nmp_call == "SB" & mhc_call == "NB", "Tier3_Discordant_NMPstrong",
    nmp_call %in% c("WB", "NB") & mhc_call == "binder",
      "Tier3_Discordant_MHCstrong",
    default = "Excluded"
  )]
  cross[, combined_score := rowMeans(
    cbind(netmhcpan_EL_score, mhcflurry_presentation_score),
    na.rm = TRUE
  )]
  cross[, peptide_length := as.integer(len)]
  retained <- cross[concordance_tier != "Excluded"]
  cat("[INFO] Tier counts:\n")
  print(retained[, .N, by = concordance_tier][order(concordance_tier)])
  retained_parts[[idx]] <- retained
  rm(cross, retained); gc(verbose = FALSE)
}

retained <- rbindlist(retained_parts, use.names = TRUE, fill = TRUE)
if (!nrow(retained)) stop("[ERROR] No retained concordance candidates")

tier1 <- retained[concordance_tier == "Tier1_HighConfidence"]
tier2 <- retained[concordance_tier == "Tier2_MediumConfidence"]
tier3 <- retained[grepl("^Tier3", concordance_tier)]

setorder(tier1, peptide_length, netmhcpan_EL_rank, mhcflurry_affinity)
setorder(tier2, peptide_length, netmhcpan_EL_rank, mhcflurry_affinity)
setorder(tier3, peptide_length, concordance_tier, netmhcpan_EL_rank,
         mhcflurry_affinity)

tier1_file <- file.path(output_dir, sprintf("canonical_concordance_tier1_%s.tsv", run_date))
tier2_file <- file.path(output_dir, sprintf("canonical_concordance_tier2_%s.tsv", run_date))
tier3_file <- file.path(output_dir, sprintf("canonical_concordance_tier3_%s.tsv", run_date))
all_file <- file.path(output_dir, sprintf("canonical_concordance_all_%s.tsv", run_date))

fwrite(tier1, tier1_file, sep = "\t", na = "NA")
fwrite(tier2, tier2_file, sep = "\t", na = "NA")
fwrite(tier3, tier3_file, sep = "\t", na = "NA")
fwrite(retained, all_file, sep = "\t", na = "NA")

cat(sprintf("\n[INFO] Tier 1 rows: %d\n", nrow(tier1)))
cat(sprintf("[INFO] Tier 2 rows: %d\n", nrow(tier2)))
cat(sprintf("[INFO] Tier 3 rows: %d\n", nrow(tier3)))

# BED is Tier 1 only, matching SSNIP Step 15d.
tier1 <- tier1[!sapply(peptide, is_low_complexity)]
cat(sprintf("[INFO] Tier 1 after low-complexity filter: %d rows\n", nrow(tier1)))

mapped_parts <- vector("list", length(lengths))
for (idx in seq_along(lengths)) {
  len <- lengths[idx]
  candidates <- tier1[peptide_length == as.integer(len)]
  if (!nrow(candidates)) next

  coord_file <- newest(sprintf("^canonical_%smers_[0-9]{8}\\.tsv$", len))
  if (is.null(coord_file)) stop("[ERROR] Missing coordinate file for ", len, "mer")
  coords <- fread(
    coord_file,
    select = c("enst", "aa_start", "aa_end", "n_mer")
  )
  setnames(coords, "n_mer", "peptide")
  coords <- unique(coords, by = c("enst", "aa_start", "aa_end", "peptide"))
  mapped_parts[[idx]] <- merge(
    candidates, coords,
    by = "peptide", allow.cartesian = TRUE
  )
  rm(coords); gc(verbose = FALSE)
}

mapped <- rbindlist(mapped_parts, use.names = TRUE, fill = TRUE)
if (!nrow(mapped)) stop("[ERROR] No Tier-1 peptides mapped to coordinates")

expected_pairs <- unique(tier1[, .(peptide, allele)])
observed_pairs <- unique(mapped[, .(peptide, allele)])
if (nrow(fsetdiff(expected_pairs, observed_pairs))) {
  stop("[ERROR] One or more Tier-1 peptide-HLA pairs lack coordinates")
}

mapped_tsv <- file.path(
  output_dir,
  sprintf("canonical_immunopeptidome_SB_%s.tsv", run_date)
)
fwrite(mapped, mapped_tsv, sep = "\t", na = "NA")

# Convert to BED coordinates, then collapse every Tier-1 allele into column 5.
mapped[, `:=`(
  bed_start = as.integer(aa_start) - 1L,
  bed_end = as.integer(aa_end),
  bed_hla = bed_allele(allele)
)]
if (anyNA(mapped[, .(enst, bed_start, bed_end, peptide, bed_hla)])) {
  stop("[ERROR] Missing BED values")
}
if (any(mapped$bed_start < 0L) ||
    any((mapped$bed_end - mapped$bed_start) != nchar(mapped$peptide))) {
  stop("[ERROR] Invalid BED coordinates")
}

mapped_pairs <- unique(mapped[, .(
  enst, bed_start, bed_end, peptide, bed_hla
)])
bed <- mapped_pairs[, .(
  HLA_ALLELES = paste(sort(unique(bed_hla)), collapse = ",")
), by = .(
  ENST_ID = enst,
  AA_START = bed_start,
  AA_END = bed_end,
  PEPTIDE = peptide
)]
setorder(bed, ENST_ID, AA_START, AA_END, PEPTIDE)

# Prove no allele association was lost during collapse.
collapsed_associations <- sum(lengths(strsplit(bed$HLA_ALLELES, ",", fixed = TRUE)))
if (collapsed_associations != nrow(mapped_pairs)) {
  stop("[ERROR] HLA allele loss detected during BED collapse")
}

bed_file <- file.path(
  output_dir,
  sprintf("canonical_immunopeptidome_SB_%s.bed", run_date)
)
fwrite(bed, bed_file, sep = "\t", col.names = FALSE, quote = FALSE)

cat("\n=== Step 4 summary ===\n")
cat(sprintf("Tier 1: %d rows\n", nrow(tier1)))
cat(sprintf("Tier 2: %d rows\n", nrow(tier2)))
cat(sprintf("Tier 3: %d rows\n", nrow(tier3)))
cat(sprintf("BED rows: %d\n", nrow(bed)))
cat(sprintf("BED allele associations: %d\n", collapsed_associations))
cat(sprintf("BED unique ENST IDs: %d\n", uniqueN(bed$ENST_ID)))
cat(sprintf("BED: %s\n", bed_file))
cat(sprintf("[DONE] Step 4 completed in %.1f minutes\n",
            (proc.time() - start_time)[3L] / 60))
