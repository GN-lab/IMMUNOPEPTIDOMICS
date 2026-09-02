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
chunk_rows <- suppressWarnings(as.integer(Sys.getenv(
  "STEP4_CHUNK_ROWS", unset = "1000000"
)))
if (is.na(chunk_rows) || chunk_rows < 10000L) {
  stop("[ERROR] STEP4_CHUNK_ROWS must be an integer >= 10000")
}

cat("[CONFIG] OUTPUT_DIR: ", output_dir, "\n", sep = "")
cat("[CONFIG] Tier 1: NMP EL rank <0.5 AND MHCflurry affinity <500 nM\n")
cat("[CONFIG] Tier 2: NMP EL rank >=0.5 and <2.0 AND MHCflurry affinity <500 nM\n")
cat("[CONFIG] Tier 3: one tool binds and the other does not\n")
cat(sprintf("[CONFIG] Chunk rows: %d\n", chunk_rows))

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

# Stream a delimited file through bounded line chunks. The callback filters each
# chunk before it is retained, so memory scales with biologically relevant rows
# rather than with the approximately one-billion-row raw prediction matrices.
read_filtered_chunks <- function(path, sep, selected, column_names,
                                 label, filter_chunk) {
  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)

  header_line <- readLines(con, n = 1L, warn = FALSE)
  first_data <- readLines(con, n = 1L, warn = FALSE)
  if (!length(header_line) || !length(first_data)) {
    stop("[ERROR] Empty ", label, " file: ", path)
  }

  pending <- first_data
  retained_parts <- list()
  input_rows <- 0
  retained_rows <- 0
  chunk_number <- 0L

  repeat {
    need <- chunk_rows - length(pending)
    new_lines <- if (need > 0L) {
      readLines(con, n = need, warn = FALSE)
    } else {
      character()
    }
    n_new <- length(new_lines)
    lines <- c(pending, new_lines)
    pending <- character()
    if (!length(lines)) break

    chunk_number <- chunk_number + 1L
    input_rows <- input_rows + length(lines)
    chunk_text <- paste(lines, collapse = "\n")
    rm(lines, new_lines)

    dt <- fread(
      text = chunk_text,
      header = FALSE,
      sep = sep,
      select = selected,
      col.names = column_names,
      na.strings = c("", "NA"),
      showProgress = FALSE
    )
    rm(chunk_text)

    kept <- filter_chunk(dt)
    rm(dt)
    if (nrow(kept)) {
      retained_parts[[length(retained_parts) + 1L]] <- kept
      retained_rows <- retained_rows + nrow(kept)
    }

    if (chunk_number %% 10L == 0L || n_new < need) {
      cat(sprintf(
        "[INFO] %s chunk %d: %.0f rows scanned; %.0f candidate rows retained\n",
        label, chunk_number, input_rows, retained_rows
      ))
    }
    if (n_new < need) break
  }

  if (!length(retained_parts)) {
    stop("[ERROR] No relevant rows retained from ", label, ": ", path)
  }
  rbindlist(retained_parts, use.names = TRUE)
}

read_netmhcpan_relevant <- function(path) {
  header_line <- readLines(path, n = 1L, warn = FALSE)
  header <- strsplit(header_line, "\t", fixed = TRUE)[[1L]]
  required <- c(
    "allele", "peptide", "netmhcpan_EL_score", "netmhcpan_EL_rank",
    "netmhcpan_BA_score", "netmhcpan_BA_rank"
  )
  selected <- match(required, header)
  if (anyNA(selected)) {
    stop("[ERROR] Unsupported NetMHCpan header: ", paste(header, collapse = "\t"))
  }

  nmp <- read_filtered_chunks(
    path = path,
    sep = "\t",
    selected = selected,
    column_names = required,
    label = "NetMHCpan",
    filter_chunk = function(dt) {
      dt[, allele := normalize_allele(allele)]
      dt[, netmhcpan_EL_rank := as.numeric(netmhcpan_EL_rank)]
      dt[, netmhcpan_EL_score := as.numeric(netmhcpan_EL_score)]
      dt <- dt[!is.na(netmhcpan_EL_rank) & netmhcpan_EL_rank < NMP_WB_RANK]
      if (nrow(dt)) {
        setorder(dt, allele, peptide, netmhcpan_EL_rank,
                 -netmhcpan_EL_score)
        dt <- unique(dt, by = c("allele", "peptide"))
      }
      dt
    }
  )

  # A peptide-allele pair can straddle a line-chunk boundary. Apply the same
  # ordering and de-duplication once globally after raw rows are filtered.
  setorder(nmp, allele, peptide, netmhcpan_EL_rank, -netmhcpan_EL_score)
  unique(nmp, by = c("allele", "peptide"))
}

# Step 3a previously wrote an 8-field header above 11-field data rows.
# In that legacy file affinity is field 7 and presentation score is field 10.
read_mhcflurry_binders <- function(path) {
  probe <- readLines(path, n = 2L, warn = FALSE)
  if (length(probe) < 2L) stop("[ERROR] Empty MHCflurry file: ", path)
  header <- strsplit(probe[1L], ",", fixed = TRUE)[[1L]]
  header <- trimws(gsub('^"|"$', "", header))
  data_n <- lengths(strsplit(probe[2L], ",", fixed = TRUE))

  if (length(header) == 8L && data_n == 11L) {
    cat("[WARN] Legacy 8-header/11-data MHCflurry CSV detected; using fields 1,2,7,10\n")
    selected <- c(1L, 2L, 7L, 10L)
  } else {
    affinity <- intersect(
      c("mhcflurry_affinity", "mhcflurry_binding_affinity"), header
    )[1L]
    required <- c("peptide", "allele", "mhcflurry_presentation_score")
    if (length(setdiff(required, header)) || is.na(affinity)) {
      stop("[ERROR] Unsupported MHCflurry header: ", paste(header, collapse = ","))
    }
    selected <- match(
      c("peptide", "allele", affinity, "mhcflurry_presentation_score"),
      header
    )
  }

  output_names <- c(
    "peptide", "allele", "mhcflurry_affinity",
    "mhcflurry_presentation_score"
  )
  mhc <- read_filtered_chunks(
    path = path,
    sep = ",",
    selected = selected,
    column_names = output_names,
    label = "MHCflurry",
    filter_chunk = function(dt) {
      dt[, allele := normalize_allele(allele)]
      dt[, mhcflurry_affinity := as.numeric(mhcflurry_affinity)]
      dt[, mhcflurry_presentation_score :=
           as.numeric(mhcflurry_presentation_score)]
      dt <- dt[!is.na(mhcflurry_affinity) &
               mhcflurry_affinity < MHC_AFF_NM]
      if (nrow(dt)) {
        setorder(dt, allele, peptide, mhcflurry_affinity,
                 -mhcflurry_presentation_score)
        dt <- unique(dt, by = c("allele", "peptide"))
      }
      dt
    }
  )

  if (all(is.na(mhc$mhcflurry_affinity)) ||
      all(is.na(mhc$mhcflurry_presentation_score))) {
    stop("[ERROR] Parsed MHCflurry score columns contain only NA")
  }
  setorder(mhc, allele, peptide, mhcflurry_affinity,
           -mhcflurry_presentation_score)
  unique(mhc, by = c("allele", "peptide"))
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
  # Only NMP SB/WB rows can enter a retained concordance tier by themselves.
  # NMP NB rows are recovered implicitly when a MHC binder has no NMP SB/WB row.
  nmp_relevant <- read_netmhcpan_relevant(nmp_file)
  cat(sprintf("[INFO] NMP relevant (EL rank <2): %d rows\n", nrow(nmp_relevant)))
  gc(verbose = FALSE)

  cat("[INFO] MHCflurry: ", basename(mhc_file), "\n", sep = "")
  mhc_binder <- read_mhcflurry_binders(mhc_file)
  cat(sprintf("[INFO] MHC binders (<500 nM): %d rows\n", nrow(mhc_binder)))
  gc(verbose = FALSE)

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
