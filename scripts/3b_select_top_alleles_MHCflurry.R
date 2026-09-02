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

# Read the full Step 3a CSV in bounded chunks and retain only the current
# best allele for each peptide. The 70-allele CSV is expected to approach a
# billion rows, so loading it with one fread() call can exceed 64 GB. Memory
# here scales with the number of unique peptides, not peptide x allele rows.
read_mhcflurry_top_chunked <- function(path) {
  chunk_rows <- suppressWarnings(as.integer(Sys.getenv(
    "MHCFLURRY_CHUNK_ROWS", unset = "1000000"
  )))
  if (is.na(chunk_rows) || chunk_rows < 1000L) {
    stop("[ERROR] MHCFLURRY_CHUNK_ROWS must be an integer >= 1000")
  }

  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)

  header_line <- readLines(con, n = 1L, warn = FALSE)
  first_data  <- readLines(con, n = 1L, warn = FALSE)
  if (length(header_line) == 0L || length(first_data) == 0L) {
    stop("[ERROR] Empty MHCflurry file: ", path)
  }

  header <- strsplit(header_line, ",", fixed = TRUE)[[1L]]
  header <- trimws(gsub('^"|"$', "", header))
  data_n <- length(strsplit(first_data, ",", fixed = TRUE)[[1L]])

  if (length(header) == 8L && data_n == 11L) {
    cat("[WARN] Detected legacy 8-header/11-data CSV; parsing fields 1,2,7,10\n")
    selected <- c(1L, 2L, 7L, 10L)
  } else {
    affinity <- intersect(
      c("mhcflurry_affinity", "mhcflurry_binding_affinity"), header
    )[1L]
    required <- c("peptide", "allele", "mhcflurry_presentation_score")
    missing <- setdiff(required, header)
    if (length(missing) || is.na(affinity)) {
      stop("[ERROR] Unsupported MHCflurry columns in ", path,
           "; header: ", paste(header, collapse = ","))
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
  pending <- first_data
  best <- NULL
  input_rows <- 0
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
    if (length(lines) == 0L) break

    chunk_number <- chunk_number + 1L
    input_rows <- input_rows + length(lines)
    chunk_text <- paste(lines, collapse = "\n")
    rm(lines, new_lines)

    dt <- fread(
      text = chunk_text, header = FALSE, select = selected,
      col.names = output_names, na.strings = c("", "NA"),
      showProgress = FALSE
    )
    rm(chunk_text)

    dt[, mhcflurry_presentation_score :=
         as.numeric(mhcflurry_presentation_score)]
    dt[, mhcflurry_affinity := as.numeric(mhcflurry_affinity)]
    dt <- dt[
      !is.na(peptide) & peptide != "" & !is.na(allele) & allele != "" &
      !is.na(mhcflurry_presentation_score)
    ]

    # Reduce duplicate peptides inside this chunk first.
    setorder(dt, peptide, -mhcflurry_presentation_score)
    chunk_top <- dt[, .SD[1L], by = peptide]
    rm(dt)

    if (is.null(best)) {
      best <- chunk_top
      setkey(best, peptide)
    } else {
      setkey(chunk_top, peptide)

      # Update existing peptides only when the new score is strictly higher;
      # ties retain the first allele, matching the previous stable ordering.
      updates <- best[chunk_top, on = .(peptide), nomatch = 0L,
        .(
          peptide,
          old_score = mhcflurry_presentation_score,
          new_allele = i.allele,
          new_affinity = i.mhcflurry_affinity,
          new_score = i.mhcflurry_presentation_score
        )
      ]
      updates <- updates[new_score > old_score]
      if (nrow(updates) > 0L) {
        best[updates, on = .(peptide), `:=`(
          allele = i.new_allele,
          mhcflurry_affinity = i.new_affinity,
          mhcflurry_presentation_score = i.new_score
        )]
      }

      new_peptides <- chunk_top[!best, on = .(peptide)]
      if (nrow(new_peptides) > 0L) {
        best <- rbindlist(list(best, new_peptides), use.names = TRUE)
        setkey(best, peptide)
      }
      rm(chunk_top, updates, new_peptides)
    }

    cat(sprintf(
      "[INFO] Chunk %d: %.0f input rows scanned; %d peptide maxima retained\n",
      chunk_number, input_rows, nrow(best)
    ))
    gc(verbose = FALSE)

    if (n_new < need) break
  }

  if (is.null(best) || nrow(best) == 0L) {
    stop("[ERROR] No valid MHCflurry rows parsed from: ", path)
  }
  setorder(best, peptide)
  best
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

  top <- read_mhcflurry_top_chunked(f)
  cat(sprintf("[INFO] %d unique peptides after best-allele selection\n", nrow(top)))

  out_file <- file.path(directory_out,
    paste0("canonical_", len, "mer_mhcflurry_top_", current_date, ".tsv"))
  fwrite(top, out_file, sep = "\t", col.names = TRUE, quote = FALSE)
  cat(sprintf("[INFO] Written: %s\n", out_file))

  rm(top); gc()
}

cat("\n[DONE] Step 3b complete.\n")
