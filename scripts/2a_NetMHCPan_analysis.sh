#!/bin/bash
#SBATCH --job-name=immunopep_netmhcpan_9mer
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=8042
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/2a_netmhcpan_%j.log
#SBATCH --error=logs/2a_netmhcpan_%j.err

set -euo pipefail

# ============================================================
# JUGNU_IMMUNOPEPTIDOME — Step 2a: NetMHCpan 4.2 Predictions
# July 2026 | Gaurav Raichand | The Institute of Cancer Research
#
# Input:  canonical_09mers_YYYYMMDD.tsv  (from Step 1, col: n_mer)
#         HLA_ALLELE_FILE                 (standard alleles, from config.sh)
# Output: canonical_09mer_netmhcpan_YYYY_MMDD.tsv
#
# NOTE: Running 9mers only for now.
#       To run all lengths, change the loop at the bottom to: for LEN in 8 9 10 11
#
# Design:
#   - Alleles batched 20 at a time (NetMHCpan -a char limit = 1024)
#   - Peptides chunked at CHUNK_SIZE lines to control memory (~12M per length)
#   - Existing raw .txt skipped on re-run (safe restart)
#   - Zero hardcoded paths -- everything from config.sh
# ============================================================

source "${WORKDIR}/config.sh"

: "${OUTPUT_DIR:?OUTPUT_DIR not set}"
: "${HLA_ALLELE_FILE:?HLA_ALLELE_FILE not set -- run scripts/setup_alleles.sh first}"
: "${NETMHCPAN_BIN:?NETMHCPAN_BIN not set}"
: "${CHUNK_SIZE:?CHUNK_SIZE not set}"

# An explicit date can still be supplied, but each peptide length also checks
# for a verified completed raw file whose TSV is missing and resumes that parse
# automatically. This prevents a next-day retry from starting predictions again.
RUN_DATE=${NETMHCPAN_RUN_DATE:-$(date +%Y_%m%d)}
echo "[CONFIG] OUTPUT_DIR:      ${OUTPUT_DIR}"
echo "[CONFIG] HLA_ALLELE_FILE: ${HLA_ALLELE_FILE}"
echo "[CONFIG] CHUNK_SIZE:      ${CHUNK_SIZE}"
echo "[CONFIG] RUN_DATE:        ${RUN_DATE}"

###########################################################################
# 0. Validate NetMHCpan
###########################################################################

if [[ ! -x "${NETMHCPAN_BIN}" ]]; then
  echo "[ERROR] NetMHCpan not executable: ${NETMHCPAN_BIN}"; exit 1
fi
echo "GILGFVFTL" > "${TMPDIR}/test_ip.pep"
if ! "${NETMHCPAN_BIN}" -a HLA-A02:01 -p "${TMPDIR}/test_ip.pep" -l 9 \
    > "${TMPDIR}/test_ip.out" 2>&1; then
  echo "[ERROR] NetMHCpan sanity check failed:"; cat "${TMPDIR}/test_ip.out"; exit 1
fi
echo "[INFO] NetMHCpan validated: ${NETMHCPAN_BIN}"

###########################################################################
# 1. Load standard HLA allele list
###########################################################################

if [[ ! -f "${HLA_ALLELE_FILE}" ]] || [[ ! -s "${HLA_ALLELE_FILE}" ]]; then
  echo "[ERROR] HLA allele file not found: ${HLA_ALLELE_FILE}"
  echo "[ERROR] Run once: bash scripts/setup_alleles.sh"
  exit 1
fi

ALLELE_LIST=$(paste -sd "," "${HLA_ALLELE_FILE}")
N_ALLELES=$(wc -l < "${HLA_ALLELE_FILE}")
echo "[INFO] ${N_ALLELES} standard HLA-A/B/C alleles loaded"

###########################################################################
# 2. Run NetMHCpan per length (sequential, chunked)
###########################################################################

run_netmhcpan_canonical() {
  local LENGTH=$1
  local LENPAD
  LENPAD=$(printf "%02d" "${LENGTH}")

  # Find newest Step 1 output for this length
  local INPUT_PEP
  INPUT_PEP=$(find "${OUTPUT_DIR}" -maxdepth 1 -type f \
    -name "canonical_${LENPAD}mer_netmhcpan_input_*.pep" | sort | tail -n 1)

  if [[ -z "${INPUT_PEP}" ]] || [[ ! -s "${INPUT_PEP}" ]]; then
    echo "[WARN] No canonical_${LENPAD}mer_netmhcpan_input_*.pep found -- skipping"
    echo "[WARN] Has Step 1 finished successfully?"
    return
  fi

  local EFFECTIVE_RUN_DATE="${RUN_DATE}"

  # If no date was explicitly requested, resume the newest verified raw file
  # that completed prediction but does not yet have its parsed TSV.
  if [[ -z "${NETMHCPAN_RUN_DATE:-}" ]]; then
    local RESUME_MARKER=""
    local CANDIDATE_MARKER CANDIDATE_RAW CANDIDATE_TSV
    while IFS= read -r CANDIDATE_MARKER; do
      [[ -z "${CANDIDATE_MARKER}" ]] && continue
      CANDIDATE_RAW="${CANDIDATE_MARKER%.complete}"
      CANDIDATE_TSV="${CANDIDATE_RAW%.txt}.tsv"
      if [[ -s "${CANDIDATE_RAW}" ]] && [[ ! -s "${CANDIDATE_TSV}" ]]; then
        RESUME_MARKER="${CANDIDATE_MARKER}"
      fi
    done < <(
      find "${OUTPUT_DIR}" -maxdepth 1 -type f \
        -name "canonical_${LENPAD}mer_netmhcpan_*.txt.complete" | sort
    )

    if [[ -n "${RESUME_MARKER}" ]]; then
      local RESUME_BASENAME
      RESUME_BASENAME=$(basename "${RESUME_MARKER%.complete}")
      EFFECTIVE_RUN_DATE=${RESUME_BASENAME#canonical_${LENPAD}mer_netmhcpan_}
      EFFECTIVE_RUN_DATE=${EFFECTIVE_RUN_DATE%.txt}
      echo "[INFO] Auto-resuming verified raw prediction date: ${EFFECTIVE_RUN_DATE}"
    fi
  fi

  local OUT_RAW="${OUTPUT_DIR}/canonical_${LENPAD}mer_netmhcpan_${EFFECTIVE_RUN_DATE}.txt"
  local OUT_TSV="${OUTPUT_DIR}/canonical_${LENPAD}mer_netmhcpan_${EFFECTIVE_RUN_DATE}.tsv"
  local RAW_COMPLETE="${OUT_RAW}.complete"

  local N_PEPS
  N_PEPS=$(wc -l < "${INPUT_PEP}")
  echo ""
  echo "[INFO] === ${LENPAD}mer | ${N_PEPS} unique peptides | $(basename ${INPUT_PEP}) ==="

  if [[ -s "${OUT_RAW}" ]] && [[ -f "${RAW_COMPLETE}" ]]; then
    echo "[INFO] Verified raw .txt already exists -- skipping NetMHCpan, re-parsing only"
  elif [[ -e "${OUT_RAW}" ]] || [[ -e "${RAW_COMPLETE}" ]]; then
    echo "[ERROR] Refusing to use or overwrite an unverified raw prediction:" >&2
    echo "[ERROR]   raw:    ${OUT_RAW}" >&2
    echo "[ERROR]   marker: ${RAW_COMPLETE}" >&2
    echo "[ERROR] A complete raw file must have its matching .txt.complete marker." >&2
    exit 1
  else

    if [[ ${N_PEPS} -eq 0 ]]; then
      echo "[WARN] No peptides for ${LENPAD}mer -- skipping"; return
    fi

    # Write predictions to a job-specific partial file. Only atomically move it
    # to OUT_RAW and create the completion marker after every chunk succeeds.
    local OUT_RAW_PARTIAL="${OUT_RAW}.partial.${SLURM_JOB_ID:-$$}"
    : > "${OUT_RAW_PARTIAL}"

    local CHUNK_PREFIX="${TMPDIR}/ip_chunk_${LENPAD}_"
    rm -f "${CHUNK_PREFIX}"*
    split -l "${CHUNK_SIZE}" "${INPUT_PEP}" "${CHUNK_PREFIX}"
    local N_CHUNKS
    N_CHUNKS=$(ls "${CHUNK_PREFIX}"* | wc -l)
    echo "[INFO] ${N_CHUNKS} chunk(s) of <= ${CHUNK_SIZE} peptides"

    local BATCH_SIZE=20
    local ALLELES_ARRAY
    IFS=',' read -ra ALLELES_ARRAY <<< "${ALLELE_LIST}"
    local TOTAL=${#ALLELES_ARRAY[@]}

    local CHUNK_NUM=0
    for CHUNK_FILE in "${CHUNK_PREFIX}"*; do
      CHUNK_NUM=$((CHUNK_NUM + 1))
      local NLINES
      NLINES=$(wc -l < "${CHUNK_FILE}")
      echo "[INFO]   Chunk ${CHUNK_NUM}/${N_CHUNKS}: ${NLINES} peptides"

      for (( i=0; i<TOTAL; i+=BATCH_SIZE )); do
        local BATCH_END=$(( i + BATCH_SIZE < TOTAL ? i + BATCH_SIZE : TOTAL ))
        local BATCH_ALLELES
        BATCH_ALLELES=$(IFS=','; echo "${ALLELES_ARRAY[*]:$i:$BATCH_SIZE}")
        "${NETMHCPAN_BIN}" \
          -a "${BATCH_ALLELES}" \
          -p "${CHUNK_FILE}" \
          -l "${LENGTH}" \
          -BA \
          >> "${OUT_RAW_PARTIAL}"
      done
      rm -f "${CHUNK_FILE}"
    done

    mv "${OUT_RAW_PARTIAL}" "${OUT_RAW}"
    {
      echo "status=complete"
      echo "run_date=${EFFECTIVE_RUN_DATE}"
      echo "alleles=${N_ALLELES}"
      echo "peptides=${N_PEPS}"
      echo "raw_bytes=$(stat -c%s "${OUT_RAW}")"
    } > "${RAW_COMPLETE}"
    echo "[INFO] Raw output: ${OUT_RAW}"
    echo "[INFO] Completion marker: ${RAW_COMPLETE}"
  fi

  # Parse raw -> complete TSV in bounded chunks.
  #
  # Do not accumulate the full 70-allele output in memory: the raw output is
  # >100 GB. Keep every valid prediction (SB, WB and NB), but flush parsed rows
  # to disk every PARSE_CHUNK_SIZE records. Step 2b remains responsible for
  # selecting the strong binders. Write atomically through a temporary file so
  # an interrupted parse cannot be mistaken for a completed TSV on restart.
  python3 - <<PYEOF
import os

infile, outfile = "${OUT_RAW}", "${OUT_TSV}"
tmpfile = outfile + ".tmp"
PARSE_CHUNK_SIZE = 100_000
raw_lines = parsed_rows = chunks_written = 0
row_chunk = []

if os.path.exists(tmpfile):
    os.remove(tmpfile)

with open(infile, "r", errors="replace", buffering=1024 * 1024) as fh, \
     open(tmpfile, "w", buffering=1024 * 1024) as out:
    out.write("allele\tpeptide\tcore\t"
              "netmhcpan_EL_score\tnetmhcpan_EL_rank\t"
              "netmhcpan_BA_score\tnetmhcpan_BA_rank\tbinder\n")

    for raw_lines, line in enumerate(fh, start=1):
        if raw_lines % 10_000_000 == 0:
            print(
                f"[INFO] Parse progress: {raw_lines:,} raw lines; "
                f"{parsed_rows:,} prediction rows written",
                flush=True,
            )

        line = line.rstrip()
        if not line or line[0] in ("#", "-") or \
           line.startswith(" Pos") or line.startswith("Protein") or \
           line.startswith("Error"):
            continue

        p = line.split()
        if len(p) < 13:
            continue

        allele, peptide, core = p[1], p[2], p[3]
        if not allele.startswith("HLA-") or not peptide.isalpha():
            continue

        try:
            el_score, el_rank = p[11], p[12]
            ba_score = p[13] if len(p) > 13 else "NA"
            ba_rank = p[14] if len(p) > 14 else "NA"
            float(el_score)
            float(el_rank)
            if ba_score != "NA":
                float(ba_score)
            if ba_rank != "NA":
                float(ba_rank)
        except (IndexError, ValueError):
            continue

        # NetMHCpan appends "<= SB" or "<= WB" only to binder rows.
        binder = (p[-1]
                  if len(p) >= 2 and p[-2] == "<=" and p[-1] in ("SB", "WB")
                  else "NB")

        row_chunk.append("\t".join([
            allele, peptide, core,
            el_score, el_rank, ba_score, ba_rank, binder,
        ]) + "\n")
        parsed_rows += 1

        if len(row_chunk) >= PARSE_CHUNK_SIZE:
            out.writelines(row_chunk)
            row_chunk.clear()
            chunks_written += 1

    if row_chunk:
        out.writelines(row_chunk)
        row_chunk.clear()
        chunks_written += 1

if parsed_rows == 0:
    if os.path.exists(tmpfile):
        os.remove(tmpfile)
    raise RuntimeError("No NetMHCpan prediction rows were parsed")

os.replace(tmpfile, outfile)
print(
    f"[INFO] Parsed and wrote {parsed_rows:,} prediction rows "
    f"in {chunks_written:,} chunks -> {outfile}",
    flush=True,
)
PYEOF

  echo "[INFO] ${LENPAD}mer TSV: ${OUT_TSV}"
}

echo ""
echo "=== Running NetMHCpan on canonical immunopeptidome (9mer only) ==="

# NOTE: Change to "for LEN in 8 9 10 11" to run all lengths
for LEN in 9; do
  run_netmhcpan_canonical ${LEN}
done

echo ""
echo "=== Summary ==="
for LEN in 9; do
  LENPAD=$(printf "%02d" "${LEN}")
  TSV="${OUTPUT_DIR}/canonical_${LENPAD}mer_netmhcpan_${RUN_DATE}.tsv"
  if [[ -f "${TSV}" ]]; then
    echo "  ${LENPAD}mer: $(wc -l < "${TSV}") rows -> $(basename ${TSV})"
  else
    echo "  ${LENPAD}mer: MISSING"
  fi
done
