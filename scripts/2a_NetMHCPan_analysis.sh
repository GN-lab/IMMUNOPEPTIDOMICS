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

RUN_DATE=$(date +%Y_%m%d)
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

  local OUT_RAW="${OUTPUT_DIR}/canonical_${LENPAD}mer_netmhcpan_${RUN_DATE}.txt"
  local OUT_TSV="${OUTPUT_DIR}/canonical_${LENPAD}mer_netmhcpan_${RUN_DATE}.tsv"

  local N_PEPS
  N_PEPS=$(wc -l < "${INPUT_PEP}")
  echo ""
  echo "[INFO] === ${LENPAD}mer | ${N_PEPS} unique peptides | $(basename ${INPUT_PEP}) ==="

  if [[ -f "${OUT_RAW}" ]] && [[ -s "${OUT_RAW}" ]]; then
    echo "[INFO] Raw .txt already exists -- skipping NetMHCpan, re-parsing only"
  else

    if [[ ${N_PEPS} -eq 0 ]]; then
      echo "[WARN] No peptides for ${LENPAD}mer -- skipping"; return
    fi

    > "${OUT_RAW}"

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
          >> "${OUT_RAW}"
      done
      rm -f "${CHUNK_FILE}"
    done
    echo "[INFO] Raw output: ${OUT_RAW}"
  fi

  # Parse raw -> TSV
  python3 - <<PYEOF
infile, outfile = "${OUT_RAW}", "${OUT_TSV}"
rows = []
with open(infile) as fh:
    for line in fh:
        line = line.rstrip()
        if not line or line[0] in ("#", "-") or \
           line.startswith(" Pos") or line.startswith("Protein") or \
           line.startswith("Error"):
            continue
        p = line.split()
        if len(p) < 11:
            continue
        try:
            allele, peptide, core = p[1], p[2], p[3]
            el_score, el_rank     = p[11], p[12]
            ba_score = p[13] if len(p) > 13 else "NA"
            ba_rank  = p[14] if len(p) > 14 else "NA"
            binder = ("SB" if len(p) >= 16 and p[-2] == "<=" and p[-1] == "SB"
                      else "WB" if len(p) >= 16 and p[-2] == "<=" and p[-1] == "WB"
                      else "NB")
            if not allele.startswith("HLA") or not peptide.isalpha():
                continue
            rows.append([allele, peptide, core,
                         el_score, el_rank, ba_score, ba_rank, binder])
        except (IndexError, ValueError):
            continue

with open(outfile, "w") as out:
    out.write("allele\tpeptide\tcore\t"
              "netmhcpan_EL_score\tnetmhcpan_EL_rank\t"
              "netmhcpan_BA_score\tnetmhcpan_BA_rank\tbinder\n")
    for r in rows:
        out.write("\t".join(r) + "\n")
print(f"[INFO] Parsed {len(rows)} rows -> {outfile}")
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