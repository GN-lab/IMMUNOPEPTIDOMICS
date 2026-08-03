#!/bin/bash
#SBATCH --job-name=immunopep_mhcflurry_09mer
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=13600
#SBATCH --gpus-per-node=1
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/3a_mhcflurry_%j.log
#SBATCH --error=logs/3a_mhcflurry_%j.err

set -euo pipefail

# ============================================================
# JUGNU_IMMUNOPEPTIDOME — Step 3a: MHCFlurry 2.0 Predictions
# July 2026 | Gaurav Raichand | The Institute of Cancer Research
#
# Strategy: generate MHCflurry input ONE ALLELE AT A TIME to avoid
# holding the full peptide x allele cross-join in memory (10M x 105
# alleles = ~1B rows = hundreds of GB if materialised all at once).
#
# Per length per allele:
#   1. Build a small CSV: peptide | allele | n_flank | c_flank
#      (unique peptides only, ~10M rows per allele)
#   2. Run mhcflurry-predict on that CSV
#   3. Append scored rows to the combined output file
#   4. Delete the temporary input CSV
#
# Output: canonical_0Xmers_flank_mhcflurry_YYYYMMDD.csv  (x4)
#   One row per (peptide, allele) with mhcflurry scores.
#
# NOTE: Running 09mer only for now.
#       To run all lengths, change the loop at the bottom to:
#         for LEN in 08 09 10 11; do
# ============================================================

source "${WORKDIR}/config.sh"
: "${OUTPUT_DIR:?OUTPUT_DIR not set}"
: "${HLA_ALLELE_FILE:?HLA_ALLELE_FILE not set}"

# Override the neojunction_viz env activated by config.sh with mhcflurry-env
# (config.sh sources conda and loads modules -- we just need to switch env)
conda deactivate || true
conda activate /data/rds/DMP/UCEC/EVOLIMMU/graichand/.conda_envs/mhcflurry-env

RUN_DATE=$(date +%Y%m%d)
echo "[CONFIG] OUTPUT_DIR:      ${OUTPUT_DIR}"
echo "[CONFIG] HLA_ALLELE_FILE: ${HLA_ALLELE_FILE}"
echo "[CONFIG] RUN_DATE:        ${RUN_DATE}"

# Load allele list
mapfile -t ALLELES < "${HLA_ALLELE_FILE}"
N_ALLELES=${#ALLELES[@]}
echo "[INFO] ${N_ALLELES} HLA alleles loaded"

###########################################################################
# Process each length
# NOTE: Change to "for LEN in 08 09 10 11" to run all lengths
###########################################################################

for LEN in 09; do

  # Find coordinate TSV from Step 1 (has n_mer + flanks)
  COORD_TSV=$(find "${OUTPUT_DIR}" -maxdepth 1 -type f \
    -name "canonical_${LEN}mers_*.tsv" | sort | tail -n 1)

  if [[ -z "${COORD_TSV}" ]] || [[ ! -s "${COORD_TSV}" ]]; then
    echo "[WARN] No canonical_${LEN}mers_*.tsv found -- skipping ${LEN}mer"
    continue
  fi

  OUT_CSV="${OUTPUT_DIR}/canonical_${LEN}mers_flank_mhcflurry_${RUN_DATE}.csv"

  if [[ -f "${OUT_CSV}" ]] && [[ -s "${OUT_CSV}" ]]; then
    echo "[INFO] ${LEN}mer: output already exists -- skipping"
    continue
  fi

  echo ""
  echo "[INFO] === ${LEN}mer: $(basename ${COORD_TSV}) ==="

  # Extract unique peptides + flanks once (col 4=n_mer, 5=ctex_up, 6=ctex_dn)
  PEPS_TSV="${TMPDIR}/ip_${LEN}mer_peps.tsv"
  echo "[INFO] Extracting unique peptides and flanks..."
  awk -F'\t' 'NR>1 && $4!="" {print $4"\t"$5"\t"$6}' "${COORD_TSV}" \
    | sort -u \
    > "${PEPS_TSV}"
  N_PEPS=$(wc -l < "${PEPS_TSV}")
  echo "[INFO] ${N_PEPS} unique peptides"

  # Start empty. The first successful prediction supplies MHCflurry's native
  # header so the header always has the same fields as the appended data rows.
  : > "${OUT_CSV}"

  ALLELE_NUM=0
  for ALLELE in "${ALLELES[@]}"; do
    ALLELE_NUM=$((ALLELE_NUM + 1))
    if (( ALLELE_NUM % 10 == 0 )); then
      echo "[INFO]   ${LEN}mer: allele ${ALLELE_NUM}/${N_ALLELES} -- ${ALLELE}"
    fi

    # Build MHCflurry input CSV for this allele
    ALLELE_CSV="${TMPDIR}/ip_${LEN}mer_${ALLELE_NUM}.csv"
    echo "peptide,allele,n_flank,c_flank" > "${ALLELE_CSV}"

    # Strip "-" padding from flanks; write peptide,allele,n_flank,c_flank
    awk -F'\t' -v allele="${ALLELE}" '{
      pep = $1
      nf  = $2; gsub(/-/, "", nf)
      cf  = $3; gsub(/-/, "", cf)
      print pep","allele","nf","cf
    }' "${PEPS_TSV}" >> "${ALLELE_CSV}"

    # Run MHCflurry
    ALLELE_OUT="${TMPDIR}/ip_${LEN}mer_${ALLELE_NUM}_out.csv"
    if ! mhcflurry-predict "${ALLELE_CSV}" --out "${ALLELE_OUT}" 2>/dev/null; then
      echo "[WARN]   ${ALLELE}: mhcflurry-predict failed -- skipping"
      rm -f "${ALLELE_CSV}" "${ALLELE_OUT}"
      continue
    fi

    # Preserve the native header once, then append scored rows.
    if [[ ! -s "${OUT_CSV}" ]]; then
      head -n 1 "${ALLELE_OUT}" > "${OUT_CSV}"
    fi
    tail -n +2 "${ALLELE_OUT}" >> "${OUT_CSV}"
    rm -f "${ALLELE_CSV}" "${ALLELE_OUT}"
  done

  rm -f "${PEPS_TSV}"
  N_ROWS=$(wc -l < "${OUT_CSV}")
  echo "[INFO] ${LEN}mer done: ${N_ROWS} rows -> $(basename ${OUT_CSV})"

done

echo ""
echo "=== MHCflurry predictions complete ==="
# NOTE: Change to "for LEN in 08 09 10 11" to summarise all lengths
for LEN in 09; do
  OUT="${OUTPUT_DIR}/canonical_${LEN}mers_flank_mhcflurry_${RUN_DATE}.csv"
  if [[ -f "${OUT}" ]]; then
    echo "  ${LEN}mer: $(wc -l < "${OUT}") rows -> $(basename ${OUT})"
  else
    echo "  ${LEN}mer: MISSING"
  fi
done
