#!/bin/bash
#SBATCH --job-name=IMMUNOPEP_JUGNU
#SBATCH --partition=compute
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=8042
#SBATCH --time=5-00:00:00
#SBATCH --output=logs/master_%j.log

set -euo pipefail
set -E
trap 'rc=$?; echo "**ERROR** rc=$rc line $LINENO: $BASH_COMMAND" >&2; exit $rc' ERR

# ============================================================
# JUGNU_IMMUNOPEPTIDOME — master.sh
# July 2026 | Gaurav Raichand | The Institute of Cancer Research
#
# Usage: sbatch master.sh          (from WORKDIR)
#        bash master.sh --dry-run  (print jobs without submitting)
#
# Pipeline steps:
#   1  Gen_nmers.R                    -- tile canonical FASTA -> 8-11mer TSVs
#   2a NetMHCPan_analysis.sh          -- score all peptides against HLA alleles
#   2b select_top_alleles_NetMHC.R    -- filter SB, best allele per peptide
#   3a mhcflurry2_analysis.sh (GPU)   -- MHCflurry presentation scores
#   3b select_top_alleles_MHCflurry.R -- filter and rank by presentation score
#   4  gen_bedfiles.R                 -- bed files, matrix, figures
#
# Each step:
#   - Is submitted as a durable SLURM job (not --wrap) so restarts are safe
#   - Validates its key output before touching the checkpoint
#   - Is skipped automatically if checkpoint already exists (safe re-run)
#   - Depends on the previous step via SLURM --dependency=afterok
# ============================================================

# WORKDIR must be hardcoded -- BASH_SOURCE[0] resolves to the SLURM spool
# path (/tmp/slurmd/job.../slurm_script) when submitted via sbatch, not the
# real file location. Keeping this hardcoded is intentional.
WORKDIR="/data/rds/DMP/UCEC/EVOLIMMU/graichand/IMMUNOPEPTIDOMICS"
DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

# MASTER_SCRIPT_PATH: durable path so child jobs can call back for validation.
# SLURM stages $0 under /tmp which gets swept -- read from config.sh instead.
if [[ -z "${MASTER_SCRIPT_PATH:-}" ]]; then
  [[ -f "${WORKDIR}/config.sh" ]] && source "${WORKDIR}/config.sh"
fi
SCRIPT_PATH="${MASTER_SCRIPT_PATH:-${WORKDIR}/master.sh}"
if [[ ! -f "$SCRIPT_PATH" ]]; then
  echo "[ERROR] master.sh not found at ${SCRIPT_PATH}" >&2
  echo "[ERROR] Set MASTER_SCRIPT_PATH in config.sh to the full path of this file." >&2
  exit 1
fi

###############################################################################
#  SELF-DISPATCH: validation callbacks from child jobs
###############################################################################

_v_single() {
  local step_name="$1" search_dir="$2" pattern="$3" min_lines="${4:-2}"
  if [[ "$search_dir" == "SKIP" || "$pattern" == "SKIP" ]]; then
    echo "[VALIDATE SKIP] $step_name"; return 0
  fi
  if [[ ! -d "$search_dir" ]]; then
    echo "[VALIDATE FAIL] $step_name: dir missing: $search_dir" >&2; return 1
  fi
  local match
  match=$(find "$search_dir" -maxdepth 1 -type f \
    -regextype posix-extended -regex ".*/${pattern}" 2>/dev/null | sort | tail -n 1)
  if [[ -z "$match" ]]; then
    echo "[VALIDATE FAIL] $step_name: no file matching '${pattern}' in ${search_dir}" >&2
    return 1
  fi
  local size
  size=$(stat -c%s "$match" 2>/dev/null || echo 0)
  if [[ "$size" -lt 50 ]]; then
    echo "[VALIDATE FAIL] $step_name: $match is empty (${size} bytes)" >&2; return 1
  fi
  case "$match" in
    *.tsv|*.csv|*.txt)
      local nlines
      nlines=$(wc -l < "$match")
      if [[ "$nlines" -lt "$min_lines" ]]; then
        echo "[VALIDATE FAIL] $step_name: $match has only $nlines line(s)" >&2; return 1
      fi
      ;;
  esac
  echo "[VALIDATE OK] $step_name: $match (${size} bytes)"
  return 0
}

if [[ "${1:-}" == "__validate_step__" ]]; then
  shift; source "${WORKDIR}/config.sh"; _v_single "$@"; exit $?
fi

###############################################################################
#  NORMAL MODE
###############################################################################

source "${WORKDIR}/config.sh"
module load R

# Validate inputs exist before submitting anything
if [[ ! -f "${FASTA_CANONICAL}" ]]; then
  echo "[ERROR] FASTA not found: ${FASTA_CANONICAL}" >&2
  echo "[ERROR] Download UP000005640_9606.fasta into ${WORKDIR}/data/" >&2
  exit 1
fi
if [[ ! -f "${HLA_ALLELE_FILE}" ]] || [[ ! -s "${HLA_ALLELE_FILE}" ]]; then
  echo "[ERROR] Standard HLA allele file not found: ${HLA_ALLELE_FILE}" >&2
  echo "[ERROR] Generate it once with: bash scripts/setup_alleles.sh" >&2
  exit 1
fi

mkdir -p "${WORKDIR}/.checkpoints" "${WORKDIR}/logs" "${WORKDIR}/.slurm_jobs"

JOBSCRIPT_DIR="${WORKDIR}/.slurm_jobs"
step_done()  { [[ -f "${WORKDIR}/.checkpoints/$1.done" ]]; }

write_job_script() {
  local name="$1" body="$2"
  local path="${JOBSCRIPT_DIR}/${name}.sh"
  printf '#!/bin/bash\nset -euo pipefail\nmodule load R\n%s\n' "$body" > "$path"
  chmod +x "$path"
  echo "$path"
}

submit_job() {
  # submit_job NAME COMMAND DEP VAL_DIR VAL_PATTERN VAL_MINLINES EXTRA_SBATCH_ARGS
  local name="$1" command="$2" dep="${3:-}"
  local val_dir="${4:-SKIP}" val_pat="${5:-SKIP}" val_min="${6:-2}"
  local extra="${7:-}"

  local body="${command}
bash \"${SCRIPT_PATH}\" __validate_step__ \"${name}\" \"${val_dir}\" \"${val_pat}\" \"${val_min}\"
touch \"${WORKDIR}/.checkpoints/${name}.done\""

  local script
  script=$(write_job_script "$name" "$body")

  local dep_flag=""
  [[ -n "$dep" ]] && dep_flag="--dependency=afterok:${dep}"

  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] Would submit: $name (dep=${dep:-none}) script=$script"
    echo "dry_run_${name}"
    return
  fi

  local jobid
  jobid=$(sbatch --parsable \
    --job-name="$name" \
    --partition=compute \
    --cpus-per-task=4 \
    --mem-per-cpu=8042 \
    --time=12:00:00 \
    --output="logs/${name}_%j.log" \
    $dep_flag $extra \
    "$script")
  echo "  Submitted ${name} -> jobid ${jobid} (script: $script)" >&2
  echo "$jobid"
}

jobPrev=""

echo "=== JUGNU_IMMUNOPEPTIDOME Pipeline Submission ==="
echo "    WORKDIR:  ${WORKDIR}"
echo "    FASTA:    ${FASTA_CANONICAL}"
echo "    HLA:      ${HLA_ALLELE_FILE}"
echo "    OUTPUT:   ${OUTPUT_DIR}"
echo ""

# -----------------------------------------------------------------------
# Step 1: Gen_nmers.R
#   Resources: compute, 16 CPUs, 8042MB/cpu (~128G), 4h
#   Input:  FASTA_CANONICAL
#   Output: canonical_0Xmers_YYYYMMDD.tsv (x4)
# -----------------------------------------------------------------------
if ! step_done "step1_gen_nmers"; then
  BODY_1="source \"${WORKDIR}/config.sh\"
Rscript \"${SCRIPTS_DIR}/1_gen_nmers.R\"
bash \"${SCRIPT_PATH}\" __validate_step__ step1_gen_nmers \"${OUTPUT_DIR}\" 'canonical_09mers_[0-9]{8}\\.tsv' 2
touch \"${WORKDIR}/.checkpoints/step1_gen_nmers.done\""
  SCRIPT_1=$(write_job_script "step1_gen_nmers" "$BODY_1")
  echo "[INFO] Submitting Step 1: Gen_nmers.R (16 CPUs, ~128G, 4h)..."
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] step1_gen_nmers"; jobPrev="dry_run_1"
  else
    jobPrev=$(sbatch --parsable \
      --job-name="step1_gen_nmers" \
      --partition=compute \
      --cpus-per-task=16 \
      --mem-per-cpu=8042 \
      --time=4:00:00 \
      --output="logs/step1_gen_nmers_%j.log" \
      ${jobPrev:+--dependency=afterok:${jobPrev}} \
      "$SCRIPT_1")
    echo "  Submitted step1_gen_nmers -> jobid ${jobPrev}"
  fi
else
  echo "  step1_gen_nmers: already done (checkpoint exists)"
fi

# -----------------------------------------------------------------------
# Step 2a: NetMHCpan predictions
#   Resources: compute, 24 CPUs, 8042MB/cpu (~192G), 5 days
#   Input:  canonical_0Xmers_*.tsv
#   Output: canonical_0Xmer_netmhcpan_YYYY_MMDD.tsv (x4)
# -----------------------------------------------------------------------
if ! step_done "step2a_netmhcpan"; then
  BODY_2A="source \"${WORKDIR}/config.sh\"
bash \"${SCRIPTS_DIR}/2a_NetMHCPan_analysis.sh\"
bash \"${SCRIPT_PATH}\" __validate_step__ step2a_netmhcpan \"${OUTPUT_DIR}\" 'canonical_09mer_netmhcpan_[0-9_]+\\.tsv' 2
touch \"${WORKDIR}/.checkpoints/step2a_netmhcpan.done\""
  SCRIPT_2A=$(write_job_script "step2a_netmhcpan" "$BODY_2A")
  echo "[INFO] Submitting Step 2a: NetMHCpan (24 CPUs, ~192G, 5 days)..."
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] step2a_netmhcpan"; jobPrev="dry_run_2a"
  else
    jobPrev=$(sbatch --parsable \
      --job-name="step2a_netmhcpan" \
      --partition=compute \
      --cpus-per-task=24 \
      --mem-per-cpu=8042 \
      --time=5-00:00:00 \
      --output="logs/step2a_netmhcpan_%j.log" \
      ${jobPrev:+--dependency=afterok:${jobPrev}} \
      "$SCRIPT_2A")
    echo "  Submitted step2a_netmhcpan -> jobid ${jobPrev}"
  fi
else
  echo "  step2a_netmhcpan: already done"
fi

# -----------------------------------------------------------------------
# Step 2b: Select top alleles (NetMHCpan)
#   Resources: compute, 8 CPUs, ~64G, 2h
#   Input:  canonical_0Xmer_netmhcpan_*.tsv
#   Output: canonical_0Xmer_netmhcpan_selected_alleles_YYYYMMDD.tsv (x4)
# -----------------------------------------------------------------------
if ! step_done "step2b_select_netmhc"; then
  echo "[INFO] Submitting Step 2b: Select top NetMHCpan alleles (8 CPUs, ~64G, 2h)..."
  BODY_2B="source \"${WORKDIR}/config.sh\"
Rscript \"${SCRIPTS_DIR}/2b_select_top_alleles_NetMHC.R\"
bash \"${SCRIPT_PATH}\" __validate_step__ step2b_select_netmhc \"${OUTPUT_DIR}\" 'canonical_09mer_netmhcpan_selected_alleles_[0-9]{8}\\.tsv' 2
touch \"${WORKDIR}/.checkpoints/step2b_select_netmhc.done\""
  SCRIPT_2B=$(write_job_script "step2b_select_netmhc" "$BODY_2B")
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] step2b_select_netmhc"; jobPrev="dry_run_2b"
  else
    jobPrev=$(sbatch --parsable \
      --job-name="step2b_select_netmhc" \
      --partition=compute \
      --cpus-per-task=8 \
      --mem-per-cpu=8042 \
      --time=2:00:00 \
      --output="logs/step2b_select_netmhc_%j.log" \
      ${jobPrev:+--dependency=afterok:${jobPrev}} \
      "$SCRIPT_2B")
    echo "  Submitted step2b_select_netmhc -> jobid ${jobPrev}"
  fi
else
  echo "  step2b_select_netmhc: already done"
fi

# -----------------------------------------------------------------------
# Step 3a: MHCflurry predictions (GPU)
#   Resources: gpu, 4 CPUs, 13600MB/cpu, 1 GPU, 7 days
#   Input:  canonical_0Xmer_mhcflurry_input_*.csv (written by 2b prep step)
#   Output: canonical_0Xmers_flank_mhcflurry_*.csv (x4)
# -----------------------------------------------------------------------
if ! step_done "step3a_mhcflurry"; then
  BODY_3A="source \"${WORKDIR}/config.sh\"
bash \"${SCRIPTS_DIR}/3a_mhcflurry2_analysis_with_flank.sh\"
bash \"${SCRIPT_PATH}\" __validate_step__ step3a_mhcflurry \"${OUTPUT_DIR}\" 'canonical_09mers_flank_mhcflurry_[0-9_]+\\.csv' 2
touch \"${WORKDIR}/.checkpoints/step3a_mhcflurry.done\""
  SCRIPT_3A=$(write_job_script "step3a_mhcflurry" "$BODY_3A")
  echo "[INFO] Submitting Step 3a: MHCflurry GPU (4 CPUs, 1 GPU, 7 days)..."
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] step3a_mhcflurry"; jobPrev="dry_run_3a"
  else
    jobPrev=$(sbatch --parsable \
      --job-name="step3a_mhcflurry" \
      --partition=gpu \
      --cpus-per-task=4 \
      --mem-per-cpu=13600 \
      --gpus-per-node=1 \
      --time=7-00:00:00 \
      --output="logs/step3a_mhcflurry_%j.log" \
      ${jobPrev:+--dependency=afterok:${jobPrev}} \
      "$SCRIPT_3A")
    echo "  Submitted step3a_mhcflurry -> jobid ${jobPrev}"
  fi
else
  echo "  step3a_mhcflurry: already done"
fi

# -----------------------------------------------------------------------
# Step 3b: Select top alleles (MHCflurry)
#   Resources: compute, 8 CPUs, ~64G, 2h
#   Input:  canonical_0Xmers_flank_mhcflurry_*.csv
#   Output: canonical_0Xmer_mhcflurry_top_YYYYMMDD.tsv (x4)
#           canonical_0Xmer_mhcflurry_per_sample_YYYYMMDD.tsv (x4)
# -----------------------------------------------------------------------
if ! step_done "step3b_select_mhcflurry"; then
  BODY_3B="source \"${WORKDIR}/config.sh\"
Rscript \"${SCRIPTS_DIR}/3b_select_top_alleles_MHCflurry.R\"
bash \"${SCRIPT_PATH}\" __validate_step__ step3b_select_mhcflurry \"${OUTPUT_DIR}\" 'canonical_09mer_mhcflurry_top_[0-9]{8}\\.tsv' 2
touch \"${WORKDIR}/.checkpoints/step3b_select_mhcflurry.done\""
  SCRIPT_3B=$(write_job_script "step3b_select_mhcflurry" "$BODY_3B")
  echo "[INFO] Submitting Step 3b: Select top MHCflurry alleles (8 CPUs, ~64G, 2h)..."
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] step3b_select_mhcflurry"; jobPrev="dry_run_3b"
  else
    jobPrev=$(sbatch --parsable \
      --job-name="step3b_select_mhcflurry" \
      --partition=compute \
      --cpus-per-task=8 \
      --mem-per-cpu=8042 \
      --time=2:00:00 \
      --output="logs/step3b_select_mhcflurry_%j.log" \
      ${jobPrev:+--dependency=afterok:${jobPrev}} \
      "$SCRIPT_3B")
    echo "  Submitted step3b_select_mhcflurry -> jobid ${jobPrev}"
  fi
else
  echo "  step3b_select_mhcflurry: already done"
fi

# -----------------------------------------------------------------------
# Step 4: Filter SB peptides, join coordinates, generate bed file
#   Resources: compute, 8 CPUs, ~64G, 2h
#   Input:  canonical_*_top_*.tsv, canonical_*_per_sample_*.tsv
#   Output: canonical_immunopeptidome_SB_YYYYMMDD.bed  +  .tsv
# -----------------------------------------------------------------------
if ! step_done "step4_gen_bedfiles"; then
  BODY_4="source \"${WORKDIR}/config.sh\"
Rscript \"${SCRIPTS_DIR}/4_gen_bedfiles.R\"
bash \"${SCRIPT_PATH}\" __validate_step__ step4_gen_bedfiles \"${OUTPUT_DIR}\" 'canonical_immunopeptidome_SB_[0-9]{8}\\.tsv' 2
touch \"${WORKDIR}/.checkpoints/step4_gen_bedfiles.done\""
  SCRIPT_4=$(write_job_script "step4_gen_bedfiles" "$BODY_4")
  echo "[INFO] Submitting Step 4: Bed files + matrix + figures (8 CPUs, ~64G, 2h)..."
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] step4_gen_bedfiles"; jobPrev="dry_run_4"
  else
    jobPrev=$(sbatch --parsable \
      --job-name="step4_gen_bedfiles" \
      --partition=compute \
      --cpus-per-task=8 \
      --mem-per-cpu=8042 \
      --time=2:00:00 \
      --output="logs/step4_gen_bedfiles_%j.log" \
      ${jobPrev:+--dependency=afterok:${jobPrev}} \
      "$SCRIPT_4")
    echo "  Submitted step4_gen_bedfiles -> jobid ${jobPrev}"
  fi
else
  echo "  step4_gen_bedfiles: already done"
fi

echo ""
echo "=== Submission complete ==="
echo "    Chain: step1 -> step2a -> step2b -> step3a -> step3b -> step4"
echo "    Monitor:  squeue -u \$(whoami)"
echo "    Re-run:   delete .checkpoints/STEPNAME.done, then sbatch master.sh"
echo "    Dry-run:  bash master.sh --dry-run"
