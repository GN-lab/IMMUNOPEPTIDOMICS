#!/bin/bash
# ============================================================
# JUGNU_IMMUNOPEPTIDOME — config.sh
# July 2026 | Gaurav Raichand | The Institute of Cancer Research
# Single source of truth. Every script sources this.
# ============================================================

export HOME="/data/rds/DMP/UCEC/EVOLIMMU/graichand/fake_home"
export CONDA_PKGS_DIRS="/data/rds/DMP/UCEC/EVOLIMMU/graichand/.conda_cache"
export CONDA_ENVS_DIRS="/data/rds/DMP/UCEC/EVOLIMMU/graichand/.conda_envs"
export CONDARC="/data/rds/DMP/UCEC/EVOLIMMU/graichand/.condarc"
export CONDA_CONFIG_DIR="/data/rds/DMP/UCEC/EVOLIMMU/graichand/.conda_config"
export TMPDIR="/data/rds/DMP/UCEC/EVOLIMMU/graichand/tmp"
mkdir -p "$HOME" "$CONDA_PKGS_DIRS" "$CONDA_ENVS_DIRS" "$CONDA_CONFIG_DIR" "$TMPDIR"

source "/data/scratch/DMP/UCEC/EVOLIMMU/graichand/miniconda3/etc/profile.d/conda.sh"
module load Mamba/23.1.0-0
module load CUDA/12.1.1

export R_LIBS_USER="/data/rds/DMP/UCEC/EVOLIMMU/graichand/R_libs"
conda activate /data/rds/DMP/UCEC/EVOLIMMU/graichand/.conda_envs/neojunction_viz/

# --- Pipeline root ------------------------------------------
export WORKDIR="/data/rds/DMP/UCEC/EVOLIMMU/graichand/IMMUNOPEPTIDOMICS"
export MASTER_SCRIPT_PATH="${WORKDIR}/master.sh"
export SCRIPTS_DIR="${WORKDIR}/scripts"

# --- Input --------------------------------------------------
# UniProt canonical human proteome (Swiss-Prot, NO isoforms)
export FASTA_CANONICAL="${WORKDIR}/data/UP000005640_9606.fasta"

# Standard HLA allele list -- generate once with:
#   bash scripts/setup_alleles.sh


export HLA_ALLELE_FILE="${WORKDIR}/data/iedb_70_alleles.txt"

# --- Tools --------------------------------------------------
export NETMHCPAN_BIN="/data/rds/DMP/UCEC/EVOLIMMU/graichand/netMHCpan-4.2/netMHCpan"

# --- Output -------------------------------------------------
export OUTPUT_DIR="${WORKDIR}/results"
export STEP13_OUTPUT_DIR="${OUTPUT_DIR}"
export STEP14_OUTPUT_DIR="${OUTPUT_DIR}"

# --- Runtime parameters -------------------------------------
export THREADS=24
export CHUNK_SIZE=5000

mkdir -p "${OUTPUT_DIR}" "${WORKDIR}/logs" "${WORKDIR}/data" \
         "${WORKDIR}/.checkpoints" "${WORKDIR}/.slurm_jobs"

echo "[CONFIG] JUGNU_IMMUNOPEPTIDOME loaded from ${BASH_SOURCE[0]}"
echo "[CONFIG] WORKDIR:          ${WORKDIR}"
echo "[CONFIG] FASTA_CANONICAL:  ${FASTA_CANONICAL}"
echo "[CONFIG] HLA_ALLELE_FILE:  ${HLA_ALLELE_FILE}"
echo "[CONFIG] NETMHCPAN_BIN:    ${NETMHCPAN_BIN}"
echo "[CONFIG] OUTPUT_DIR:       ${OUTPUT_DIR}"
