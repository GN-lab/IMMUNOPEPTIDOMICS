
# IMMUNOPEPTIDOMICS

A SLURM-based pipeline for generating a canonical human immunopeptidome by integrating peptide–HLA predictions from **NetMHCpan** and **MHCflurry**.

The pipeline tiles protein-coding Ensembl transcripts into peptides, predicts HLA class I presentation, performs cross-tool concordance analysis, and produces transcript-coordinate BED files of high-confidence candidates.

## Pipeline overview

1. **Peptide generation**
   - Extracts 8–11-mer peptides from protein-coding transcripts.
   - Uses `EnsDb.Hsapiens.v86` to retain accurate ENST-relative amino-acid coordinates.
   - Generates peptide sequences and 30-amino-acid flanking contexts.

2. **NetMHCpan prediction**
   - Scores peptides against a supplied HLA class I allele list.
   - Processes peptides in chunks for large-scale cluster execution.

3. **MHCflurry prediction**
   - Calculates binding affinity and presentation scores using peptide flanking sequences.
   - Supports GPU execution through SLURM.

4. **Cross-tool concordance**
   - Combines NetMHCpan and MHCflurry predictions by peptide and HLA allele.
   - Separates candidates into confidence tiers.

5. **BED generation**
   - Maps high-confidence peptides back to all ENST coordinate occurrences.
   - Collapses all associated HLA alleles into the fifth BED column.

## Concordance tiers

| Tier | NetMHCpan | MHCflurry |
|---|---|---|
| Tier 1: high confidence | EL rank `< 0.5` | Affinity `< 500 nM` |
| Tier 2: medium confidence | EL rank `>= 0.5` and `< 2.0` | Affinity `< 500 nM` |
| Tier 3: discordant | Predicted by only one tool | Non-binder or missing in the other tool |

The final canonical immunopeptidome BED contains **Tier 1 candidates only**.

## Output BED format

```text
ENST_ID    AA_START    AA_END    PEPTIDE    HLA_ALLELES
```

Example:

```text
ENST00000000233    5    14    SALFSRIFG    HLA-A0301,HLA-A1101
```

Coordinates use BED convention:

- `AA_START`: zero-based
- `AA_END`: end-exclusive
- `AA_END - AA_START`: peptide length
- Multiple Tier 1 HLA alleles are retained as a comma-separated list

## Repository structure

```text
IMMUNOPEPTIDOMICS/
├── master.sh
├── config.example.sh
├── setup_alleles.sh
├── scripts/
│   ├── 1_gen_nmers.R
│   ├── 2a_NetMHCPan_analysis.sh
│   ├── 2b_select_top_alleles_NetMHC.R
│   ├── 3a_mhcflurry2_analysis_with_flank.sh
│   ├── 3b_select_top_alleles_MHCflurry.R
│   └── 4_gen_bedfiles.R
├── analysis/
│   └── venn_exact_matches.R
├── data/
├── results/
├── logs/
└── .checkpoints/
```

Generated results, logs, checkpoints, reference FASTA files and local configuration files should not be committed to Git.

## Requirements

### Cluster environment

- Linux
- SLURM
- Bash
- R
- Python/Conda or Mamba
- CUDA-compatible GPU for MHCflurry, if using GPU execution

### Prediction software

- NetMHCpan 4.2
- MHCflurry

NetMHCpan must be downloaded and licensed separately according to its distribution terms.

### R packages

```r
data.table
stringr
ensembldb
EnsDb.Hsapiens.v86
```

Install the CRAN packages with:

```r
install.packages(c("data.table", "stringr"))
```

Install the Bioconductor packages with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install(c(
    "ensembldb",
    "EnsDb.Hsapiens.v86"
))
```

## Configuration

Create a local configuration file:

```bash
cp config.example.sh config.sh
```

Edit `config.sh` and provide:

```bash
WORKDIR="/path/to/IMMUNOPEPTIDOMICS"
NETMHCPAN_BIN="/path/to/netMHCpan"
HLA_ALLELE_FILE="${WORKDIR}/data/hla_alleles.txt"
OUTPUT_DIR="${WORKDIR}/results"
THREADS=24
CHUNK_SIZE=5000
```

The real `config.sh` should remain untracked because it may contain user-specific paths.

Load the configuration:

```bash
source config.sh
```

## HLA allele list

Provide one NetMHCpan-compatible allele per line:

```text
HLA-A01:01
HLA-A02:01
HLA-A03:01
HLA-B07:02
HLA-C07:01
```

Alternatively, generate an HLA-A, HLA-B and HLA-C list from the local NetMHCpan installation:

```bash
bash setup_alleles.sh
```

## Running the pipeline

Run a dry-run to inspect the jobs without submitting them:

```bash
bash master.sh --dry-run
```

Submit the pipeline:

```bash
sbatch master.sh
```

The master script submits dependent SLURM jobs and validates each major output before creating its checkpoint.

## Checkpoints and restarting

Successful steps create files under:

```text
.checkpoints/
```

Example:

```text
.checkpoints/step1_gen_nmers.done
.checkpoints/step2a_netmhcpan.done
.checkpoints/step2b_select_netmhc.done
.checkpoints/step3a_mhcflurry.done
.checkpoints/step3b_select_mhcflurry.done
.checkpoints/step4_gen_bedfiles.done
```

Rerunning `master.sh` skips completed steps.

To rerun a particular step and its downstream steps, remove only the relevant checkpoint files after confirming that this is intended.

## Main outputs

### Peptide coordinate tables

```text
canonical_09mers_YYYYMMDD.tsv
```

Columns:

```text
enst
aa_start
aa_end
n_mer
ctex_up
ctex_dn
```

These intermediate coordinates are one-based and inclusive.

### Concordance tables

```text
canonical_concordance_tier1_YYYYMMDD.tsv
canonical_concordance_tier2_YYYYMMDD.tsv
canonical_concordance_tier3_YYYYMMDD.tsv
canonical_concordance_all_YYYYMMDD.tsv
```

### Final high-confidence immunopeptidome

```text
canonical_immunopeptidome_SB_YYYYMMDD.tsv
canonical_immunopeptidome_SB_YYYYMMDD.bed
```

The BED file contains Tier 1 candidates with all associated strong-binding HLA alleles retained in column 5.

## Candidate-level comparisons

For comparisons where HLA identity is not required, use the first four BED columns:

```bash
cut -f1-4 input.bed |
LC_ALL=C sort -u \
> candidates.bed
```

The comparison key is then:

```text
ENST + AA start + AA end + peptide
```

Exact candidate matches can be identified with BEDTools:

```bash
bedtools intersect \
  -sorted \
  -a dataset_A.sorted.bed \
  -b dataset_B.sorted.bed \
  -f 1.0 -r \
  -wa -wb |
awk 'BEGIN {FS=OFS="\t"}
     $1==$5 && $2==$6 && $3==$7 && $4==$8 {
         print $1,$2,$3,$4
     }' |
LC_ALL=C sort -u \
> exact_candidate_matches.bed
```

## Notes

- Peptide tiling supports 8–11-mers.
- The current prediction workflow is configured primarily for 9-mers.
- Large generated files are intentionally excluded from the repository.
- Reference proteomes and licensed prediction software must be obtained separately.

## Author

**Gaurav Raichand**  
The Institute of Cancer Research
