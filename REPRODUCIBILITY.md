# Reproducibility Guide

## Research overview

This repository investigates the effects of fragment length and post-mortem damage (PMD) on mapping ancient-human-mitochondrial-DNA-like reads near artificial linear-reference boundaries.

The study compares BWA-aln and Bowtie2 using the original human mitochondrial revised Cambridge Reference Sequence (rCRS) and a circularly rotated representation shifted by 8,000 bp.

## 1. Software environment

The experiments were conducted under Linux using a Conda environment.

Key software versions:

| Software | Version |
|---|---|
| Python | 3.11 |
| BWA | 0.7.19-r1273 |
| Bowtie2 | 2.5.5 |
| SAMtools | 1.24 |
| Gargammel | 1.1.4 |
| pysam | 0.24.0 |
| NumPy | 2.4.6 |
| pandas | 3.1.0 |
| Matplotlib | 3.11.2 |

To create the portable Conda environment:

```bash
conda env create -f environment.yml
conda activate adna_mapping
```

The environment.yml file records the principal dependencies, without pinning every package to its original build.

For stronger software-environment reproducibility on compatible linux-64 systems, the repository also includes environment-linux-64.txt, an explicit Conda package specification.

The explicit specification can be installed using:

```bash
conda create -n adna_mapping_exact --file environment-linux-64.txt
conda activate adna_mapping_exact
```

The explicit package specification is platform-specific and may depend on continued availability of the original Conda packages.

## 2. Reference sequences

The original mitochondrial reference is:

reference/rCRS_NC_012920.1.fasta

Accession: NC_012920.1
Reference length: 16,569 bp

The shifted reference is:

reference/rCRS_shifted_8000.fasta

Its sequence is an exact 8,000-base rotation of the original:

shifted = original[8000:] + original[:8000]

Both reference representations contain the same nucleotide sequence in different linear coordinate systems.

Generated aligner indexes are excluded from Git and can be rebuilt:

```bash
bash scripts/build_reference_indexes.sh
```

Run all commands in this guide from the repository root.

## 3. Experimental design

Fragment lengths:

- 30 bp
- 40 bp
- 50 bp
- 70 bp

Sampling groups:

- Original-reference boundary
- Shifted-reference boundary
- Interior controls

Each length contains L-1 genomic fragments per group, represented on both forward and reverse strands.

Across the four fragment lengths, the baseline contains 1,116 reads.

Baseline read selection uses deterministic seeds:

seed = 1601 + fragment_length

The original and shifted boundary groups are enumerated from their genomic positions. Interior fragments are sampled without replacement.

## 4. Generate undamaged reads and ground truth

```bash
python scripts/generate_circular_length_baseline.py
```

This creates baseline FASTA files in:

simulation/circular_length_baseline/

And ground-truth TSV files in:

ground_truth/circular_length_baseline/

This generator writes its outputs directly. Do not run it over the archived files unless regeneration is intended.

## 5. Generate reference alignments

After building the reference indexes:

```bash
bash scripts/run_circular_length_baseline_alignments.sh
```

The script processes 16 undamaged alignment configurations.

Alignment parameters:

BWA-aln:

```bash
bwa aln -n 0.04
bwa samse
```

Bowtie2:

```bash
bowtie2 --end-to-end --sensitive -f
```

Existing SAM files are reused only if they pass primary-record count validation. Invalid existing SAM files are not overwritten.

Expected baseline total:

- 16 alignment configurations
- 4,464 primary alignment records

## 6. Generate PMD datasets

Three PMD parameter settings were used:

| Severity | deamSim -damage parameters |
|---|---|
| Low | 0.03,0.30,0.005,0.05 |
| Moderate | 0.03,0.30,0.01,0.10 |
| High | 0.03,0.30,0.03,0.30 |

The original PMD realization uses the parameters and seeds recorded in:

simulation/circular_pmd_main/damage_manifest.tsv

Generate or validate these 12 datasets:

```bash
python scripts/generate_circular_pmd_main.py
```

Generate or validate the additional 108 datasets:

```bash
python scripts/generate_circular_pmd_replicates.py
```

The latter script records its configuration in:

simulation/circular_pmd_replicates/damage_replicates_manifest.tsv

For replicates 2 through 10, the seed formula is:

seed = 20000 + 10000 * replicate + 10 * length + level_code

Where level_code is 1 for low, 2 for moderate, and 3 for high.

All ten PMD realizations use the same underlying genomic fragments. They represent independent stochastic damage realizations, not independent biological replicates.

The simulation validation permits C-to-T and G-to-A substitutions without insertions or deletions.

Existing damaged FASTA files are validated and preserved rather than automatically overwritten.

## 7. Generate PMD alignments

Main PMD realization:

```bash
bash scripts/run_circular_pmd_alignments.sh
```

Expected:

- 48 alignment configurations
- 13,392 primary alignment records

Additional replicates 2 through 10:

```bash
bash scripts/run_circular_pmd_replicate_alignments.sh
```

Expected:

- 432 alignment configurations
- 120,528 primary alignment records

Across all PMD realizations:

- 480 alignment configurations
- 133,920 primary alignment records

The alignment scripts use the same BWA-aln and Bowtie2 parameters as the undamaged baseline.

The main PMD alignment script writes its outputs directly. Regenerating that stage may overwrite existing alignment outputs. Use a separate working copy for a clean rerun.

## 8. Evaluate alignments

Run the following commands in order:

```bash
python scripts/evaluate_circular_length_baseline.py
python scripts/evaluate_circular_pmd_all.py
python scripts/evaluate_circular_pmd_replicates.py
python scripts/summarize_circular_pmd_replicates.py
python scripts/compare_circular_reference_shifting.py
python scripts/evaluate_circular_mapq_filtering.py
```

The principal evaluation outputs are stored in results/.

The reference-shifting comparison uses the replicate-level recall table produced by the summary script.

The MAPQ analysis examines alignment files separately.

## 9. Generate final manuscript figures

```bash
python scripts/redesign_manuscript_figures.py
```

This script first executes generate_all_manuscript_figures.py, then generates the final redesigned manuscript figures.

The six figures are saved in figures/ as PNG and PDF files:

1. figure1_pmd_length_recall
2. figure2_reference_boundary
3. figure3_shifting_effect
4. figure4_damage_baseline
5. figure5_mapping_errors
6. figure6_mapq_filtering

Running this command overwrites figure files with matching names. Use a separate working copy if the archived figures must remain unchanged.

## 10. Reproducibility notes

- The underlying genomic fragments are fixed across PMD replicates.
- PMD replicates are stochastic simulation replicates, not biological replicates.
- Mapping is evaluated against known synthetic ground-truth origins.
- The reference shift relocates the artificial linear boundary; it does not eliminate the boundary.
- Primary-record counts exclude secondary and supplementary SAM alignments using samtools view -c -F 0x900.
- Successful record-count validation does not independently establish mapping correctness.
- The repository preserves experimental outputs for inspection.
- Exact regenerated outputs may depend on software builds and deterministic behavior of the underlying tools.

## 11. Repository conventions

- reference/: original and shifted mitochondrial references
- simulation/: undamaged and PMD-simulated reads
- ground_truth/: known read origins
- alignments/: alignment files and intermediates
- results/: evaluation tables and summaries
- figures/: manuscript figures
- logs/: recorded experiment logs
- scripts/: simulation, alignment, evaluation, and plotting scripts

Generated BWA, Bowtie2, and FASTA index files are excluded from Git.

A local experimental copy of Gargammel source is also excluded. The documented PMD generation pipeline invokes the Conda-installed deamSim executable.

## 12. Scope

The project uses controlled synthetic mitochondrial-DNA-like reads and known ground truth. Its findings should not be interpreted as direct estimates of performance across all authentic ancient-DNA samples, sequencing protocols, or reference genomes.
