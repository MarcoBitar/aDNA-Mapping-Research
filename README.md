# Ancient Human Mitochondrial DNA Mapping: Fragment Length, Post-Mortem Damage, and Reference-Boundary Effects

### A Comparative Computational Study Using BWA-aln and Bowtie2

**Authors:** Marco Bitar and Matéo Bitar  
**Institution:** American University of Beirut (AUB)  
**Course:** CMPS 355, Bioinformatics Algorithms  
**Academic Year:** Fall 2026  
**Project Type:** Graduate Research Project

---

## Overview

Ancient DNA (aDNA) provides valuable information about human ancestry, population history, migration, and evolution. However, ancient DNA fragments are typically short and frequently contain post-mortem nucleotide damage, making accurate alignment to reference genomes challenging.

These difficulties are particularly relevant to mitochondrial DNA (mtDNA), which is biologically circular but is commonly represented as a linear sequence by conventional read aligners. This linear representation introduces an artificial boundary at which reads spanning the circular junction cannot be represented as ordinary contiguous alignments.

This research investigates how fragment length, post-mortem damage (PMD), reference representation, and aligner choice affect the mapping of simulated ancient human mitochondrial DNA reads.

Using BWA-aln and Bowtie2, we conduct controlled experiments against the standard human mitochondrial revised Cambridge Reference Sequence (rCRS) and an alternative representation shifted by 8,000 nucleotides.

The project evaluates whether shifting the reference restores correct mapping at the original boundary, how the resulting boundary relocation affects other genomic regions, and how these effects vary with fragment length and damage severity.

## Research Question

**How do fragment length and post-mortem damage affect the mapping accuracy and sensitivity of ancient human mitochondrial DNA reads near the artificial boundary of a linearized reference genome, and to what extent can reference shifting mitigate these effects in BWA-aln and Bowtie2?**

## Research Objectives

1. Quantify mapping failures affecting reads that span an artificial mitochondrial reference boundary.
2. Compare BWA-aln and Bowtie2 under controlled fragment-length and post-mortem damage conditions.
3. Evaluate reference shifting as a strategy for recovering reads spanning the original rCRS boundary.
4. Determine whether shifting the reference introduces equivalent mapping failures at the relocated boundary.
5. Examine incorrect alignments, alignment operations, and mapping-quality scores.
6. Provide reproducible datasets, scripts, ground-truth annotations, alignment outputs, and experimental results.

## Experimental Design

### Reference Genome

The experiments use the human mitochondrial revised Cambridge Reference Sequence (rCRS), accession **NC_012920.1**, containing 16,569 base pairs.

Two reference representations are evaluated:

| Reference | Description |
|---|---|
| Original rCRS | Standard linearized representation of the mitochondrial genome |
| Shifted rCRS | The same circular sequence rotated by 8,000 nucleotides |

Reference shifting preserves the underlying mitochondrial sequence while changing the location of the artificial linearization boundary.

### Fragment Lengths

Four read lengths are evaluated:

| Fragment Length | Distinct Baseline Reads |
|---|---:|
| 30 bp | 174 |
| 40 bp | 234 |
| 50 bp | 294 |
| 70 bp | 414 |
| **Total** | **1,116** |

For each length, reads are generated from three genomic groups:

- **Original-boundary reads:** Fragments crossing the artificial boundary of the standard rCRS reference.
- **Shifted-boundary reads:** Fragments crossing the artificial boundary introduced by the shifted reference.
- **Interior reads:** Control fragments that do not cross either reference boundary.

Both forward and reverse strand orientations are represented.

### Post-Mortem Damage

The study includes an undamaged baseline and three levels of simulated post-mortem damage:

- Low PMD
- Moderate PMD
- High PMD

Damage simulation models nucleotide substitutions associated with ancient DNA degradation, including C-to-T and G-to-A changes.

Ten independently seeded damage realizations are evaluated for each damaged fragment-length condition. The underlying genomic fragment positions remain fixed across these realizations, allowing the effects of stochastic damage to be assessed without introducing changes in fragment sampling.

### Alignment Methods

| Tool | Version | Configuration |
|---|---|---|
| BWA-aln | 0.7.19 | `bwa aln -n 0.04`, followed by `bwa samse` |
| Bowtie2 | 2.5.5 | `--end-to-end --sensitive -f` |
| SAMtools | 1.24 | Alignment-file inspection and processing |

Both aligners are evaluated against both reference representations under the same experimental conditions.

### Evaluation Metrics

Each simulated read has a known genomic origin, allowing alignment outcomes to be compared with ground truth.

The evaluation distinguishes:

- **Correctly mapped reads:** Primary alignments satisfying the exact genomic placement, strand, and full-length CIGAR criterion.
- **Incorrectly mapped reads:** Reported primary alignments that do not satisfy the predefined correctness criterion.
- **Unmapped reads:** Reads without a reported mapped primary alignment.

The principal metric is strict mapping recall:

\[
\text{Strict Recall} =
\frac{\text{Correctly Mapped Reads}}{\text{Total Simulated Reads}}
\]

Additional analyses consider mapping precision, incorrect-mapping counts, unmapped reads, CIGAR operations, and mapping-quality (MAPQ) filtering.

For shifted-reference alignments, reported coordinates are translated back to the original mitochondrial coordinate system before evaluation.

## Experimental Scope

The completed study includes:

| Component | Total |
|---|---:|
| Fragment lengths | 4 |
| Reference representations | 2 |
| Aligners | 2 |
| Damaged-read realizations per condition | 10 |
| Baseline alignment configurations | 16 |
| Damaged alignment configurations | 480 |
| Total alignment configurations | 496 |
| Total primary alignment records evaluated | 138,384 |

## Main Findings

### 1. Boundary-Associated Mapping Failures

In the undamaged baseline, both aligners achieved 0% strict recall for reads spanning the active artificial reference boundary.

For reads that did not cross the active boundary, strict recall reached 100% in the undamaged baseline.

### 2. Recovery Through Reference Shifting

Shifting the reference restored high strict recall for reads spanning the original rCRS boundary.

Under high PMD, mean recall across ten damage realizations was:

| Fragment Length | BWA-aln | Bowtie2 |
|---|---:|---:|
| 30 bp | 98.97% | 95.17% |
| 40 bp | 97.69% | 96.92% |
| 50 bp | 98.78% | 98.88% |
| 70 bp | 99.35% | 99.57% |

However, shifting also relocated the artificial boundary. Reads spanning the newly introduced boundary exhibited corresponding mapping failures.

Thus, reference shifting changed which genomic fragments were affected rather than eliminating the underlying linearization limitation.

### 3. Fragment Length and Damage Severity

Increasing post-mortem damage generally reduced strict recall for reads that could otherwise be mapped across a contiguous reference region.

The magnitude of this reduction varied with fragment length and aligner.

The largest observed mean recall loss relative to the undamaged baseline was 4.83 percentage points for Bowtie2 at 30 bp under high PMD.

### 4. Different Alignment Failure Modes

Under high PMD, Bowtie2 produced 721 incorrectly mapped active-boundary alignments across both reference representations and ten damage realizations.

All 721 contained insertion operations in their CIGAR strings, despite the absence of true insertions or deletions in the simulated reads.

BWA-aln produced no incorrectly mapped active-boundary alignments in these conditions, instead reporting the corresponding reads as unmapped.

### 5. Mapping-Quality Filtering

For Bowtie2 on the original reference under high PMD, applying a minimum MAPQ threshold of 10 reduced the 408 incorrectly mapped original-boundary alignments to one.

A threshold of 30 eliminated all 408 observed incorrect alignments but also excluded more correctly mapped interior reads.

These findings document a tradeoff between suppressing incorrect alignments and retaining correctly mapped fragments.

## Repository Structure

```text
aDNA_Mapping_Research/
├── scripts/
│   ├── generate_circular_length_baseline.py
│   ├── generate_circular_pmd_replicates.py
│   ├── run_circular_pmd_alignments.sh
│   ├── run_circular_pmd_replicate_alignments.sh
│   ├── evaluate_circular_pmd_replicates.py
│   ├── summarize_circular_pmd_replicates.py
│   ├── compare_circular_reference_shifting.py
│   ├── evaluate_circular_mapq_filtering.py
│   └── redesign_manuscript_figures.py
├── reference/
├── simulation/
│   ├── circular_length_baseline/
│   ├── circular_pmd_main/
│   └── circular_pmd_replicates/
├── ground_truth/
├── alignments/
├── results/
├── figures/
└── logs/
```

The repository also contains calibration scripts, validation scripts, pilot experiments, and intermediate outputs used during development and verification.

## Reproducibility

The project was developed and evaluated in a Linux environment using WSL2 Ubuntu and a Conda environment named `adna_mapping`.

The main software dependencies include:

- Python 3
- BWA 0.7.19
- Bowtie2 2.5.5
- SAMtools 1.24
- Gargammel-based ancient DNA damage simulation
- pysam
- NumPy, pandas, and Matplotlib for analysis and visualization

### Experimental Workflow

The computational workflow consists of the following stages:

1. Prepare the original and shifted mitochondrial reference sequences.
2. Generate controlled boundary-spanning and interior fragments.
3. Generate the undamaged baseline and PMD-affected read datasets.
4. Align each dataset against both references using BWA-aln and Bowtie2.
5. Evaluate primary alignments against the known genomic origins.
6. Aggregate results across independently seeded damage realizations.
7. Compare aligners, fragment lengths, damage levels, and reference representations.
8. Generate manuscript figures and comparative result tables.

The repository preserves the generated experimental inputs and outputs so that the reported findings can be independently inspected.

## Results and Figures

The `results/` directory contains evaluation outputs in CSV and TSV formats, including:

- `circular_length_baseline.csv`
- `circular_pmd_replicates_all.csv`
- `circular_pmd_replicate_recall.csv`
- `circular_pmd_replicate_summary.csv`
- `circular_reference_shifting_comparison.csv`

The `figures/` directory contains the final manuscript visualizations, including analyses of fragment length, PMD severity, reference shifting, incorrect alignments, and MAPQ filtering.

## Scope and Limitations

This research uses simulated ancient human mitochondrial DNA-like reads with known genomic origins.

The experiments assess a defined set of fragment lengths, PMD conditions, aligner configurations, and reference-boundary locations.

The ten damage realizations are independent stochastic damage simulations on fixed fragments, not independent biological samples.

Consequently, the findings characterize the tested computational conditions and should not be interpreted as demonstrating universal superiority of either aligner or performance across all empirical ancient DNA datasets.

## Authors and Contributions

This research was conducted collaboratively by:

**Marco Bitar**  
Master's Student, Computer Science  
American University of Beirut  
GitHub: [MarcoBitar](https://github.com/MarcoBitar)

**Matéo Bitar**  
Research Collaborator and Co-author  
GitHub: [MateoBitar](https://github.com/MateoBitar)

Both authors contributed to the research project, including its experimental development, analysis, interpretation, and scientific reporting.

## Academic Context

Developed as part of **CMPS 355: Bioinformatics Algorithms**, American University of Beirut, Fall 2026.

The project investigates a focused computational biology question through controlled simulation, comparative alignment benchmarking, and reproducible quantitative analysis.

## Citation

If you use or build upon this research, please credit both authors:

Bitar, M., & Bitar, M. (2026). *Effects of fragment length and post-mortem damage on boundary-associated mapping accuracy in ancient human mitochondrial DNA using BWA-aln and Bowtie2* [Research project]. American University of Beirut.

## License

Licensing information for the original project code and third-party components will be specified separately.
