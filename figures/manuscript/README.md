# Final Manuscript Figures

This directory contains the four figures included in the final CMPS355 research manuscript.

## Figure Correspondence

| Manuscript figure | Repository file | Description |
|---|---|---|
| Figure 1 | `figure1_boundary_groups_30bp.png` | Strict recall for original-boundary, shifted-boundary, and interior groups using 30-bp fragments |
| Figure 2 | `figure2_length_recall_high_pmd.png` | Strict recall by fragment length for BWA-aln and Bowtie2 under high PMD |
| Figure 3 | `figure3_bowtie2_boundary_errors.png` | Incorrect Bowtie2 active-boundary alignments by fragment length for original and shifted references |
| Figure 4 | `figure4_bowtie2_mapq_filtering.png` | Bowtie2 MAPQ filtering: retained boundary errors and interior strict recall |

## Additional Analysis Figures

The parent `figures/` directory contains six additional analysis figures in PNG and PDF formats. Their numbering follows the analysis-generation workflow and is separate from the four figures in the final manuscript.

## Interpretation

Strict correctness requires the expected original-rCRS start position, strand, and full-length match CIGAR. Therefore, strict recall should not be interpreted as coordinate-only mapping accuracy.

The reported PMD seed repetitions represent stochastic damage realizations on fixed source fragments, not independent biological replicates.

## Figure Generation and Reproducibility

The four images in this directory are the figures included in the final manuscript.

The repository also provides two figure-generation scripts:

- `scripts/generate_all_manuscript_figures.py`
- `scripts/redesign_manuscript_figures.py`

These scripts generate six analysis figures in the parent `figures/` directory. They do not currently reproduce the four final manuscript figures with identical numbering and presentation.

The manuscript figures are provided separately to establish an unambiguous correspondence between the published figures and the repository files.

The underlying experimental results are stored in the repository's `results/` directory.
