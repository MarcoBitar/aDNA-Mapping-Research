import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pysam

# ============================================================
# CONFIGURATION
# ============================================================

ROOT = Path(".")
FIGURES = ROOT / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

SUMMARY_CSV = ROOT / "results/circular_pmd_replicate_summary.csv"
REPLICATE_CSV = ROOT / "results/circular_pmd_replicate_recall.csv"
ALL_CSV = ROOT / "results/circular_pmd_replicates_all.csv"
SHIFT_CSV = ROOT / "results/circular_reference_shifting_comparison.csv"

TRUTH_DIR = ROOT / "ground_truth/circular_length_baseline"
MAIN_DIR = ROOT / "alignments/circular_pmd_main"
REPLICATE_DIR = ROOT / "alignments/circular_pmd_replicates"

N = 16569
SHIFT = 8000

LENGTHS = [30, 40, 50, 70]
DAMAGE = ["low", "moderate", "high"]
ALIGNERS = ["bwa", "bowtie2"]
REFERENCES = ["original", "shifted"]
GROUPS = ["original_boundary", "shifted_boundary", "interior"]
THRESHOLDS = [0, 1, 5, 10, 20, 30]

COLORS = {
    "low": "#228B78",
    "moderate": "#D89A2B",
    "high": "#C84B55",
    "bwa": "#2878B5",
    "bowtie2": "#D97732",
    "original": "#2878B5",
    "shifted": "#D97732",
}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

# ============================================================
# UTILITIES
# ============================================================

def read_csv(path):
    assert path.exists(), f"Missing input: {path}"
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def save_figure(fig, filename):
    png = FIGURES / f"{filename}.png"
    pdf = FIGURES / f"{filename}.pdf"

    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    assert png.exists() and png.stat().st_size > 0
    assert pdf.exists() and pdf.stat().st_size > 0

    print(f"  Saved: {png}")
    print(f"  Saved: {pdf}")


def style_axis(ax, ylabel=None):
    ax.grid(axis="y", alpha=0.22, linestyle="--")
    ax.set_axisbelow(True)
    if ylabel:
        ax.set_ylabel(ylabel)


summary_rows = read_csv(SUMMARY_CSV)
replicate_rows = read_csv(REPLICATE_CSV)
all_rows = read_csv(ALL_CSV)
shift_rows = read_csv(SHIFT_CSV)

assert len(summary_rows) == 144
assert len(replicate_rows) == 1440
assert len(all_rows) == 2976
assert len(shift_rows) == 72

summary = {}
for r in summary_rows:
    key = (
        int(r["length"]),
        r["damage"],
        r["aligner"],
        r["reference"],
        r["group"],
    )
    assert key not in summary
    summary[key] = r

shift_summary = {}
for r in shift_rows:
    key = (
        int(r["length"]),
        r["damage"],
        r["aligner"],
        r["group"],
    )
    assert key not in shift_summary
    shift_summary[key] = r

print("Validated all four input CSV files.")

# ============================================================
# FIGURE 1 — FRAGMENT LENGTH × PMD SEVERITY
# ============================================================

print("\nFIGURE 1 — Fragment length and PMD severity")

fig, axes = plt.subplots(
    1, 2, figsize=(11, 4.8),
    sharey=True, constrained_layout=True
)

for ax, aligner in zip(axes, ALIGNERS):
    for damage, marker in zip(DAMAGE, ["o", "s", "^"]):
        means = []
        sds = []

        for length in LENGTHS:
            r = summary[
                (length, damage, aligner,
                 "shifted", "original_boundary")
            ]
            assert int(r["replicates"]) == 10

            means.append(100 * float(r["mean_recall"]))
            sds.append(100 * float(r["sd_recall"]))

        ax.errorbar(
            LENGTHS, means, yerr=sds,
            marker=marker, markersize=6,
            capsize=3, linewidth=1.8,
            color=COLORS[damage],
            label=damage.capitalize() + " PMD",
        )

    ax.set_title("BWA-aln" if aligner == "bwa" else "Bowtie2")
    ax.set_xlabel("Fragment length (bp)")
    ax.set_xticks(LENGTHS)
    ax.set_ylim(90, 101)
    style_axis(ax)

axes[0].set_ylabel("Strict full-length recall (%)")
axes[0].legend(frameon=False, loc="lower right")

fig.suptitle(
    "Fragment length and post-mortem damage",
    fontweight="bold"
)

save_figure(fig, "figure1_pmd_length_recall")

# ============================================================
# FIGURE 2 — ORIGINAL VS SHIFTED REFERENCE
# ============================================================

print("\nFIGURE 2 — Reference-boundary relocation")

fig, axes = plt.subplots(
    1, 2, figsize=(12, 4.8),
    sharey=True, constrained_layout=True
)

group_labels = [
    "Original\nboundary",
    "Shifted\nboundary",
    "Interior",
]

for ax, aligner in zip(axes, ALIGNERS):
    x = np.arange(3)
    width = 0.36

    for offset, reference in [
        (-width / 2, "original"),
        (width / 2, "shifted"),
    ]:
        means = []
        sds = []

        for group in GROUPS:
            r = summary[
                (30, "high", aligner, reference, group)
            ]
            means.append(100 * float(r["mean_recall"]))
            sds.append(100 * float(r["sd_recall"]))

        ax.bar(
            x + offset, means,
            width=width,
            yerr=sds,
            capsize=3,
            label=reference.capitalize() + " reference",
            color=COLORS[reference],
            alpha=0.9,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(group_labels)
    ax.set_ylim(0, 110)
    ax.set_title("BWA-aln" if aligner == "bwa" else "Bowtie2")
    style_axis(ax)

axes[0].set_ylabel("Strict full-length recall (%)")
axes[0].legend(frameon=False, loc="upper center")

fig.suptitle(
    "Reference shifting relocates the mapping boundary",
    fontweight="bold"
)

save_figure(fig, "figure2_reference_boundary")

# ============================================================
# FIGURE 3 — PAIRED REFERENCE-SHIFTING EFFECT
# ============================================================

print("\nFIGURE 3 — Paired reference-shifting effect")

fig, axes = plt.subplots(
    1, 2, figsize=(11, 4.8),
    sharey=True, constrained_layout=True
)

for ax, aligner in zip(axes, ALIGNERS):
    for group, marker, color in [
        ("original_boundary", "o", "#228B78"),
        ("shifted_boundary", "s", "#C84B55"),
        ("interior", "^", "#777777"),
    ]:
        means = []
        sds = []

        for length in LENGTHS:
            r = shift_summary[
                (length, "high", aligner, group)
            ]
            means.append(float(r["mean_delta_pp"]))
            sds.append(float(r["sd_delta_pp"]))

        ax.errorbar(
            LENGTHS, means, yerr=sds,
            marker=marker,
            linewidth=1.8,
            markersize=6,
            capsize=3,
            color=color,
            label=group.replace("_", " ").capitalize(),
        )

    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(LENGTHS)
    ax.set_xlabel("Fragment length (bp)")
    ax.set_ylim(-110, 110)
    ax.set_title("BWA-aln" if aligner == "bwa" else "Bowtie2")
    style_axis(ax)

axes[0].set_ylabel(
    "Recall change: shifted − original (pp)"
)
axes[0].legend(frameon=False, fontsize=8)

fig.suptitle(
    "Paired effect of reference shifting under high PMD",
    fontweight="bold"
)

save_figure(fig, "figure3_shifting_effect")

# ============================================================
# FIGURE 4 — UNDAMAGED BASELINE VS PMD
# ============================================================

print("\nFIGURE 4 — Undamaged baseline versus PMD")

baseline = defaultdict(lambda: {
    "correct": 0,
    "total": 0,
})

for r in all_rows:
    if r["damage"] != "none":
        continue
    if r["reference"] != "shifted":
        continue
    if r["group"] != "original_boundary":
        continue

    key = (int(r["length"]), r["aligner"])
    baseline[key]["correct"] += int(r["correct"])
    baseline[key]["total"] += (
        int(r["correct"])
        + int(r["incorrect"])
        + int(r["unmapped"])
    )

assert len(baseline) == 8

fig, axes = plt.subplots(
    1, 2, figsize=(11, 4.8),
    sharey=True, constrained_layout=True
)

for ax, aligner in zip(axes, ALIGNERS):
    x = np.arange(len(LENGTHS))
    width = 0.19

    for index, damage in enumerate(
        ["none", "low", "moderate", "high"]
    ):
        values = []

        for length in LENGTHS:
            if damage == "none":
                b = baseline[(length, aligner)]
                assert b["total"] == 2 * (length - 1)
                recall = 100 * b["correct"] / b["total"]
            else:
                r = summary[
                    (length, damage, aligner,
                     "shifted", "original_boundary")
                ]
                recall = 100 * float(r["mean_recall"])

            values.append(100 - recall)

        color = (
            "#777777"
            if damage == "none"
            else COLORS[damage]
        )

        ax.bar(
            x + (index - 1.5) * width,
            values,
            width,
            label=damage.capitalize(),
            color=color,
        )

    ax.set_xticks(x)
    ax.set_xticklabels([f"{n} bp" for n in LENGTHS])
    ax.set_title("BWA-aln" if aligner == "bwa" else "Bowtie2")
    style_axis(ax)

axes[0].set_ylabel(
    "Strict recall loss relative to 100% (pp)"
)
axes[0].legend(frameon=False)

fig.suptitle(
    "Damage-associated recall loss after boundary relocation",
    fontweight="bold"
)

save_figure(fig, "figure4_damage_baseline")

# ============================================================
# READ TRUTH AND SAM FILES FOR FIGURES 5 AND 6
# ============================================================

print("\nReading high-PMD Bowtie2 SAM alignments...")

truth_by_length = {}

for length in LENGTHS:
    path = TRUTH_DIR / f"reads_{length}bp_none.tsv"
    with path.open() as f:
        truth_by_length[length] = {
            r["read_id"]: r
            for r in csv.DictReader(f, delimiter="\t")
        }

# Store aggregate statistics only; no need to retain all reads.
boundary_errors = Counter()
mapq_errors = Counter()
mapq_filter = defaultdict(
    lambda: {
        "correct": 0,
        "incorrect": 0,
        "unmapped_or_filtered": 0,
    }
)

for length in LENGTHS:
    truth = truth_by_length[length]

    for reference in REFERENCES:
        active_group = (
            "original_boundary"
            if reference == "original"
            else "shifted_boundary"
        )

        for rep in range(1, 11):
            if rep == 1:
                path = (
                    MAIN_DIR /
                    f"bowtie2_{reference}_{length}bp_high.sam"
                )
            else:
                path = (
                    REPLICATE_DIR /
                    f"bowtie2_{reference}_{length}bp_high_rep{rep}.sam"
                )

            assert path.exists(), path

            with pysam.AlignmentFile(str(path), "r") as sam:
                for read in sam.fetch(until_eof=True):
                    if read.is_secondary or read.is_supplementary:
                        continue

                    expected = truth[read.query_name]
                    group = expected["group"]

                    start_correct = False

                    if not read.is_unmapped:
                        start = read.reference_start + 1

                        if reference == "shifted":
                            start = (
                                (start - 1 + SHIFT) % N
                            ) + 1

                        start_correct = (
                            start == int(expected["start_1based"])
                        )

                    strand_correct = (
                        ("-" if read.is_reverse else "+")
                        == expected["strand"]
                    )

                    cigar_correct = (
                        read.cigarstring == f"{length}M"
                    )

                    correct = (
                        not read.is_unmapped
                        and start_correct
                        and strand_correct
                        and cigar_correct
                    )

                    # Active-boundary error audit
                    if (
                        group == active_group
                        and not read.is_unmapped
                        and not correct
                    ):
                        boundary_errors["incorrect"] += 1
                        boundary_errors[
                            f"length_{length}"
                        ] += 1

                        if start_correct:
                            boundary_errors[
                                "correct_start"
                            ] += 1

                        if strand_correct:
                            boundary_errors[
                                "correct_strand"
                            ] += 1

                        operations = {
                            op for op, size
                            in (read.cigartuples or [])
                        }

                        if 1 in operations:
                            boundary_errors[
                                "insertion"
                            ] += 1

                        if 2 in operations:
                            boundary_errors[
                                "deletion"
                            ] += 1

                        if 4 in operations:
                            boundary_errors[
                                "softclip"
                            ] += 1

                        mapq_errors[read.mapping_quality] += 1

                    # MAPQ filtering on original reference
                    if reference == "original":
                        for threshold in THRESHOLDS:
                            key = (length, group, threshold)
                            counts = mapq_filter[key]

                            if (
                                read.is_unmapped
                                or read.mapping_quality < threshold
                            ):
                                counts[
                                    "unmapped_or_filtered"
                                ] += 1
                            elif correct:
                                counts["correct"] += 1
                            else:
                                counts["incorrect"] += 1

assert boundary_errors["incorrect"] == 721
assert boundary_errors["insertion"] == 721
assert boundary_errors["deletion"] == 0
assert boundary_errors["softclip"] == 0
assert boundary_errors["correct_start"] == 381
assert boundary_errors["correct_strand"] == 721

assert mapq_errors == Counter({
    0: 564,
    3: 108,
    8: 48,
    23: 1,
})

assert len(mapq_filter) == 72

# Validate threshold-zero counts against the main evaluator.
expected_counts = defaultdict(
    lambda: {"correct": 0, "incorrect": 0, "unmapped": 0}
)

for r in all_rows:
    if (
        r["damage"] == "high"
        and r["aligner"] == "bowtie2"
        and r["reference"] == "original"
    ):
        key = (int(r["length"]), r["group"])
        for metric in ["correct", "incorrect", "unmapped"]:
            expected_counts[key][metric] += int(r[metric])

for key, expected in expected_counts.items():
    observed = mapq_filter[(key[0], key[1], 0)]
    assert observed["correct"] == expected["correct"]
    assert observed["incorrect"] == expected["incorrect"]
    assert (
        observed["unmapped_or_filtered"]
        == expected["unmapped"]
    )

print("  All 721 boundary errors verified.")
print("  MAPQ distribution verified.")
print("  MAPQ >= 0 matches the main evaluator.")

# ============================================================
# FIGURE 5 — BOWTIE2 ALIGNMENT ARTIFACTS
# ============================================================

print("\nFIGURE 5 — Bowtie2 alignment artifacts")

fig, axes = plt.subplots(
    1, 3, figsize=(13.5, 4.2),
    constrained_layout=True
)

# Panel A: Incorrect alignments by fragment length
length_counts = [
    boundary_errors[f"length_{n}"]
    for n in LENGTHS
]

axes[0].bar(
    [str(n) for n in LENGTHS],
    length_counts,
    color="#C84B55",
)
axes[0].set_title("A. Incorrect alignments")
axes[0].set_xlabel("Fragment length (bp)")
style_axis(axes[0], "Number of alignments")

# Panel B: Correct vs incorrect genomic start
axes[1].bar(
    ["Correct start", "Incorrect start"],
    [
        boundary_errors["correct_start"],
        (
            boundary_errors["incorrect"]
            - boundary_errors["correct_start"]
        ),
    ],
    color=["#228B78", "#C84B55"],
)
axes[1].set_title("B. Genomic start")
style_axis(axes[1], "Number of alignments")

# Panel C: MAPQ distribution
mapq_values = sorted(mapq_errors)

axes[2].bar(
    [str(v) for v in mapq_values],
    [mapq_errors[v] for v in mapq_values],
    color="#2878B5",
)
axes[2].set_title("C. Mapping quality")
axes[2].set_xlabel("MAPQ")
style_axis(axes[2], "Number of alignments")

fig.suptitle(
    "Bowtie2 active-boundary alignment artifacts under high PMD",
    fontweight="bold"
)

save_figure(fig, "figure5_mapping_errors")

# ============================================================
# FIGURE 6 — MAPQ FILTERING TRADEOFF
# ============================================================

print("\nFIGURE 6 — MAPQ filtering tradeoff")

fig, axes = plt.subplots(
    1, 2, figsize=(11, 4.8),
    constrained_layout=True
)

# Panel A: Interior strict recall by fragment length
for length in LENGTHS:
    recalls = []

    for threshold in THRESHOLDS:
        counts = mapq_filter[
            (length, "interior", threshold)
        ]

        total = sum(counts.values())
        assert total == 10 * 2 * (length - 1)

        recalls.append(
            100 * counts["correct"] / total
        )

    axes[0].plot(
        THRESHOLDS,
        recalls,
        marker="o",
        linewidth=1.7,
        label=f"{length} bp",
    )

axes[0].set_title("A. Interior-read strict recall")
axes[0].set_xlabel("Minimum MAPQ threshold")
axes[0].set_xticks(THRESHOLDS)
axes[0].set_ylim(50, 102)
axes[0].legend(frameon=False)
style_axis(axes[0], "Strict recall (%)")

# Panel B: Incorrect active-boundary alignments retained
incorrect_retained = []

for threshold in THRESHOLDS:
    total = sum(
        mapq_filter[
            (length, "original_boundary", threshold)
        ]["incorrect"]
        for length in LENGTHS
    )
    incorrect_retained.append(total)

assert incorrect_retained == [
    408, 86, 31, 1, 1, 0
]

axes[1].bar(
    [str(t) for t in THRESHOLDS],
    incorrect_retained,
    color="#C84B55",
)
axes[1].set_title(
    "B. Incorrect boundary alignments retained"
)
axes[1].set_xlabel("Minimum MAPQ threshold")
style_axis(axes[1], "Incorrect alignments")

fig.suptitle(
    "MAPQ filtering: error suppression versus recall loss",
    fontweight="bold"
)

save_figure(fig, "figure6_mapq_filtering")

# ============================================================
# MANUSCRIPT FIGURE CAPTIONS
# ============================================================

captions = """
FIGURE 1 — Fragment length and post-mortem damage

Mean strict full-length mapping recall for simulated ancient
human mitochondrial DNA reads of 30, 40, 50, and 70 bp.
Original-boundary-spanning reads were aligned to the reference
shifted by 8,000 bp. Curves represent low, moderate, and high
simulated post-mortem damage for BWA-aln and Bowtie2.
Error bars indicate ±1 sample standard deviation across ten
independently seeded damage realizations of fixed fragments.


FIGURE 2 — Reference-boundary relocation

Strict full-length mapping recall for 30 bp reads under high
simulated post-mortem damage, grouped by genomic position
relative to the original and shifted linearization boundaries.
Bars show mean recall across ten damage realizations;
error bars indicate ±1 sample standard deviation.
Reference shifting restores mapping at the original boundary
while relocating the boundary-associated failure to a new
genomic position.


FIGURE 3 — Paired effect of reference shifting

Mean paired difference in strict full-length recall between
shifted and original reference alignments under high simulated
post-mortem damage. Positive values indicate improved recall
after shifting. Negative values indicate reduced recall.
Each difference compares the same simulated reads and damage
realization. Error bars indicate ±1 sample standard deviation
of paired differences across ten realizations.


FIGURE 4 — Damage-associated recall loss

Loss of strict full-length mapping recall relative to 100%
for original-boundary reads aligned to the shifted reference.
Undamaged reads constitute a single deterministic baseline.
Low, moderate, and high PMD bars represent mean recall loss
across ten stochastic damage realizations. No uncertainty
bars are shown in this figure; replicate variation is
displayed in Figure 1.


FIGURE 5 — Bowtie2 alignment artifacts

Characteristics of 721 incorrect Bowtie2 alignments involving
reads crossing the active reference boundary under high PMD.
Results combine ten damage realizations and both reference
configurations. All 721 incorrect alignments contained
insertion operations, while 381 retained the correct genomic
start. Mapping quality was predominantly low. Incorrectness
was defined using exact start, strand, and full-length CIGAR.


FIGURE 6 — MAPQ filtering tradeoff

Effect of minimum mapping-quality thresholds on Bowtie2
alignments under high PMD using the original reference.
Panel A shows strict recall for interior reads by fragment
length. Panel B shows the number of incorrect original-boundary
alignments retained after filtering. Counts are aggregated
across ten damage realizations. MAPQ filtering suppresses
incorrect alignments but can also remove correctly mapped reads.
"""

caption_path = FIGURES / "figure_captions.txt"
caption_path.write_text(captions.strip() + "\n")

# ============================================================
# FINAL CHECK
# ============================================================

for number in range(1, 7):
    png_files = list(FIGURES.glob(f"figure{number}_*.png"))
    pdf_files = list(FIGURES.glob(f"figure{number}_*.pdf"))

    assert len(png_files) == 1
    assert len(pdf_files) == 1

assert caption_path.exists()

print("\n" + "=" * 55)
print("ALL SIX MANUSCRIPT FIGURES GENERATED SUCCESSFULLY")
print("=" * 55)
print("PNG figures: 6")
print("PDF figures: 6")
print("Caption file:", caption_path)
print("Output directory:", FIGURES.resolve())
