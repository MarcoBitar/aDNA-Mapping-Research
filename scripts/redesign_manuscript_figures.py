from pathlib import Path
import runpy

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

# Run the original validated analysis first.
# This makes its checked results available for the redesigned figures.
original = Path("scripts/generate_all_manuscript_figures.py")
assert original.exists(), f"Missing original script: {original}"

data = runpy.run_path(str(original))

summary = data["summary"]
shift_summary = data["shift_summary"]
baseline = data["baseline"]
boundary_errors = data["boundary_errors"]
mapq_counts = data["mapq_errors"]
filter_counts = data["mapq_filter"]

def sm(length, damage, aligner, reference, group, metric):
    row = summary[(length, damage, aligner, reference, group)]
    return float(row[metric])

def pa(length, damage, aligner, group, metric):
    row = shift_summary[(length, damage, aligner, group)]
    return float(row[metric])

error_type = {
    "correct_start": boundary_errors["correct_start"],
    "wrong_start": (
        boundary_errors["incorrect"]
        - boundary_errors["correct_start"]
    ),
}

errors_by_length = {
    length: boundary_errors[f"length_{length}"]
    for length in [30, 40, 50, 70]
}

OUT = Path("figures")
LENGTHS = [30, 40, 50, 70]
DAMAGE = ["low", "moderate", "high"]
GROUPS = ["original_boundary", "shifted_boundary", "interior"]

BLUE = "#3977B8"
ORANGE = "#E59A43"
GREEN = "#36A184"
RED = "#D66565"
PURPLE = "#8574C5"
GRAY = "#788494"

DAMAGE_COLORS = {
    "low": GREEN,
    "moderate": ORANGE,
    "high": RED,
}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 12,
    "axes.titlesize": 15,
    "axes.labelsize": 12,
    "axes.titleweight": "bold",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": "#BCC5CE",
    "axes.labelcolor": "#344054",
    "xtick.color": "#475467",
    "ytick.color": "#475467",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "pdf.fonttype": 42,
})

def finish(fig, name):
    fig.savefig(OUT / f"{name}.png",
                dpi=300, bbox_inches="tight",
                facecolor="white")
    fig.savefig(OUT / f"{name}.pdf",
                bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    print("REDESIGNED:", name)

def clean(ax):
    ax.grid(axis="y", color="#E6EAF0", linewidth=0.9)
    ax.set_axisbelow(True)
    ax.tick_params(length=0, pad=7)

def title(fig, heading, subtitle):
    # Reserve a dedicated header area above the subplots.
    # This prevents the title, subtitle, and subplot titles
    # from overlapping.
    fig.set_layout_engine(None)

    fig.subplots_adjust(
        top=0.69,
        bottom=0.18,
        left=0.09,
        right=0.96,
        wspace=0.25,
    )

    fig.suptitle(
        heading,
        x=0.09,
        y=0.97,
        ha="left",
        va="top",
        fontsize=19,
        fontweight="bold",
    )

    fig.text(
        0.09,
        0.865,
        subtitle,
        ha="left",
        va="top",
        fontsize=11,
        color="#667085",
    )

# ============================================================
# FIGURE 1 — CLEAN LINE GRAPHS
# ============================================================

fig, axes = plt.subplots(
    1, 2, figsize=(14, 8.0),
    sharey=True, constrained_layout=True
)

for ax, aligner in zip(axes, ["bwa", "bowtie2"]):
    for damage, marker in zip(DAMAGE, ["o", "s", "^"]):
        means = [
            100 * sm(L, damage, aligner, "shifted",
                     "original_boundary", "mean_recall")
            for L in LENGTHS
        ]
        errors = [
            100 * sm(L, damage, aligner, "shifted",
                     "original_boundary", "sd_recall")
            for L in LENGTHS
        ]

        ax.errorbar(
            LENGTHS, means, yerr=errors,
            color=DAMAGE_COLORS[damage],
            marker=marker, markersize=8,
            linewidth=2.5, capsize=4,
            label=damage.capitalize() + " PMD"
        )

    ax.set_title(
        "BWA-aln" if aligner == "bwa" else "Bowtie2",
        pad=24
    )
    ax.set_xlabel("Fragment length (bp)")
    ax.set_xticks(LENGTHS)
    ax.set_ylim(89, 102)
    clean(ax)

axes[0].set_ylabel("Strict mapping recall (%)")

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(
    handles, labels, loc="outside lower center",
    ncol=3, frameon=False, fontsize=12
)

title(
    fig,
    "01  Fragment length and DNA damage",
    "Original-boundary reads mapped to shifted reference | Mean ± SD, 10 damage realizations"
)

# Dedicated Figure 1 layout with enough room for everything.
fig.set_layout_engine(None)

fig.set_size_inches(14, 8)

fig.subplots_adjust(
    left=0.09,
    right=0.96,
    top=0.75,
    bottom=0.19,
    wspace=0.22,
)

# Main title and subtitle, separated from subplot headings.
fig.suptitle(
    "01  Fragment length and DNA damage",
    x=0.09,
    y=0.97,
    ha="left",
    va="top",
    fontsize=20,
    fontweight="bold",
)

# Remove the previous subtitle and replace it.
for text in list(fig.texts):
    if "Original-boundary reads mapped" in text.get_text():
        text.remove()

fig.text(
    0.09,
    0.90,
    "Original-boundary reads mapped to shifted reference | Mean ± SD, 10 damage realizations",
    ha="left",
    va="top",
    fontsize=12,
    color="#667085",
)

# Subplot headings remain directly above their graphs.
for ax in axes:
    ax.set_title(
        ax.get_title(),
        fontsize=15,
        fontweight="bold",
        pad=12,
    )

# Put the legend below the graphs.
for legend in list(fig.legends):
    legend.remove()

handles, labels = axes[0].get_legend_handles_labels()

fig.legend(
    handles,
    labels,
    loc="lower center",
    bbox_to_anchor=(0.5, 0.055),
    ncol=3,
    frameon=False,
    fontsize=12,
)

finish(fig, "figure1_pmd_length_recall")

# ============================================================
# FIGURE 2 — LARGE GROUPED BAR CHARTS
# ============================================================

fig, axes = plt.subplots(
    1, 2, figsize=(14, 5.8),
    sharey=True, constrained_layout=True
)

labels = ["Original\nboundary", "Shifted\nboundary", "Interior"]
x = np.arange(3)
width = 0.35

for ax, aligner in zip(axes, ["bwa", "bowtie2"]):
    for j, reference in enumerate(["original", "shifted"]):
        means = [
            100 * sm(30, "high", aligner, reference,
                     group, "mean_recall")
            for group in GROUPS
        ]
        errors = [
            100 * sm(30, "high", aligner, reference,
                     group, "sd_recall")
            for group in GROUPS
        ]

        bars = ax.bar(
            x + (j - 0.5) * width,
            means, width,
            yerr=errors, capsize=4,
            color=BLUE if j == 0 else ORANGE,
            label=reference.capitalize() + " reference",
            edgecolor="white"
        )

        for bar, value in zip(bars, means):
            if value < 0.05:
                ax.text(
                    bar.get_x() + bar.get_width()/2,
                    3, "0%",
                    ha="center", fontsize=10,
                    color="#475467"
                )

    ax.set_title(
        "BWA-aln" if aligner == "bwa" else "Bowtie2",
        pad=24
    )
    ax.set_xticks(x, labels)
    ax.set_ylim(0, 115)
    clean(ax)

axes[0].set_ylabel("Strict mapping recall (%)")

handles, labels_legend = axes[0].get_legend_handles_labels()
fig.legend(
    handles, labels_legend,
    loc="outside lower center",
    ncol=2, frameon=False
)

title(
    fig,
    "02  Shifting the reference moves the boundary problem",
    "30 bp reads | High PMD | Mean ± SD across 10 damage realizations"
)
finish(fig, "figure2_reference_boundary")

# ============================================================
# FIGURE 3 — CLEAR POSITIVE/NEGATIVE LINE GRAPHS
# ============================================================

fig, axes = plt.subplots(
    1, 2, figsize=(14, 5.8),
    sharey=True, constrained_layout=True
)

group_specs = [
    ("original_boundary", GREEN, "o", "Original boundary"),
    ("shifted_boundary", RED, "s", "Shifted boundary"),
    ("interior", GRAY, "^", "Interior"),
]

for ax, aligner in zip(axes, ["bwa", "bowtie2"]):
    for group, color, marker, label in group_specs:
        means = [
            pa(L, "high", aligner, group, "mean_delta_pp")
            for L in LENGTHS
        ]
        errors = [
            pa(L, "high", aligner, group, "sd_delta_pp")
            for L in LENGTHS
        ]

        ax.errorbar(
            LENGTHS, means, yerr=errors,
            marker=marker, markersize=8,
            linewidth=2.5, capsize=4,
            color=color, label=label
        )

    ax.axhline(0, color="#475467", linewidth=1.3)
    ax.set_title(
        "BWA-aln" if aligner == "bwa" else "Bowtie2",
        pad=24
    )
    ax.set_xlabel("Fragment length (bp)")
    ax.set_xticks(LENGTHS)
    ax.set_ylim(-110, 110)
    clean(ax)

axes[0].set_ylabel("Change in recall (percentage points)")

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(
    handles, labels,
    loc="outside lower center",
    ncol=3, frameon=False
)

title(
    fig,
    "03  Where reference shifting helps — and hurts",
    "Shifted minus original reference | High PMD | Paired mean ± SD"
)
finish(fig, "figure3_shifting_effect")

# ============================================================
# FIGURE 4 — RECALL LOSS AS GROUPED BARS
# ============================================================

fig, axes = plt.subplots(
    1, 2, figsize=(14, 5.8),
    sharey=True, constrained_layout=True
)

x = np.arange(4)
width = 0.24

for ax, aligner in zip(axes, ["bwa", "bowtie2"]):
    for j, damage in enumerate(DAMAGE):
        losses = [
            100 * (
                1 - sm(
                    L, damage, aligner,
                    "shifted", "original_boundary",
                    "mean_recall"
                )
            )
            for L in LENGTHS
        ]

        bars = ax.bar(
            x + (j - 1) * width,
            losses, width,
            color=DAMAGE_COLORS[damage],
            label=damage.capitalize() + " PMD",
            edgecolor="white"
        )

        for bar, value in zip(bars, losses):
            if value >= 0.5:
                ax.text(
                    bar.get_x() + bar.get_width()/2,
                    value + 0.12,
                    f"{value:.1f}",
                    ha="center", fontsize=9
                )

    ax.set_title(
        "BWA-aln" if aligner == "bwa" else "Bowtie2",
        pad=24
    )
    ax.set_xticks(x, [f"{L} bp" for L in LENGTHS])
    ax.set_ylim(0, 7)
    clean(ax)

axes[0].set_ylabel("Recall loss relative to 100% (pp)")

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(
    handles, labels,
    loc="outside lower center",
    ncol=3, frameon=False
)

title(
    fig,
    "04  How much recall is lost as DNA damage increases?",
    "Shifted reference | Undamaged baseline: 100% | PMD means across 10 realizations"
)
finish(fig, "figure4_damage_baseline")

# ============================================================
# FIGURE 5 — READABLE ERROR-ANALYSIS BAR CHARTS
# ============================================================

fig, axes = plt.subplots(
    1, 3, figsize=(16, 5.5),
    constrained_layout=True
)

def labeled_bars(ax, labels, values, colors, heading):
    bars = ax.bar(
        labels, values, color=colors,
        width=0.62, edgecolor="white"
    )
    for bar, value in zip(bars, values):
        ax.annotate(
            f"{value:,}",
            (bar.get_x() + bar.get_width()/2,
             bar.get_height()),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center", fontsize=11,
            fontweight="bold"
        )
    ax.set_ylim(0, max(1, max(values)) * 1.18)
    ax.set_title(heading, pad=24)
    clean(ax)

labeled_bars(
    axes[0],
    [str(L) for L in LENGTHS],
    [errors_by_length[L] for L in LENGTHS],
    [BLUE] * 4,
    "Errors by read length"
)
axes[0].set_xlabel("Fragment length (bp)")
axes[0].set_ylabel("Incorrect alignments")

labeled_bars(
    axes[1],
    ["Correct\nstart", "Wrong\nstart"],
    [
        error_type["correct_start"],
        error_type["wrong_start"]
    ],
    [GREEN, RED],
    "Reported genomic start"
)

qs = sorted(mapq_counts)
labeled_bars(
    axes[2],
    [str(q) for q in qs],
    [mapq_counts[q] for q in qs],
    [PURPLE] * len(qs),
    "Mapping-quality scores"
)
axes[2].set_xlabel("MAPQ")

title(
    fig,
    "05  What goes wrong in Bowtie2 boundary alignments?",
    "721 incorrect alignments | All contain insertion CIGAR operations | High PMD"
)
finish(fig, "figure5_mapping_errors")

# ============================================================
# FIGURE 6 — MAPQ TRADEOFF
# ============================================================

fig, axes = plt.subplots(
    1, 2, figsize=(14, 5.8),
    constrained_layout=True
)

THRESHOLDS = [0, 1, 5, 10, 20, 30]
LENGTH_COLORS = [BLUE, GREEN, ORANGE, PURPLE]

for L, color in zip(LENGTHS, LENGTH_COLORS):
    recalls = []
    for t in THRESHOLDS:
        counts = filter_counts[L, "interior", t]
        total = 10 * 2 * (L - 1)
        recalls.append(100 * counts["correct"] / total)

    axes[0].plot(
        THRESHOLDS, recalls,
        marker="o", markersize=7,
        linewidth=2.4, color=color,
        label=f"{L} bp"
    )

axes[0].set_title("Correct interior reads retained", pad=24)
axes[0].set_xlabel("Minimum MAPQ threshold")
axes[0].set_ylabel("Strict mapping recall (%)")
axes[0].set_xticks(THRESHOLDS)
axes[0].set_ylim(55, 103)
axes[0].legend(frameon=False, ncol=2)
clean(axes[0])

incorrect = [
    sum(
        filter_counts[L, "original_boundary", t]["incorrect"]
        for L in LENGTHS
    )
    for t in THRESHOLDS
]

assert incorrect == [408, 86, 31, 1, 1, 0]

labeled_bars(
    axes[1],
    [str(t) for t in THRESHOLDS],
    incorrect,
    [RED] * len(THRESHOLDS),
    "Incorrect boundary mappings retained"
)
axes[1].set_xlabel("Minimum MAPQ threshold")
axes[1].set_ylabel("Incorrect alignments")

title(
    fig,
    "06  MAPQ filtering removes errors, but also correct reads",
    "Bowtie2 | Original reference | High PMD | 10 damage realizations"
)
finish(fig, "figure6_mapq_filtering")

print("\n" + "=" * 58)
print("ALL SIX FIGURES REDESIGNED SUCCESSFULLY")
print("=" * 58)
print("Updated 6 PNG files and 6 PDF files in figures/")
