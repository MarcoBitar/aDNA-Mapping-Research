import csv
from collections import defaultdict
from pathlib import Path

import pysam

N = 16569
SHIFT = 8000

TRUTH_PATH = Path(
    "ground_truth/reads_50bp_none_rep1.tsv"
)

ALIGNMENT_DIR = Path(
    "alignments/circular_boundary_main"
)

OUTPUT_PATH = Path(
    "results/circular_boundary_50bp_none_rep1.csv"
)

FILES = {
    ("BWA-aln", "original"):
        "bwa_original_50bp_none_rep1.sam",

    ("BWA-aln", "shifted"):
        "bwa_shifted_50bp_none_rep1.sam",

    ("Bowtie2", "original"):
        "bowtie2_original_50bp_none_rep1.sam",

    ("Bowtie2", "shifted"):
        "bowtie2_shifted_50bp_none_rep1.sam",
}

GROUPS = [
    "original_boundary",
    "shifted_boundary",
    "interior",
]

def convert_position(position, reference_type):
    if reference_type == "shifted":
        return ((position - 1 + SHIFT) % N) + 1
    return position

with TRUTH_PATH.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

truth = {row["read_id"]: row for row in rows}

assert len(rows) == 3000
assert len(truth) == 3000

summary = []

for (aligner, reference_type), filename in FILES.items():

    path = ALIGNMENT_DIR / filename

    results = defaultdict(
        lambda: {
            "correct": 0,
            "incorrect": 0,
            "unmapped": 0,
        }
    )

    seen = set()

    with pysam.AlignmentFile(str(path), "r") as sam:

        for read in sam.fetch(until_eof=True):

            if read.is_secondary or read.is_supplementary:
                continue

            read_id = read.query_name

            assert read_id in truth
            assert read_id not in seen
            seen.add(read_id)

            row = truth[read_id]
            group = row["group"]

            expected_start = int(row["start_1based"])
            expected_length = int(row["length"])

            if read.is_unmapped:
                results[group]["unmapped"] += 1
                continue

            mapped_start = convert_position(
                read.reference_start + 1,
                reference_type
            )

            correct = (
                mapped_start == expected_start
                and not read.is_reverse
                and read.cigarstring == f"{expected_length}M"
            )

            if correct:
                results[group]["correct"] += 1
            else:
                results[group]["incorrect"] += 1

    assert seen == set(truth)

    for group in GROUPS:
        counts = results[group]

        correct = counts["correct"]
        incorrect = counts["incorrect"]
        unmapped = counts["unmapped"]

        total = correct + incorrect + unmapped
        mapped = correct + incorrect

        assert total == 1000

        recall = correct / total
        precision = (
            correct / mapped
            if mapped > 0
            else None
        )

        summary.append({
            "aligner": aligner,
            "reference": reference_type,
            "group": group,
            "total": total,
            "correct": correct,
            "incorrect": incorrect,
            "unmapped": unmapped,
            "recall": round(recall, 6),
            "precision": (
                round(precision, 6)
                if precision is not None
                else "NA"
            ),
        })

with OUTPUT_PATH.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "aligner",
            "reference",
            "group",
            "total",
            "correct",
            "incorrect",
            "unmapped",
            "recall",
            "precision",
        ]
    )

    writer.writeheader()
    writer.writerows(summary)

print("========== BOUNDARY-SPECIFIC EVALUATION ==========")

for row in summary:
    print(
        f"{row['aligner']:8s} "
        f"{row['reference']:8s} "
        f"{row['group']:20s} "
        f"correct={row['correct']:4d} "
        f"incorrect={row['incorrect']:4d} "
        f"unmapped={row['unmapped']:4d} "
        f"recall={row['recall']:.2%} "
        f"precision={row['precision']}"
    )

print("\nResults saved to:", OUTPUT_PATH)
print("Evaluation: COMPLETED")
