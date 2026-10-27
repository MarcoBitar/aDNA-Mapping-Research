import csv
from collections import defaultdict
from pathlib import Path

import pysam

N = 16569
SHIFT = 8000
LENGTHS = [30, 40, 50, 70]

ALIGN_DIR = Path("alignments/circular_length_baseline")
TRUTH_DIR = Path("ground_truth/circular_length_baseline")
OUTPUT = Path("results/circular_length_baseline.csv")

GROUPS = [
    "original_boundary",
    "shifted_boundary",
    "interior"
]

STRANDS = ["+", "-"]

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

results = []

for length in LENGTHS:

    truth_path = TRUTH_DIR / f"reads_{length}bp_none.tsv"

    with truth_path.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    truth = {row["read_id"]: row for row in rows}

    expected_total = 6 * (length - 1)

    assert len(rows) == expected_total
    assert len(truth) == expected_total

    for aligner in ["bwa", "bowtie2"]:
        for ref_type in ["original", "shifted"]:

            filename = (
                f"{aligner}_{ref_type}_{length}bp_none.sam"
            )

            counts = defaultdict(
                lambda: {
                    "correct": 0,
                    "incorrect": 0,
                    "unmapped": 0
                }
            )

            seen = set()

            with pysam.AlignmentFile(
                ALIGN_DIR / filename, "r"
            ) as sam:

                for read in sam.fetch(until_eof=True):

                    if read.is_secondary or read.is_supplementary:
                        continue

                    read_id = read.query_name

                    assert read_id in truth, (
                        f"Unknown read: {read_id}"
                    )

                    assert read_id not in seen, (
                        f"Duplicate primary record: {read_id}"
                    )

                    seen.add(read_id)

                    row = truth[read_id]

                    group = row["group"]
                    strand = row["strand"]
                    key = (group, strand)

                    if read.is_unmapped:
                        counts[key]["unmapped"] += 1
                        continue

                    expected_start = int(row["start_1based"])

                    observed_start = read.reference_start + 1

                    if ref_type == "shifted":
                        observed_start = (
                            (observed_start - 1 + SHIFT) % N
                        ) + 1

                    observed_strand = (
                        "-" if read.is_reverse else "+"
                    )

                    correct = (
                        observed_start == expected_start
                        and observed_strand == strand
                        and read.cigarstring == f"{length}M"
                    )

                    if correct:
                        counts[key]["correct"] += 1
                    else:
                        counts[key]["incorrect"] += 1

            assert seen == set(truth), (
                f"Missing primary records in {filename}"
            )

            for group in GROUPS:
                for strand in STRANDS:

                    c = counts[(group, strand)]

                    correct = c["correct"]
                    incorrect = c["incorrect"]
                    unmapped = c["unmapped"]

                    total = correct + incorrect + unmapped

                    assert total == length - 1

                    mapped = correct + incorrect

                    recall = correct / total

                    precision = (
                        correct / mapped if mapped else None
                    )

                    results.append({
                        "length": length,
                        "aligner": aligner,
                        "reference": ref_type,
                        "group": group,
                        "strand": strand,
                        "correct": correct,
                        "incorrect": incorrect,
                        "unmapped": unmapped,
                        "recall": recall,
                        "precision": precision
                    })

with OUTPUT.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "length",
            "aligner",
            "reference",
            "group",
            "strand",
            "correct",
            "incorrect",
            "unmapped",
            "recall",
            "precision"
        ]
    )
    writer.writeheader()
    writer.writerows(results)

print("========== LENGTH BASELINE EVALUATION ==========")

for length in LENGTHS:

    print(f"\n========== {length} BP ==========")

    for row in results:
        if row["length"] != length:
            continue

        precision = (
            "NA"
            if row["precision"] is None
            else f"{row['precision']:.2%}"
        )

        print(
            f"{row['aligner']:7s} "
            f"{row['reference']:8s} "
            f"{row['group']:19s} "
            f"strand={row['strand']} "
            f"correct={row['correct']:3d} "
            f"incorrect={row['incorrect']:3d} "
            f"unmapped={row['unmapped']:3d} "
            f"recall={row['recall']:.2%} "
            f"precision={precision}"
        )

print("\nResults saved to:", OUTPUT)
print("Total evaluation rows:", len(results))
print("EVALUATION: COMPLETED")
