import csv
from collections import defaultdict
from pathlib import Path

import pysam

N = 16569
SHIFT = 8000

LENGTHS = [30, 40, 50, 70]
DAMAGE_LEVELS = ["none", "low", "moderate", "high"]
ALIGNERS = ["bwa", "bowtie2"]
REFERENCES = ["original", "shifted"]
GROUPS = ["original_boundary", "shifted_boundary", "interior"]
STRANDS = ["+", "-"]

TRUTH_DIR = Path("ground_truth/circular_length_baseline")
BASELINE_DIR = Path("alignments/circular_length_baseline")
PMD_DIR = Path("alignments/circular_pmd_main")

OUTPUT = Path("results/circular_pmd_all.csv")
OUTPUT.parent.mkdir(parents=True, exist_ok=True)

results = []

for length in LENGTHS:

    truth_path = TRUTH_DIR / f"reads_{length}bp_none.tsv"

    with truth_path.open() as f:
        truth_rows = list(csv.DictReader(f, delimiter="\t"))

    truth = {row["read_id"]: row for row in truth_rows}

    expected_total = 6 * (length - 1)

    assert len(truth_rows) == expected_total
    assert len(truth) == expected_total

    for damage in DAMAGE_LEVELS:
        for aligner in ALIGNERS:
            for ref_type in REFERENCES:

                align_dir = (
                    BASELINE_DIR if damage == "none"
                    else PMD_DIR
                )

                filename = (
                    f"{aligner}_{ref_type}_{length}bp_{damage}.sam"
                )

                counts = defaultdict(
                    lambda: {
                        "correct": 0,
                        "incorrect": 0,
                        "unmapped": 0,
                        "correct_mapq": [],
                        "incorrect_mapq": []
                    }
                )

                seen = set()

                with pysam.AlignmentFile(
                    align_dir / filename, "r"
                ) as sam:

                    for read in sam.fetch(until_eof=True):

                        if read.is_secondary or read.is_supplementary:
                            continue

                        read_id = read.query_name

                        assert read_id in truth, (
                            f"Unknown read: {read_id} in {filename}"
                        )

                        assert read_id not in seen, (
                            f"Duplicate primary record: "
                            f"{read_id} in {filename}"
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
                            counts[key]["correct_mapq"].append(
                                read.mapping_quality
                            )
                        else:
                            counts[key]["incorrect"] += 1
                            counts[key]["incorrect_mapq"].append(
                                read.mapping_quality
                            )

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

                        assert total == length - 1, (
                            f"Unexpected group size: "
                            f"{filename}, {group}, {strand}"
                        )

                        mapped = correct + incorrect

                        recall = correct / total

                        precision = (
                            correct / mapped if mapped else None
                        )

                        mean_correct_mapq = (
                            sum(c["correct_mapq"])
                            / len(c["correct_mapq"])
                            if c["correct_mapq"] else None
                        )

                        mean_incorrect_mapq = (
                            sum(c["incorrect_mapq"])
                            / len(c["incorrect_mapq"])
                            if c["incorrect_mapq"] else None
                        )

                        results.append({
                            "length": length,
                            "damage": damage,
                            "aligner": aligner,
                            "reference": ref_type,
                            "group": group,
                            "strand": strand,
                            "correct": correct,
                            "incorrect": incorrect,
                            "unmapped": unmapped,
                            "recall": recall,
                            "precision": precision,
                            "mean_correct_mapq": mean_correct_mapq,
                            "mean_incorrect_mapq": mean_incorrect_mapq
                        })

fieldnames = [
    "length",
    "damage",
    "aligner",
    "reference",
    "group",
    "strand",
    "correct",
    "incorrect",
    "unmapped",
    "recall",
    "precision",
    "mean_correct_mapq",
    "mean_incorrect_mapq"
]

with OUTPUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)

assert len(results) == 384

print("========== CONSOLIDATED EVALUATION ==========")
print("Evaluated configurations:", 64)
print("Evaluation rows:", len(results))
print("Results saved:", OUTPUT)

print("\n========== GLOBAL VALIDATION ==========")

total_correct = sum(row["correct"] for row in results)
total_incorrect = sum(row["incorrect"] for row in results)
total_unmapped = sum(row["unmapped"] for row in results)

print("Correct:", total_correct)
print("Incorrect:", total_incorrect)
print("Unmapped:", total_unmapped)
print(
    "Total:",
    total_correct + total_incorrect + total_unmapped
)

assert (
    total_correct + total_incorrect + total_unmapped
    == 17856
)

print("\n========== BASELINE CONSISTENCY CHECK ==========")

with open("results/circular_length_baseline.csv") as f:
    baseline_rows = list(csv.DictReader(f))

new_baseline = {
    (
        str(row["length"]),
        row["aligner"],
        row["reference"],
        row["group"],
        row["strand"]
    ): row
    for row in results
    if row["damage"] == "none"
}

assert len(baseline_rows) == 96
assert len(new_baseline) == 96

for old in baseline_rows:

    key = (
        old["length"],
        old["aligner"],
        old["reference"],
        old["group"],
        old["strand"]
    )

    new = new_baseline[key]

    for metric in ["correct", "incorrect", "unmapped"]:
        assert int(old[metric]) == new[metric], (
            f"Baseline mismatch: {key}, {metric}"
        )

    assert abs(float(old["recall"]) - new["recall"]) < 1e-12

    if old["precision"] == "":
        assert new["precision"] is None
    else:
        assert (
            abs(float(old["precision"]) - new["precision"])
            < 1e-12
        )

print("All 96 baseline rows match the original evaluator.")
print("CONSOLIDATED EVALUATION COMPLETED SUCCESSFULLY")
