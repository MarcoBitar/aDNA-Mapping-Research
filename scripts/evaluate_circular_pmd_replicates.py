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
PMD_MAIN_DIR = Path("alignments/circular_pmd_main")
REPLICATE_DIR = Path("alignments/circular_pmd_replicates")

ORIGINAL_RESULTS = Path("results/circular_pmd_all.csv")
OUTPUT = Path("results/circular_pmd_replicates_all.csv")

FIELDNAMES = [
    "replicate",
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
    "mean_incorrect_mapq",
]


def get_alignment_path(length, damage, aligner, reference, replicate):
    if damage == "none":
        assert replicate == 0
        return (
            BASELINE_DIR
            / f"{aligner}_{reference}_{length}bp_none.sam"
        )

    if replicate == 1:
        return (
            PMD_MAIN_DIR
            / f"{aligner}_{reference}_{length}bp_{damage}.sam"
        )

    return (
        REPLICATE_DIR
        / f"{aligner}_{reference}_{length}bp_{damage}_rep{replicate}.sam"
    )


def evaluate_alignment(sam_path, truth, length, reference):
    counts = defaultdict(
        lambda: {
            "correct": 0,
            "incorrect": 0,
            "unmapped": 0,
            "correct_mapq": [],
            "incorrect_mapq": [],
        }
    )

    seen = set()

    with pysam.AlignmentFile(sam_path, "r") as sam:
        for read in sam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue

            read_id = read.query_name

            assert read_id in truth, (
                f"Unknown read: {read_id} in {sam_path}"
            )

            assert read_id not in seen, (
                f"Duplicate primary record: {read_id} in {sam_path}"
            )

            seen.add(read_id)

            truth_row = truth[read_id]

            group = truth_row["group"]
            strand = truth_row["strand"]
            key = (group, strand)

            if read.is_unmapped:
                counts[key]["unmapped"] += 1
                continue

            expected_start = int(truth_row["start_1based"])
            observed_start = read.reference_start + 1

            if reference == "shifted":
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
        f"Missing primary records in {sam_path}"
    )

    output_rows = []

    for group in GROUPS:
        for strand in STRANDS:
            c = counts[(group, strand)]

            correct = c["correct"]
            incorrect = c["incorrect"]
            unmapped = c["unmapped"]

            total = correct + incorrect + unmapped

            assert total == length - 1, (
                f"Unexpected group size in {sam_path}: "
                f"{group}, {strand}, total={total}"
            )

            mapped = correct + incorrect

            recall = correct / total

            precision = (
                correct / mapped if mapped else None
            )

            mean_correct_mapq = (
                sum(c["correct_mapq"]) / len(c["correct_mapq"])
                if c["correct_mapq"] else None
            )

            mean_incorrect_mapq = (
                sum(c["incorrect_mapq"]) / len(c["incorrect_mapq"])
                if c["incorrect_mapq"] else None
            )

            output_rows.append({
                "group": group,
                "strand": strand,
                "correct": correct,
                "incorrect": incorrect,
                "unmapped": unmapped,
                "recall": recall,
                "precision": precision,
                "mean_correct_mapq": mean_correct_mapq,
                "mean_incorrect_mapq": mean_incorrect_mapq,
            })

    return output_rows


results = []
configurations = 0

for length in LENGTHS:
    truth_path = TRUTH_DIR / f"reads_{length}bp_none.tsv"

    with truth_path.open() as f:
        truth_rows = list(csv.DictReader(f, delimiter="\t"))

    truth = {row["read_id"]: row for row in truth_rows}

    expected_total = 6 * (length - 1)

    assert len(truth_rows) == expected_total
    assert len(truth) == expected_total

    for damage in DAMAGE_LEVELS:
        replicates = [0] if damage == "none" else range(1, 11)

        for replicate in replicates:
            for aligner in ALIGNERS:
                for reference in REFERENCES:
                    sam_path = get_alignment_path(
                        length, damage, aligner,
                        reference, replicate
                    )

                    assert sam_path.is_file(), (
                        f"Missing SAM file: {sam_path}"
                    )

                    evaluated_rows = evaluate_alignment(
                        sam_path, truth, length, reference
                    )

                    for row in evaluated_rows:
                        results.append({
                            "replicate": replicate,
                            "length": length,
                            "damage": damage,
                            "aligner": aligner,
                            "reference": reference,
                            **row,
                        })

                    configurations += 1

assert configurations == 496
assert len(results) == 2976

total_correct = sum(r["correct"] for r in results)
total_incorrect = sum(r["incorrect"] for r in results)
total_unmapped = sum(r["unmapped"] for r in results)

assert (
    total_correct + total_incorrect + total_unmapped
    == 138384
)

# Independently check the original 384 evaluation rows.
with ORIGINAL_RESULTS.open() as f:
    original_rows = list(csv.DictReader(f))

old_lookup = {
    (
        int(r["length"]),
        r["damage"],
        r["aligner"],
        r["reference"],
        r["group"],
        r["strand"],
    ): r
    for r in original_rows
}

assert len(original_rows) == 384
assert len(old_lookup) == 384

checked = 0

for row in results:
    if row["replicate"] not in (0, 1):
        continue

    key = (
        row["length"],
        row["damage"],
        row["aligner"],
        row["reference"],
        row["group"],
        row["strand"],
    )

    old = old_lookup[key]

    for metric in ["correct", "incorrect", "unmapped"]:
        assert row[metric] == int(old[metric]), (
            f"Original result mismatch: {key}, {metric}"
        )

    for metric in [
        "recall",
        "precision",
        "mean_correct_mapq",
        "mean_incorrect_mapq",
    ]:
        previous = old[metric]
        current = row[metric]

        if previous == "":
            assert current is None, (
                f"Original result mismatch: {key}, {metric}"
            )
        else:
            assert current is not None
            assert abs(float(previous) - current) < 1e-12, (
                f"Original result mismatch: {key}, {metric}"
            )

    checked += 1

assert checked == 384

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

with OUTPUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
    writer.writeheader()
    writer.writerows(results)

print("========== REPLICATED PMD EVALUATION ==========")
print("Alignment configurations:", configurations)
print("Evaluation rows:", len(results))
print("Correct:", total_correct)
print("Incorrect:", total_incorrect)
print("Unmapped:", total_unmapped)
print(
    "Total:",
    total_correct + total_incorrect + total_unmapped
)
print("Original evaluation rows verified:", checked)
print("Results saved:", OUTPUT)
print("REPLICATED EVALUATION COMPLETED SUCCESSFULLY")
