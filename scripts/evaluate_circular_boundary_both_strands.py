import csv
from collections import defaultdict
from pathlib import Path
import pysam

N = 16569
SHIFT = 8000

TRUTH = Path(
    "ground_truth/reads_50bp_none_both_strands_rep1.tsv"
)

ALIGN_DIR = Path("alignments/circular_boundary_both_strands")

OUTPUT = Path(
    "results/circular_boundary_50bp_none_both_strands_rep1.csv"
)

FILES = {
    ("BWA-aln", "original"):
        "bwa_original_50bp_none_both_strands_rep1.sam",
    ("BWA-aln", "shifted"):
        "bwa_shifted_50bp_none_both_strands_rep1.sam",
    ("Bowtie2", "original"):
        "bowtie2_original_50bp_none_both_strands_rep1.sam",
    ("Bowtie2", "shifted"):
        "bowtie2_shifted_50bp_none_both_strands_rep1.sam",
}

GROUPS = [
    "original_boundary",
    "shifted_boundary",
    "interior",
]

STRANDS = ["+", "-"]

with TRUTH.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

truth = {row["read_id"]: row for row in rows}

assert len(rows) == 3000
assert len(truth) == 3000

results = []

for (aligner, ref_type), filename in FILES.items():

    counts = defaultdict(
        lambda: {"correct": 0, "incorrect": 0, "unmapped": 0}
    )

    seen = set()

    with pysam.AlignmentFile(ALIGN_DIR / filename, "r") as sam:
        for read in sam.fetch(until_eof=True):

            if read.is_secondary or read.is_supplementary:
                continue

            read_id = read.query_name

            assert read_id in truth, f"Unknown read: {read_id}"
            assert read_id not in seen, f"Duplicate read: {read_id}"
            seen.add(read_id)

            row = truth[read_id]
            group = row["group"]
            strand = row["strand"]
            key = (group, strand)

            if read.is_unmapped:
                counts[key]["unmapped"] += 1
                continue

            expected_start = int(row["start_1based"])
            length = int(row["length"])

            observed_start = read.reference_start + 1

            if ref_type == "shifted":
                observed_start = (
                    (observed_start - 1 + SHIFT) % N
                ) + 1

            observed_strand = "-" if read.is_reverse else "+"

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
        f"Missing reads in {filename}: "
        f"{len(set(truth) - seen)}"
    )

    for group in GROUPS:
        for strand in STRANDS:

            c = counts[(group, strand)]

            correct = c["correct"]
            incorrect = c["incorrect"]
            unmapped = c["unmapped"]

            total = correct + incorrect + unmapped
            assert total == 500

            mapped = correct + incorrect

            recall = correct / total

            precision = (
                correct / mapped if mapped else None
            )

            results.append({
                "aligner": aligner,
                "reference": ref_type,
                "group": group,
                "strand": strand,
                "correct": correct,
                "incorrect": incorrect,
                "unmapped": unmapped,
                "recall": recall,
                "precision": precision,
            })

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

with OUTPUT.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "aligner", "reference", "group", "strand",
            "correct", "incorrect", "unmapped",
            "recall", "precision"
        ]
    )
    writer.writeheader()
    writer.writerows(results)

print("========== STRAND-AWARE EVALUATION ==========")

for row in results:
    precision = (
        "NA"
        if row["precision"] is None
        else f"{row['precision']:.2%}"
    )

    print(
        f"{row['aligner']:8s} "
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
print("Evaluation: COMPLETED")
