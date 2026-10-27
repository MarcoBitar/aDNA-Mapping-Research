import csv
from collections import defaultdict
from pathlib import Path

import pysam

N = 16569
LENGTHS = [30, 40, 50, 70]
THRESHOLDS = [0, 1, 5, 10, 20, 30]

RESULTS = defaultdict(
    lambda: {
        "correct": 0,
        "incorrect": 0,
        "unmapped_or_filtered": 0,
    }
)

for length in LENGTHS:
    truth_path = Path(
        f"ground_truth/circular_length_baseline/"
        f"reads_{length}bp_none.tsv"
    )

    with truth_path.open() as f:
        truth = {
            row["read_id"]: row
            for row in csv.DictReader(f, delimiter="\t")
        }

    for rep in range(1, 11):
        if rep == 1:
            path = Path(
                f"alignments/circular_pmd_main/"
                f"bowtie2_original_{length}bp_high.sam"
            )
        else:
            path = Path(
                f"alignments/circular_pmd_replicates/"
                f"bowtie2_original_{length}bp_high_rep{rep}.sam"
            )

        assert path.exists(), path

        with pysam.AlignmentFile(str(path), "r") as sam:
            for read in sam.fetch(until_eof=True):
                if read.is_secondary or read.is_supplementary:
                    continue

                expected = truth[read.query_name]
                group = expected["group"]

                start_correct = (
                    not read.is_unmapped
                    and read.reference_start + 1
                    == int(expected["start_1based"])
                )

                strand_correct = (
                    ("-" if read.is_reverse else "+")
                    == expected["strand"]
                )

                cigar_correct = (
                    read.cigarstring == f"{length}M"
                )

                correct = (
                    start_correct
                    and strand_correct
                    and cigar_correct
                )

                for threshold in THRESHOLDS:
                    key = (length, group, threshold)
                    counts = RESULTS[key]

                    if (
                        read.is_unmapped
                        or read.mapping_quality < threshold
                    ):
                        counts["unmapped_or_filtered"] += 1
                    elif correct:
                        counts["correct"] += 1
                    else:
                        counts["incorrect"] += 1

print("========== BOWTIE2 MAPQ FILTERING ==========")
print("High PMD, original reference, 10 replicates")

for length in LENGTHS:
    print(f"\n{length} bp")

    for group in [
        "original_boundary",
        "shifted_boundary",
        "interior",
    ]:
        print(f"  {group}")

        for threshold in THRESHOLDS:
            counts = RESULTS[(length, group, threshold)]

            correct = counts["correct"]
            incorrect = counts["incorrect"]
            total = sum(counts.values())

            assert total == 10 * 2 * (length - 1)

            recall = 100 * correct / total

            mapped = correct + incorrect

            precision = (
                100 * correct / mapped
                if mapped else None
            )

            precision_text = (
                f"{precision:6.2f}%"
                if precision is not None
                else "undefined"
            )

            print(
                f"    MAPQ >= {threshold:2d} "
                f"recall={recall:6.2f}% "
                f"precision={precision_text:>9s} "
                f"correct={correct:4d} "
                f"incorrect={incorrect:4d}"
            )

assert len(RESULTS) == 72

print("\nMAPQ FILTERING EVALUATION COMPLETED SUCCESSFULLY")
