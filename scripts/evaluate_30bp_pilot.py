import csv
from pathlib import Path

import pysam


ALIGNMENTS = {
    "BWA-aln": Path("alignments/bwa/pilot_30bp_undamaged.sam"),
    "Bowtie2": Path("alignments/bowtie2/pilot_30bp_undamaged.sam"),
}

OUTPUT = Path("results/pilot_30bp_ground_truth_evaluation.csv")


def parse_ground_truth(read_name):
    """
    Example:
    NC_012920.1:-:14671:14701:30e-1

    Returns the original reference, strand, start, and end.
    """
    parts = read_name.split(":")

    if len(parts) != 5:
        raise ValueError(f"Unexpected read identifier: {read_name}")

    reference = parts[0]
    strand = parts[1]
    start = int(parts[2])
    end = int(parts[3])

    if strand not in ("+", "-"):
        raise ValueError(f"Invalid strand: {read_name}")

    return reference, strand, start, end


def evaluate_alignment(sam_path, mapper):
    total = 0
    mapped = 0
    correct = 0
    incorrect = 0
    unmapped = 0

    with pysam.AlignmentFile(str(sam_path), "r") as sam:
        for read in sam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue

            total += 1

            reference, strand, start, end = parse_ground_truth(
                read.query_name
            )

            if read.is_unmapped:
                unmapped += 1
                continue

            mapped += 1

            predicted_strand = "-" if read.is_reverse else "+"

            is_correct = (
                read.reference_name == reference
                and predicted_strand == strand
                and read.reference_start == start
                and read.reference_end == end
            )

            if is_correct:
                correct += 1
            else:
                incorrect += 1
                print(
                    f"Incorrect alignment ({mapper}): "
                    f"{read.query_name}, "
                    f"predicted={read.reference_name}:"
                    f"{read.reference_start}-{read.reference_end}"
                    f"({predicted_strand})"
                )

    if total == 0:
        raise ValueError(f"No primary reads found in {sam_path}")

    return {
        "mapper": mapper,
        "total_reads": total,
        "mapped": mapped,
        "correct": correct,
        "incorrect": incorrect,
        "unmapped": unmapped,
        "mapping_rate": mapped / total,
        "correct_mapping_rate": correct / total,
        "precision_among_mapped": correct / mapped if mapped else 0,
    }


results = []

for mapper, path in ALIGNMENTS.items():
    result = evaluate_alignment(path, mapper)
    results.append(result)

    print(f"\n========== {mapper} ==========")
    print(f"Total reads: {result['total_reads']}")
    print(f"Mapped: {result['mapped']}")
    print(f"Correctly mapped: {result['correct']}")
    print(f"Incorrectly mapped: {result['incorrect']}")
    print(f"Unmapped: {result['unmapped']}")
    print(f"Mapping rate: {result['mapping_rate']:.2%}")
    print(f"Correct mapping rate: {result['correct_mapping_rate']:.2%}")
    print(
        "Precision among mapped reads: "
        f"{result['precision_among_mapped']:.2%}"
    )

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

with OUTPUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=results[0].keys())
    writer.writeheader()
    writer.writerows(results)

print(f"\nResults saved to: {OUTPUT}")
