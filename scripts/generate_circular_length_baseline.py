import csv
import random
from collections import Counter
from pathlib import Path

N = 16569
SHIFT = 8000
SEED = 1601
LENGTHS = [30, 40, 50, 70]

REFERENCE = Path("reference/rCRS_NC_012920.1.fasta")
OUTPUT_DIR = Path("simulation/circular_length_baseline")
TRUTH_DIR = Path("ground_truth/circular_length_baseline")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
TRUTH_DIR.mkdir(parents=True, exist_ok=True)

with REFERENCE.open() as f:
    next(f)
    reference = "".join(line.strip() for line in f).upper()

assert len(reference) == N
assert set(reference) <= set("ACGTN")
assert reference[3106] == "N"

complement = str.maketrans("ACGT", "TGCA")

def reverse_complement(sequence):
    return sequence.translate(complement)[::-1]

def circular_fragment(start, length):
    return "".join(
        reference[(start - 1 + offset) % N]
        for offset in range(length)
    )

def crosses_boundary(start, length, boundary):
    return any(
        ((start - 1 + offset) % N) + 1 == boundary
        for offset in range(length - 1)
    )

print("========== UNDAMAGED LENGTH BASELINE ==========")

for length in LENGTHS:

    rng = random.Random(SEED + length)

    original_starts = []
    shifted_starts = []
    interior_starts = []

    for start in range(1, N + 1):

        original = crosses_boundary(start, length, N)
        shifted = crosses_boundary(start, length, SHIFT)

        if original:
            original_starts.append(start)
        elif shifted:
            shifted_starts.append(start)
        else:
            interior_starts.append(start)

    assert len(original_starts) == length - 1
    assert len(shifted_starts) == length - 1

    # Select distinct interior positions.
    # Exclude interior fragments containing ambiguous reference bases.
    eligible_interior = [
        start
        for start in interior_starts
        if set(circular_fragment(start, length)) <= set("ACGT")
    ]

    assert len(eligible_interior) >= length - 1

    selected_interior = rng.sample(
        eligible_interior,
        length - 1
    )

    groups = {
        "original_boundary": original_starts,
        "shifted_boundary": shifted_starts,
        "interior": selected_interior
    }

    fasta_path = OUTPUT_DIR / f"reads_{length}bp_none.fasta"
    truth_path = TRUTH_DIR / f"reads_{length}bp_none.tsv"

    counts = Counter()
    counter = 0

    with fasta_path.open("w") as fasta, \
         truth_path.open("w", newline="") as truth:

        writer = csv.writer(truth, delimiter="\t")

        writer.writerow([
            "read_id",
            "group",
            "start_1based",
            "end_1based",
            "length",
            "strand",
            "replicate"
        ])

        for group, starts in groups.items():

            for start in starts:

                forward = circular_fragment(start, length)
                end = ((start + length - 2) % N) + 1

                for strand in ["+", "-"]:

                    counter += 1

                    sequence = (
                        forward if strand == "+"
                        else reverse_complement(forward)
                    )

                    read_id = f"read_{counter:06d}"

                    fasta.write(
                        f">{read_id}\n{sequence}\n"
                    )

                    writer.writerow([
                        read_id,
                        group,
                        start,
                        end,
                        length,
                        strand,
                        1
                    ])

                    counts[(group, strand)] += 1

    expected_per_group = 2 * (length - 1)
    expected_total = 6 * (length - 1)

    assert counter == expected_total

    for group in groups:
        assert counts[(group, "+")] == length - 1
        assert counts[(group, "-")] == length - 1

    print(f"\n{length} bp")
    print(f"  Original-boundary reads: {expected_per_group}")
    print(f"  Shifted-boundary reads:  {expected_per_group}")
    print(f"  Interior-control reads: {expected_per_group}")
    print(f"  Total reads:            {counter}")
    print(f"  FASTA: {fasta_path}")
    print(f"  Ground truth: {truth_path}")
    print("  Generation: PASSED")

print("\nALL LENGTH DATASETS GENERATED")
