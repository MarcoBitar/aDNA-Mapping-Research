import csv
import random
from pathlib import Path

SEED = 1601
LENGTH = 50
PER_STRAND = 500
SHIFT = 8000
REPLICATE = 1

rng = random.Random(SEED)

reference_path = Path(
    "reference/rCRS_NC_012920.1.fasta"
)

fasta_path = Path(
    "simulation/circular_boundary_main/"
    "reads_50bp_none_both_strands_rep1.fasta"
)

truth_path = Path(
    "ground_truth/"
    "reads_50bp_none_both_strands_rep1.tsv"
)

with reference_path.open() as f:
    next(f)
    reference = "".join(
        line.strip() for line in f
    ).upper()

N = len(reference)
assert N == 16569
assert 1 < LENGTH < SHIFT < N - LENGTH

complement = str.maketrans("ACGT", "TGCA")

def reverse_complement(sequence):
    return sequence.translate(complement)[::-1]

def circular_fragment(start):
    return "".join(
        reference[(start - 1 + i) % N]
        for i in range(LENGTH)
    )

def crosses_boundary(start, boundary):
    return any(
        ((start - 1 + i) % N) + 1 == boundary
        for i in range(LENGTH - 1)
    )

groups = {
    "original_boundary": [],
    "shifted_boundary": [],
    "interior": [],
}

for start in range(1, N + 1):
    original = crosses_boundary(start, N)
    shifted = crosses_boundary(start, SHIFT)

    if original:
        groups["original_boundary"].append(start)
    elif shifted:
        groups["shifted_boundary"].append(start)
    else:
        groups["interior"].append(start)

assert len(groups["original_boundary"]) == LENGTH - 1
assert len(groups["shifted_boundary"]) == LENGTH - 1
assert len(groups["interior"]) == N - 2 * (LENGTH - 1)

fasta_path.parent.mkdir(parents=True, exist_ok=True)
truth_path.parent.mkdir(parents=True, exist_ok=True)

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

    for group, eligible_starts in groups.items():
        for strand in ["+", "-"]:
            for _ in range(PER_STRAND):

                counter += 1
                start = rng.choice(eligible_starts)

                forward_sequence = circular_fragment(start)

                sequence = (
                    forward_sequence
                    if strand == "+"
                    else reverse_complement(forward_sequence)
                )

                end = ((start + LENGTH - 2) % N) + 1
                read_id = f"read_{counter:06d}"

                assert len(sequence) == LENGTH

                fasta.write(
                    f">{read_id}\n{sequence}\n"
                )

                writer.writerow([
                    read_id,
                    group,
                    start,
                    end,
                    LENGTH,
                    strand,
                    REPLICATE
                ])

assert counter == 3000

print("========== TWO-STRAND DATASET ==========")
print("Fragment length:", LENGTH)
print("Damage:", "none")
print("Seed:", SEED)
print("Replicate:", REPLICATE)

for group in groups:
    print(
        f"{group}: "
        f"{PER_STRAND} forward + "
        f"{PER_STRAND} reverse"
    )

print("Total reads:", counter)
print("Forward reads:", 3 * PER_STRAND)
print("Reverse reads:", 3 * PER_STRAND)
print("FASTA:", fasta_path)
print("Ground truth:", truth_path)
print("Generation: PASSED")
