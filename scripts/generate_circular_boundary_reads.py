import csv
import random
from pathlib import Path

SEED = 1601
LENGTH = 50
PER_GROUP = 1000
SHIFT = 8000

rng = random.Random(SEED)

reference_path = Path(
    "reference/rCRS_NC_012920.1.fasta"
)

fasta_path = Path(
    "simulation/circular_boundary_main/"
    "reads_50bp_none_rep1.fasta"
)

truth_path = Path(
    "ground_truth/"
    "reads_50bp_none_rep1.tsv"
)

with reference_path.open() as f:
    next(f)
    reference = "".join(
        line.strip() for line in f
    ).upper()

N = len(reference)
assert N == 16569
assert 1 < LENGTH < SHIFT < N - LENGTH

def extract_circular(start):
    return "".join(
        reference[(start - 1 + i) % N]
        for i in range(LENGTH)
    )

def crosses_boundary(start, boundary):
    # boundary is the original coordinate of the
    # final base before the artificial cut
    return any(
        ((start - 1 + i) % N) + 1 == boundary
        for i in range(LENGTH - 1)
    )

original_boundary = N
shifted_boundary = SHIFT

groups = {
    "original_boundary": [],
    "shifted_boundary": [],
    "interior": [],
}

for start in range(1, N + 1):
    crosses_original = crosses_boundary(
        start, original_boundary
    )
    crosses_shifted = crosses_boundary(
        start, shifted_boundary
    )

    if crosses_original:
        groups["original_boundary"].append(start)
    elif crosses_shifted:
        groups["shifted_boundary"].append(start)
    else:
        groups["interior"].append(start)

assert len(groups["original_boundary"]) == LENGTH - 1
assert len(groups["shifted_boundary"]) == LENGTH - 1
assert len(groups["interior"]) == N - 2 * (LENGTH - 1)

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

    counter = 0

    for group, starts in groups.items():
        for _ in range(PER_GROUP):
            counter += 1
            start = rng.choice(starts)

            sequence = extract_circular(start)
            end = ((start - 1 + LENGTH - 1) % N) + 1

            read_id = f"read_{counter:06d}"

            fasta.write(
                f">{read_id}\n{sequence}\n"
            )

            writer.writerow([
                read_id,
                group,
                start,
                end,
                LENGTH,
                "+",
                1
            ])

print("========== DATASET SUMMARY ==========")
print("Fragment length:", LENGTH)
print("Damage condition: none")
print("Random seed:", SEED)

for group, starts in groups.items():
    print(
        f"{group}: {PER_GROUP} reads "
        f"({len(starts)} eligible start positions)"
    )

print("Total reads:", counter)
print("FASTA:", fasta_path)
print("Ground truth:", truth_path)

assert counter == 3000
print("Validation: PASSED")
