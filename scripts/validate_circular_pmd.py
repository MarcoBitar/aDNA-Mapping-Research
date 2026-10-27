import csv
import gzip
from collections import Counter
from pathlib import Path

LENGTHS = [30, 40, 50, 70]
LEVELS = ["low", "moderate", "high"]

BASE = Path("simulation/circular_length_baseline")
DAMAGED = Path("simulation/circular_pmd_main")
OUTPUT = Path("results/circular_pmd_validation.csv")

OUTPUT.parent.mkdir(parents=True, exist_ok=True)


def read_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open

    records = {}
    current_id = None
    parts = []

    with opener(path, "rt") as f:
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if current_id is not None:
                    assert current_id not in records
                    records[current_id] = "".join(parts)

                current_id = line[1:].split()[0]
                parts = []

            elif line:
                assert current_id is not None
                parts.append(line)

    if current_id is not None:
        assert current_id not in records
        records[current_id] = "".join(parts)

    return records


results = []

print("========== PMD DATASET VALIDATION ==========")

for length in LENGTHS:

    original = read_fasta(
        BASE / f"reads_{length}bp_none.fasta"
    )

    expected_count = 6 * (length - 1)

    assert len(original) == expected_count

    print(f"\n========== {length} BP ==========")

    for level in LEVELS:

        path = DAMAGED / f"reads_{length}bp_{level}.fasta.gz"

        damaged = read_fasta(path)

        assert len(damaged) == expected_count
        assert set(damaged) == set(original)

        substitutions = Counter()
        damaged_reads = 0
        total_substitutions = 0

        for read_id, before in original.items():

            after = damaged[read_id]

            assert len(before) == length
            assert len(after) == length

            assert set(before) <= set("ACGT")
            assert set(after) <= set("ACGT")

            changed = False

            for a, b in zip(before, after):

                if a != b:
                    mutation = f"{a}>{b}"

                    assert mutation in {"C>T", "G>A"}, (
                        f"Unexpected mutation {mutation} "
                        f"in {read_id}, {length}bp {level}"
                    )

                    substitutions[mutation] += 1
                    total_substitutions += 1
                    changed = True

            if changed:
                damaged_reads += 1

        fraction = damaged_reads / expected_count

        results.append({
            "length": length,
            "damage_level": level,
            "reads": expected_count,
            "damaged_reads": damaged_reads,
            "damaged_fraction": fraction,
            "substitutions": total_substitutions,
            "C_to_T": substitutions["C>T"],
            "G_to_A": substitutions["G>A"]
        })

        print(
            f"{level:8s} "
            f"reads={expected_count:3d} "
            f"damaged={damaged_reads:3d} "
            f"fraction={fraction:6.2%} "
            f"C>T={substitutions['C>T']:3d} "
            f"G>A={substitutions['G>A']:3d} "
            f"total_changes={total_substitutions:3d} "
            f"PASS"
        )

with OUTPUT.open("w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "length",
            "damage_level",
            "reads",
            "damaged_reads",
            "damaged_fraction",
            "substitutions",
            "C_to_T",
            "G_to_A"
        ]
    )

    writer.writeheader()
    writer.writerows(results)

assert len(results) == 12

print("\nValidation rows:", len(results))
print("Results saved:", OUTPUT)
print("ALL 12 PMD DATASETS VALIDATED")
