import csv
import gzip
import os
import subprocess
from pathlib import Path

LENGTHS = [30, 40, 50, 70]

DAMAGE = {
    "low": "0.03,0.30,0.005,0.05",
    "moderate": "0.03,0.30,0.01,0.10",
    "high": "0.03,0.30,0.03,0.30",
}

LEVEL_CODES = {
    "low": 1,
    "moderate": 2,
    "high": 3,
}

REPLICATES = range(2, 11)

INPUT_DIR = Path("simulation/circular_length_baseline")
OUTPUT_DIR = Path("simulation/circular_pmd_replicates")
MANIFEST = OUTPUT_DIR / "damage_replicates_manifest.tsv"

FIELDS = [
    "replicate", "length", "damage_level",
    "parameters", "seed", "input", "output"
]

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def read_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    records = {}

    with opener(path, "rt") as handle:
        name = None
        sequence = []

        for line in handle:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    if name in records:
                        raise ValueError(f"Duplicate ID: {name}")
                    records[name] = "".join(sequence).upper()

                name = line[1:].split()[0]
                sequence = []
            else:
                if name is None:
                    raise ValueError(f"Invalid FASTA: {path}")
                sequence.append(line)

        if name is not None:
            if name in records:
                raise ValueError(f"Duplicate ID: {name}")
            records[name] = "".join(sequence).upper()

    return records


def validate_output(original, output_path, length):
    damaged = read_fasta(output_path)

    assert set(original) == set(damaged), (
        f"Read identifiers differ: {output_path}"
    )

    damaged_count = 0
    substitutions = 0

    for read_id, before in original.items():
        after = damaged[read_id]

        assert len(before) == len(after) == length
        assert set(after).issubset(set("ACGT"))

        changed = False

        for a, b in zip(before, after):
            if a != b:
                assert (a, b) in {
                    ("C", "T"),
                    ("G", "A")
                }, (
                    f"Unexpected substitution in "
                    f"{output_path}: {a}>{b}"
                )

                substitutions += 1
                changed = True

        damaged_count += int(changed)

    return len(damaged), damaged_count, substitutions


def write_manifest(rows):
    temporary = MANIFEST.with_suffix(".tsv.tmp")

    with temporary.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=FIELDS, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

    os.replace(temporary, MANIFEST)


# Preflight: establish all planned seeds and paths.
plan = []

for replicate in REPLICATES:
    for length in LENGTHS:
        for level, parameters in DAMAGE.items():
            seed = (
                20000
                + replicate * 10000
                + length * 10
                + LEVEL_CODES[level]
            )

            plan.append({
                "replicate": replicate,
                "length": length,
                "damage_level": level,
                "parameters": parameters,
                "seed": seed,
                "input": str(
                    INPUT_DIR / f"reads_{length}bp_none.fasta"
                ),
                "output": str(
                    OUTPUT_DIR
                    / f"reads_{length}bp_{level}_rep{replicate}.fasta.gz"
                ),
            })

assert len(plan) == 108
assert len({row["seed"] for row in plan}) == 108
assert len({row["output"] for row in plan}) == 108

# Validate existing manifest before resuming.
if MANIFEST.exists():
    with MANIFEST.open() as f:
        previous = list(csv.DictReader(f, delimiter="\t"))

    expected_by_output = {
        row["output"]: row for row in plan
    }

    for row in previous:
        expected = expected_by_output.get(row["output"])
        assert expected is not None, (
            f"Unexpected manifest entry: {row['output']}"
        )

        for field in FIELDS:
            assert str(row[field]) == str(expected[field]), (
                f"Manifest mismatch: {row['output']}, {field}"
            )

completed = []
generated = 0
reused = 0

for item in plan:
    input_path = Path(item["input"])
    output_path = Path(item["output"])
    length = item["length"]

    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    original = read_fasta(input_path)
    assert len(original) == 6 * (length - 1)
    assert all(len(seq) == length for seq in original.values())

    if output_path.exists():
        status = "REUSED"
        reused += 1
    else:
        status = "GENERATED"

        print(
            f"Generating rep={item['replicate']} "
            f"length={length} "
            f"damage={item['damage_level']} "
            f"seed={item['seed']}",
            flush=True
        )

        subprocess.run(
            [
                "deamSim",
                "--seed", str(item["seed"]),
                "-damage", item["parameters"],
                "-o", str(output_path),
                str(input_path),
            ],
            check=True
        )

        generated += 1

    reads, damaged, substitutions = validate_output(
        original, output_path, length
    )

    completed.append(item)
    write_manifest(completed)

    print(
        f"{status:9s} "
        f"rep={item['replicate']:2d} "
        f"{length:2d}bp "
        f"{item['damage_level']:8s} "
        f"reads={reads:3d} "
        f"damaged={damaged:3d} "
        f"substitutions={substitutions:3d}",
        flush=True
    )

assert len(completed) == 108

print("\n========== REPLICATION COMPLETE ==========")
print("Generated:", generated)
print("Reused:", reused)
print("Validated datasets:", len(completed))
print("Manifest:", MANIFEST)
print("ALL REPLICATES VALIDATED")
