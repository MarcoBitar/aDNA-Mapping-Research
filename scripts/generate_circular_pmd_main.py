import csv
import gzip
import subprocess
from pathlib import Path

MANIFEST = Path("simulation/circular_pmd_main/damage_manifest.tsv")
BASELINE_DIR = Path("simulation/circular_length_baseline")

def read_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    records = {}

    with opener(path, "rt") as handle:
        name = None
        sequence = []

        for line in handle:
            line = line.strip()

            if line.startswith(">"):
                if name is not None:
                    if name in records:
                        raise ValueError(f"Duplicate read ID: {name}")
                    records[name] = "".join(sequence).upper()

                name = line[1:].split()[0]
                sequence = []
            elif line:
                if name is None:
                    raise ValueError(f"Invalid FASTA: {path}")
                sequence.append(line)

        if name is not None:
            if name in records:
                raise ValueError(f"Duplicate read ID: {name}")
            records[name] = "".join(sequence).upper()

    return records

with MANIFEST.open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

assert len(rows) == 12
assert len({int(row["seed"]) for row in rows}) == 12

for row in rows:
    length = int(row["length"])
    level = row["damage_level"]
    seed = int(row["seed"])
    parameters = row["parameters"]

    input_path = BASELINE_DIR / f"reads_{length}bp_none.fasta"
    output_path = Path(row["output"])

    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    if not output_path.exists():
        output_path.parent.mkdir(parents=True, exist_ok=True)

        subprocess.run(
            [
                "deamSim",
                "--seed", str(seed),
                "-damage", parameters,
                "-o", str(output_path),
                str(input_path),
            ],
            check=True,
        )

        print(f"GENERATED: {output_path}")
    else:
        print(f"EXISTING: {output_path}")

    original = read_fasta(input_path)
    damaged = read_fasta(output_path)

    assert len(original) == 6 * (length - 1)
    assert set(original) == set(damaged)

    for read_id, before in original.items():
        after = damaged[read_id]

        assert len(before) == len(after) == length
        assert set(after).issubset(set("ACGT"))

        for a, b in zip(before, after):
            if a != b:
                assert (a, b) in {("C", "T"), ("G", "A")}

    print(f"VALIDATED: {length}bp | {level} | seed={seed}")

print("ALL 12 ORIGINAL PMD DATASETS VALIDATED")
