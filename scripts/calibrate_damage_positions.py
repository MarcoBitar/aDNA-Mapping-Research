import argparse
import gzip
from collections import defaultdict
from itertools import zip_longest


def read_fasta(path):
    """Yield FASTA records in their original order."""

    with gzip.open(path, "rt") as handle:
        name = None
        parts = []

        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()

                name = line[1:].split()[0]
                parts = []
            else:
                if name is None:
                    raise ValueError("Sequence before FASTA header")

                parts.append(line)

        if name is not None:
            yield name, "".join(parts).upper()


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--original", required=True)
    parser.add_argument("--damaged", required=True)

    args = parser.parse_args()

    stats = defaultdict(
        lambda: {
            "C": 0,
            "CT": 0,
            "G": 0,
            "GA": 0,
        }
    )

    fragment_count = 0

    original = read_fasta(args.original)
    damaged = read_fasta(args.damaged)

    for record_number, (original_record, damaged_record) in enumerate(
        zip_longest(original, damaged),
        start=1
    ):
        if original_record is None or damaged_record is None:
            raise ValueError(
                f"Different record counts at record {record_number}"
            )

        original_id, before = original_record
        damaged_id, after = damaged_record

        if original_id != damaged_id:
            raise ValueError(
                f"Identifier mismatch at record {record_number}"
            )

        if len(before) != len(after):
            raise ValueError(
                f"Sequence length mismatch at record {record_number}"
            )

        fragment_count += 1

        for i, (a, b) in enumerate(zip(before, after)):
            position = i + 1

            if a == "C":
                stats[position]["C"] += 1

                if b == "T":
                    stats[position]["CT"] += 1

            if a == "G":
                stats[position]["G"] += 1

                if b == "A":
                    stats[position]["GA"] += 1

    if fragment_count == 0:
        raise ValueError("No fragments found")

    print(
        "Position\tC_count\tC_to_T\tC_to_T_rate\t"
        "G_count\tG_to_A\tG_to_A_rate"
    )

    for position in sorted(stats):
        s = stats[position]

        ct_rate = (
            s["CT"] / s["C"]
            if s["C"] else 0
        )

        ga_rate = (
            s["GA"] / s["G"]
            if s["G"] else 0
        )

        print(
            f"{position}\t"
            f"{s['C']}\t"
            f"{s['CT']}\t"
            f"{ct_rate:.6f}\t"
            f"{s['G']}\t"
            f"{s['GA']}\t"
            f"{ga_rate:.6f}"
        )


if __name__ == "__main__":
    main()
