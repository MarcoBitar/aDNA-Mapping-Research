import argparse
import gzip
from collections import Counter
from itertools import zip_longest


def read_fasta(path):
    """Yield FASTA records in their original order."""

    with gzip.open(path, "rt") as handle:
        name = None
        sequence_parts = []

        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(sequence_parts).upper()

                name = line[1:].split()[0]
                sequence_parts = []

            else:
                if name is None:
                    raise ValueError(
                        "FASTA sequence found before a header"
                    )

                sequence_parts.append(line)

        if name is not None:
            yield name, "".join(sequence_parts).upper()


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--original", required=True)
    parser.add_argument("--damaged", required=True)

    args = parser.parse_args()

    original = read_fasta(args.original)
    damaged = read_fasta(args.damaged)

    substitutions = Counter()

    fragment_count = 0
    total_bases = 0
    changed_fragments = 0

    for record_number, pair in enumerate(
        zip_longest(original, damaged),
        start=1
    ):
        original_record, damaged_record = pair

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
        fragment_changed = False

        for a, b in zip(before, after):
            total_bases += 1

            if a != b:
                substitutions[f"{a}>{b}"] += 1
                fragment_changed = True

        if fragment_changed:
            changed_fragments += 1

    if fragment_count == 0:
        raise ValueError("No fragments found")

    total_substitutions = sum(substitutions.values())

    print("========== DAMAGE CALIBRATION ==========")
    print("Fragments:", fragment_count)
    print("Total bases:", total_bases)
    print("Fragments with substitutions:", changed_fragments)
    print("Total substitutions:", total_substitutions)

    print("\n========== SUBSTITUTION COUNTS ==========")

    if substitutions:
        for change, count in sorted(substitutions.items()):
            print(f"{change}: {count}")
    else:
        print("No substitutions detected")

    print("\n========== SUMMARY ==========")

    print(
        "Overall substitution fraction:",
        f"{total_substitutions / total_bases:.4%}"
    )

    print(
        "Fragments with substitutions:",
        f"{changed_fragments / fragment_count:.4%}"
    )

    print("C>T substitutions:", substitutions["C>T"])
    print("G>A substitutions:", substitutions["G>A"])


if __name__ == "__main__":
    main()
