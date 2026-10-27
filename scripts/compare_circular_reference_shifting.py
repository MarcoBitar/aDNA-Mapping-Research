import csv
import statistics
from collections import defaultdict
from pathlib import Path

INPUT = Path("results/circular_pmd_replicate_recall.csv")
OUTPUT = Path("results/circular_reference_shifting_comparison.csv")

with INPUT.open() as f:
    rows = list(csv.DictReader(f))

assert len(rows) == 1440

paired = defaultdict(dict)

for row in rows:
    key = (
        int(row["length"]),
        row["damage"],
        row["aligner"],
        row["group"],
        int(row["replicate"]),
    )

    reference = row["reference"]
    assert reference not in paired[key]
    paired[key][reference] = float(row["recall"])

assert len(paired) == 720

comparisons = defaultdict(list)

for key, values in paired.items():
    assert set(values) == {"original", "shifted"}

    length, damage, aligner, group, replicate = key

    difference = (
        values["shifted"] - values["original"]
    )

    comparisons[
        (length, damage, aligner, group)
    ].append((replicate, difference))

assert len(comparisons) == 72

summary = []

for key, values in sorted(comparisons.items()):
    length, damage, aligner, group = key

    values.sort()

    replicates = [rep for rep, _ in values]
    differences = [difference for _, difference in values]

    assert replicates == list(range(1, 11))

    summary.append({
        "length": length,
        "damage": damage,
        "aligner": aligner,
        "group": group,
        "replicates": 10,
        "mean_delta_pp": 100 * statistics.mean(differences),
        "sd_delta_pp": 100 * statistics.stdev(differences),
        "min_delta_pp": 100 * min(differences),
        "max_delta_pp": 100 * max(differences),
    })

assert len(summary) == 72

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

with OUTPUT.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(summary[0])
    )
    writer.writeheader()
    writer.writerows(summary)

print("========== PAIRED REFERENCE COMPARISON ==========")
print("Paired replicate comparisons:", len(paired))
print("Summary conditions:", len(summary))
print("Output:", OUTPUT)

print("\n========== HIGH-PMD SHIFTING EFFECT ==========")

for length in [30, 40, 50, 70]:
    print(f"\n{length} bp")

    for group in [
        "original_boundary",
        "shifted_boundary",
        "interior",
    ]:
        for aligner in ["bwa", "bowtie2"]:
            matching = [
                row for row in summary
                if row["length"] == length
                and row["damage"] == "high"
                and row["group"] == group
                and row["aligner"] == aligner
            ]

            assert len(matching) == 1
            row = matching[0]

            print(
                f"{group:18s} "
                f"{aligner:7s} "
                f"delta={row['mean_delta_pp']:+7.2f} pp "
                f"SD={row['sd_delta_pp']:5.2f} pp"
            )

print("\nPAIRED REFERENCE COMPARISON COMPLETED SUCCESSFULLY")
