import csv
import statistics
from collections import defaultdict
from pathlib import Path

INPUT = Path("results/circular_pmd_replicates_all.csv")
OUTPUT = Path("results/circular_pmd_replicate_summary.csv")
REPLICATE_OUTPUT = Path("results/circular_pmd_replicate_recall.csv")

LENGTHS = [30, 40, 50, 70]
DAMAGE = ["low", "moderate", "high"]
ALIGNERS = ["bwa", "bowtie2"]
REFERENCES = ["original", "shifted"]
GROUPS = ["original_boundary", "shifted_boundary", "interior"]

with INPUT.open() as f:
    rows = list(csv.DictReader(f))

assert len(rows) == 2976

# Combine forward and reverse strands within each replicate.
grouped = defaultdict(lambda: {
    "correct": 0,
    "incorrect": 0,
    "unmapped": 0,
    "strands": set(),
})

for row in rows:
    if row["damage"] == "none":
        continue

    key = (
        int(row["length"]),
        row["damage"],
        row["aligner"],
        row["reference"],
        row["group"],
        int(row["replicate"]),
    )

    item = grouped[key]

    assert row["strand"] not in item["strands"]
    item["strands"].add(row["strand"])

    for metric in ["correct", "incorrect", "unmapped"]:
        item[metric] += int(row[metric])

assert len(grouped) == 1440

replicate_rows = []
summary_groups = defaultdict(list)

for key, item in sorted(grouped.items()):
    length, damage, aligner, reference, group, replicate = key

    assert item["strands"] == {"+", "-"}

    total = (
        item["correct"]
        + item["incorrect"]
        + item["unmapped"]
    )

    assert total == 2 * (length - 1)

    recall = item["correct"] / total

    replicate_rows.append({
        "length": length,
        "damage": damage,
        "aligner": aligner,
        "reference": reference,
        "group": group,
        "replicate": replicate,
        "correct": item["correct"],
        "incorrect": item["incorrect"],
        "unmapped": item["unmapped"],
        "total": total,
        "recall": recall,
    })

    summary_groups[key[:-1]].append((replicate, recall))

assert len(replicate_rows) == 1440

summary_rows = []

for key, values in sorted(summary_groups.items()):
    length, damage, aligner, reference, group = key

    values.sort()
    replicates = [rep for rep, _ in values]
    recalls = [recall for _, recall in values]

    assert replicates == list(range(1, 11))

    mean_recall = statistics.mean(recalls)
    sd_recall = statistics.stdev(recalls)

    summary_rows.append({
        "length": length,
        "damage": damage,
        "aligner": aligner,
        "reference": reference,
        "group": group,
        "replicates": len(recalls),
        "mean_recall": mean_recall,
        "sd_recall": sd_recall,
        "min_recall": min(recalls),
        "max_recall": max(recalls),
    })

assert len(summary_rows) == 144

REPLICATE_OUTPUT.parent.mkdir(parents=True, exist_ok=True)

with REPLICATE_OUTPUT.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(replicate_rows[0])
    )
    writer.writeheader()
    writer.writerows(replicate_rows)

with OUTPUT.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(summary_rows[0])
    )
    writer.writeheader()
    writer.writerows(summary_rows)

print("========== REPLICATE STATISTICS ==========")
print("Replicate-level rows:", len(replicate_rows))
print("Summary conditions:", len(summary_rows))
print("Replicates per condition: 10")
print("Replicate-level CSV:", REPLICATE_OUTPUT)
print("Summary CSV:", OUTPUT)

print("\n========== HIGH-PMD RECALL SUMMARY ==========")

for length in LENGTHS:
    print(f"\n{length} bp")

    for group in GROUPS:
        print(f"  {group}")

        for reference in REFERENCES:
            for aligner in ALIGNERS:
                matching = [
                    r for r in summary_rows
                    if r["length"] == length
                    and r["damage"] == "high"
                    and r["group"] == group
                    and r["reference"] == reference
                    and r["aligner"] == aligner
                ]

                assert len(matching) == 1
                r = matching[0]

                print(
                    f"    {reference:8s} "
                    f"{aligner:7s} "
                    f"mean={100*r['mean_recall']:6.2f}% "
                    f"SD={100*r['sd_recall']:5.2f}%"
                )

print("\nREPLICATE STATISTICS COMPLETED SUCCESSFULLY")
