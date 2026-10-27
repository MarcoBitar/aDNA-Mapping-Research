import sys
import pysam

def evaluate(sam_path):
    total = 0
    mapped = 0
    correct = 0
    incorrect = 0

    with pysam.AlignmentFile(sam_path, "r") as sam:
        for read in sam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue

            total += 1

            if read.is_unmapped:
                print(f"{read.query_name}: UNMAPPED")
                continue

            mapped += 1

            # Gargammel identifier:
            # NC_012920.1:+:10891:10941:50e-1
            parts = read.query_name.split(":")

            true_reference = parts[0]
            true_strand = parts[1]
            true_start = int(parts[2])
            true_end = int(parts[3])

            # pysam uses zero-based reference_start.
            observed_reference = read.reference_name
            observed_start = read.reference_start
            observed_strand = "-" if read.is_reverse else "+"

            expected_length = true_end - true_start

            is_correct = (
                observed_reference == true_reference
                and observed_start == true_start
                and observed_strand == true_strand
                and read.cigarstring == f"{expected_length}M"
            )

            if is_correct:
                correct += 1
            else:
                incorrect += 1

            print(
                f"{read.query_name}: "
                f"position={observed_start + 1}, "
                f"strand={observed_strand}, "
                f"correct={is_correct}"
            )

    print("\nEVALUATION SUMMARY")
    print("------------------")
    print("Total reads:", total)
    print("Mapped reads:", mapped)
    print("Correct alignments:", correct)
    print("Incorrect mapped alignments:", incorrect)
    print("Unmapped reads:", total - mapped)
    print(
        "Mapping rate:",
        f"{100 * mapped / total:.2f}%" if total else "N/A"
    )
    print(
        "Correct mapping rate:",
        f"{100 * correct / total:.2f}%" if total else "N/A"
    )

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python scripts/evaluate_pilot.py <SAM file>")

    evaluate(sys.argv[1])
