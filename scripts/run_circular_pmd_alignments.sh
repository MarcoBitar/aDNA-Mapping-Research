#!/usr/bin/env bash
set -euo pipefail

mkdir -p alignments/circular_pmd_main
mkdir -p logs/circular_pmd_main

LENGTHS=(30 40 50 70)
LEVELS=(low moderate high)
REFERENCES=(original shifted)

RUNS=0
TOTAL_RECORDS=0

for LENGTH in "${LENGTHS[@]}"; do

    EXPECTED=$((6 * (LENGTH - 1)))

    for DAMAGE in "${LEVELS[@]}"; do

        INPUT="simulation/circular_pmd_main/reads_${LENGTH}bp_${DAMAGE}.fasta.gz"

        for REF in "${REFERENCES[@]}"; do

            if [[ "$REF" == "original" ]]; then
                FASTA="reference/rCRS_NC_012920.1.fasta"
                BT2_INDEX="reference/rCRS_bowtie2"
            else
                FASTA="reference/rCRS_shifted_8000.fasta"
                BT2_INDEX="reference/rCRS_shifted_8000_bowtie2"
            fi

            PREFIX="alignments/circular_pmd_main"
            LOGDIR="logs/circular_pmd_main"

            BWA_PREFIX="${PREFIX}/bwa_${REF}_${LENGTH}bp_${DAMAGE}"
            BT2_PREFIX="${PREFIX}/bowtie2_${REF}_${LENGTH}bp_${DAMAGE}"

            echo "========================================"
            echo "${LENGTH}bp | ${DAMAGE} | ${REF}"
            echo "========================================"

            # BWA-aln: same mismatch setting as baseline.
            bwa aln -n 0.04 \
                "$FASTA" "$INPUT" \
                > "${BWA_PREFIX}.sai" \
                2> "${LOGDIR}/bwa_${REF}_${LENGTH}bp_${DAMAGE}_aln.log"

            bwa samse \
                "$FASTA" \
                "${BWA_PREFIX}.sai" \
                "$INPUT" \
                > "${BWA_PREFIX}.sam" \
                2> "${LOGDIR}/bwa_${REF}_${LENGTH}bp_${DAMAGE}_samse.log"

            # Bowtie2: same mode and preset as baseline.
            bowtie2 \
                --end-to-end \
                --sensitive \
                -f \
                -x "$BT2_INDEX" \
                -U "$INPUT" \
                -S "${BT2_PREFIX}.sam" \
                2> "${LOGDIR}/bowtie2_${REF}_${LENGTH}bp_${DAMAGE}.log"

            # Count primary records, excluding supplementary
            # and secondary alignments.
            for SAM in "${BWA_PREFIX}.sam" "${BT2_PREFIX}.sam"; do

                COUNT=$(samtools view -c -F 0x900 "$SAM")

                if [[ "$COUNT" -ne "$EXPECTED" ]]; then
                    echo "ERROR: $SAM has $COUNT primary records"
                    echo "Expected: $EXPECTED"
                    exit 1
                fi

                echo "PASS: $SAM ($COUNT primary records)"

                RUNS=$((RUNS + 1))
                TOTAL_RECORDS=$((TOTAL_RECORDS + COUNT))
            done

        done
    done
done

echo
echo "========== FINAL VALIDATION =========="
echo "Alignment configurations: $RUNS"
echo "Primary alignment records: $TOTAL_RECORDS"

[[ "$RUNS" -eq 48 ]]
[[ "$TOTAL_RECORDS" -eq 13392 ]]

echo "ALL 48 PMD ALIGNMENTS COMPLETED SUCCESSFULLY"
