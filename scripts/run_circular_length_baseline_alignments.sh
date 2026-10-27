#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR="simulation/circular_length_baseline"
OUTPUT_DIR="alignments/circular_length_baseline"
LOG_DIR="logs/circular_length_baseline"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

LENGTHS=(30 40 50 70)
REFERENCES=(original shifted)

RUNS=0
TOTAL_RECORDS=0
GENERATED=0
REUSED=0

validate_sam() {
    local SAM="$1"
    local EXPECTED="$2"
    local COUNT

    [[ -s "$SAM" ]] || return 1
    COUNT=$(samtools view -c -F 0x900 "$SAM") || return 1
    [[ "$COUNT" -eq "$EXPECTED" ]]
}

for LENGTH in "${LENGTHS[@]}"; do
    INPUT="${INPUT_DIR}/reads_${LENGTH}bp_none.fasta"
    EXPECTED=$((6 * (LENGTH - 1)))

    if [[ ! -s "$INPUT" ]]; then
        echo "ERROR: Missing input: $INPUT"
        exit 1
    fi

    for REF in "${REFERENCES[@]}"; do
        if [[ "$REF" == "original" ]]; then
            FASTA="reference/rCRS_NC_012920.1.fasta"
            BT2_INDEX="reference/rCRS_bowtie2"
        else
            FASTA="reference/rCRS_shifted_8000.fasta"
            BT2_INDEX="reference/rCRS_shifted_8000_bowtie2"
        fi

        LABEL="${REF}_${LENGTH}bp_none"

        BWA_PREFIX="${OUTPUT_DIR}/bwa_${LABEL}"
        BT2_PREFIX="${OUTPUT_DIR}/bowtie2_${LABEL}"

        echo "=== ${LENGTH}bp | ${REF} ==="

        if validate_sam "${BWA_PREFIX}.sam" "$EXPECTED"; then
            echo "REUSED: ${BWA_PREFIX}.sam"
            REUSED=$((REUSED + 1))
        else
            if [[ -e "${BWA_PREFIX}.sam" ]]; then
                echo "ERROR: Invalid existing SAM: ${BWA_PREFIX}.sam"
                exit 1
            fi

            TEMP_SAI="${BWA_PREFIX}.sai.tmp"
            TEMP_SAM="${BWA_PREFIX}.sam.tmp"

            bwa aln -n 0.04 "$FASTA" "$INPUT" \
                > "$TEMP_SAI" \
                2> "${LOG_DIR}/bwa_${LABEL}_aln.log"

            bwa samse "$FASTA" "$TEMP_SAI" "$INPUT" \
                > "$TEMP_SAM" \
                2> "${LOG_DIR}/bwa_${LABEL}_samse.log"

            if ! validate_sam "$TEMP_SAM" "$EXPECTED"; then
                echo "ERROR: BWA validation failed: $TEMP_SAM"
                exit 1
            fi

            mv "$TEMP_SAI" "${BWA_PREFIX}.sai"
            mv "$TEMP_SAM" "${BWA_PREFIX}.sam"

            GENERATED=$((GENERATED + 1))
            echo "GENERATED: ${BWA_PREFIX}.sam"
        fi

        if validate_sam "${BT2_PREFIX}.sam" "$EXPECTED"; then
            echo "REUSED: ${BT2_PREFIX}.sam"
            REUSED=$((REUSED + 1))
        else
            if [[ -e "${BT2_PREFIX}.sam" ]]; then
                echo "ERROR: Invalid existing SAM: ${BT2_PREFIX}.sam"
                exit 1
            fi

            TEMP_SAM="${BT2_PREFIX}.sam.tmp"

            bowtie2 \
                --end-to-end \
                --sensitive \
                -f \
                -x "$BT2_INDEX" \
                -U "$INPUT" \
                -S "$TEMP_SAM" \
                2> "${LOG_DIR}/bowtie2_${LABEL}.log"

            if ! validate_sam "$TEMP_SAM" "$EXPECTED"; then
                echo "ERROR: Bowtie2 validation failed: $TEMP_SAM"
                exit 1
            fi

            mv "$TEMP_SAM" "${BT2_PREFIX}.sam"

            GENERATED=$((GENERATED + 1))
            echo "GENERATED: ${BT2_PREFIX}.sam"
        fi

        RUNS=$((RUNS + 2))
        TOTAL_RECORDS=$((TOTAL_RECORDS + 2 * EXPECTED))
    done
done

echo
echo "=== BASELINE ALIGNMENT VALIDATION ==="
echo "Configurations: $RUNS"
echo "Generated: $GENERATED"
echo "Reused: $REUSED"
echo "Primary records: $TOTAL_RECORDS"

[[ "$RUNS" -eq 16 ]]
[[ "$TOTAL_RECORDS" -eq 4464 ]]
[[ "$((GENERATED + REUSED))" -eq 16 ]]

echo "ALL 16 BASELINE ALIGNMENTS VALIDATED"
