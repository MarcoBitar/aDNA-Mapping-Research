#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR="simulation/circular_pmd_replicates"
OUTPUT_DIR="alignments/circular_pmd_replicates"
LOG_DIR="logs/circular_pmd_replicates"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

LENGTHS=(30 40 50 70)
LEVELS=(low moderate high)
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

    # Check SAM readability and primary-record count.
    COUNT=$(samtools view -c -F 0x900 "$SAM") || return 1

    [[ "$COUNT" -eq "$EXPECTED" ]] || return 1
}

run_bwa() {
    local FASTA="$1"
    local INPUT="$2"
    local PREFIX="$3"
    local LOG_PREFIX="$4"
    local EXPECTED="$5"

    if validate_sam "${PREFIX}.sam" "$EXPECTED"; then
        echo "REUSED: ${PREFIX}.sam"
        REUSED=$((REUSED + 1))
        return
    fi

    # Do not silently replace an existing invalid alignment.
    if [[ -e "${PREFIX}.sam" ]]; then
        echo "ERROR: Existing BWA SAM failed validation:"
        echo "${PREFIX}.sam"
        exit 1
    fi

    local TEMP_SAI="${PREFIX}.sai.tmp"
    local TEMP_SAM="${PREFIX}.sam.tmp"

    rm -f "$TEMP_SAI" "$TEMP_SAM"

    bwa aln -n 0.04 \
        "$FASTA" "$INPUT" \
        > "$TEMP_SAI" \
        2> "${LOG_PREFIX}_aln.log"

    bwa samse \
        "$FASTA" "$TEMP_SAI" "$INPUT" \
        > "$TEMP_SAM" \
        2> "${LOG_PREFIX}_samse.log"

    if ! validate_sam "$TEMP_SAM" "$EXPECTED"; then
        echo "ERROR: BWA output failed validation:"
        echo "$TEMP_SAM"
        exit 1
    fi

    mv "$TEMP_SAI" "${PREFIX}.sai"
    mv "$TEMP_SAM" "${PREFIX}.sam"

    echo "GENERATED: ${PREFIX}.sam"
    GENERATED=$((GENERATED + 1))
}

run_bowtie2() {
    local INDEX="$1"
    local INPUT="$2"
    local PREFIX="$3"
    local LOG_FILE="$4"
    local EXPECTED="$5"

    if validate_sam "${PREFIX}.sam" "$EXPECTED"; then
        echo "REUSED: ${PREFIX}.sam"
        REUSED=$((REUSED + 1))
        return
    fi

    if [[ -e "${PREFIX}.sam" ]]; then
        echo "ERROR: Existing Bowtie2 SAM failed validation:"
        echo "${PREFIX}.sam"
        exit 1
    fi

    local TEMP_SAM="${PREFIX}.sam.tmp"

    rm -f "$TEMP_SAM"

    bowtie2 \
        --end-to-end \
        --sensitive \
        -f \
        -x "$INDEX" \
        -U "$INPUT" \
        -S "$TEMP_SAM" \
        2> "$LOG_FILE"

    if ! validate_sam "$TEMP_SAM" "$EXPECTED"; then
        echo "ERROR: Bowtie2 output failed validation:"
        echo "$TEMP_SAM"
        exit 1
    fi

    mv "$TEMP_SAM" "${PREFIX}.sam"

    echo "GENERATED: ${PREFIX}.sam"
    GENERATED=$((GENERATED + 1))
}

for REP in {2..10}; do
    for LENGTH in "${LENGTHS[@]}"; do
        EXPECTED=$((6 * (LENGTH - 1)))

        for DAMAGE in "${LEVELS[@]}"; do
            INPUT="${INPUT_DIR}/reads_${LENGTH}bp_${DAMAGE}_rep${REP}.fasta.gz"

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

                LABEL="${REF}_${LENGTH}bp_${DAMAGE}_rep${REP}"

                echo "========================================"
                echo "Replicate $REP | ${LENGTH}bp | $DAMAGE | $REF"
                echo "========================================"

                run_bwa \
                    "$FASTA" "$INPUT" \
                    "${OUTPUT_DIR}/bwa_${LABEL}" \
                    "${LOG_DIR}/bwa_${LABEL}" \
                    "$EXPECTED"

                run_bowtie2 \
                    "$BT2_INDEX" "$INPUT" \
                    "${OUTPUT_DIR}/bowtie2_${LABEL}" \
                    "${LOG_DIR}/bowtie2_${LABEL}.log" \
                    "$EXPECTED"

                RUNS=$((RUNS + 2))
                TOTAL_RECORDS=$((TOTAL_RECORDS + 2 * EXPECTED))
            done
        done
    done
done

echo
echo "========== REPLICATED ALIGNMENT VALIDATION =========="
echo "Alignment configurations: $RUNS"
echo "Newly generated: $GENERATED"
echo "Reused: $REUSED"
echo "Primary alignment records: $TOTAL_RECORDS"

[[ "$RUNS" -eq 432 ]]
[[ "$TOTAL_RECORDS" -eq 120528 ]]
[[ "$((GENERATED + REUSED))" -eq 432 ]]

echo "ALL 432 REPLICATED ALIGNMENTS COMPLETED SUCCESSFULLY"
