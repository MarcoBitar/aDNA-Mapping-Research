#!/usr/bin/env bash
set -euo pipefail

ORIGINAL="reference/rCRS_NC_012920.1.fasta"
SHIFTED="reference/rCRS_shifted_8000.fasta"

for FASTA in "$ORIGINAL" "$SHIFTED"; do
    if [[ ! -s "$FASTA" ]]; then
        echo "ERROR: Missing reference: $FASTA"
        exit 1
    fi
done

echo "Building BWA indexes..."
bwa index "$ORIGINAL"
bwa index "$SHIFTED"

echo "Building Bowtie2 indexes..."
bowtie2-build "$ORIGINAL" reference/rCRS_bowtie2
bowtie2-build "$SHIFTED" reference/rCRS_shifted_8000_bowtie2

echo "Building FASTA indexes..."
samtools faidx "$ORIGINAL"
samtools faidx "$SHIFTED"

echo "Reference indexes built successfully."
