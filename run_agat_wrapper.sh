#!/bin/bash
set -euo pipefail

# Directories
GENOME_DIR="/path/to/genome/directory/"
GFF_DIR="/path/to/gff3/directory/"
CDS_DIR="/path/to/output/cds/directory/"
SINGLEISO_GFF_DIR="/path/to/output/singleiso_gff/directory/"
SINGLEISO_CDS_DIR="/path/to/output/singleiso_cds/directory/"
LOG_DIR="/path/to/output/logs/directory/"

# Make sure output dirs exist
mkdir -p "$CDS_DIR" "$SINGLEISO_GFF_DIR" "$SINGLEISO_CDS_DIR" "$LOG_DIR"

# Loop through all genome files
for genome_file in "$GENOME_DIR"/*.fa; do
    species=$(basename "$genome_file" .genome.fa)
    gff_file="$GFF_DIR/$species.gff3"

    if [[ ! -f "$gff_file" ]]; then
        echo "No matching GFF for $species" | tee -a "$LOG_DIR/errors.log"
        continue
    fi

    echo "Processing $species ..."
    {
        # Step 1: Extract CDS from original GFF
        agat_sp_extract_sequences.pl \
            -g "$gff_file" \
            -f "$genome_file" \
            -o "$CDS_DIR/$species.cds.fa"

        # Step 2: Keep only longest isoform
        agat_sp_keep_longest_isoform.pl \
            -g "$gff_file" \
            -o "$SINGLEISO_GFF_DIR/$species.gff3"

        # Step 3: Extract CDS again from single isoform GFF
        agat_sp_extract_sequences.pl \
            -g "$SINGLEISO_GFF_DIR/$species.gff3" \
            -f "$genome_file" \
            -o "$SINGLEISO_CDS_DIR/$species.cds.fa"

        echo "Finished $species"
    } >> "$LOG_DIR/${species}.log" 2>&1 || {
        echo "Failed $species (see $LOG_DIR/${species}.log)" | tee -a "$LOG_DIR/errors.log"
    }
done
