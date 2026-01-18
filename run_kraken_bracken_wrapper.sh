#!/bin/bash

# Define directories
input_dir="path/to/coding/sequences/directory/"
output_base="/path/to/output/directory/"
db_path="/path/to/database/"

# Step 1: Run Kraken2 sequentially on each .cds.fasta file
for file in "$input_dir"/*.cds.fasta; do
    # Extract the base name (without path and extension)
    base_name=$(basename "$file" .cds.fasta)
    
    # Define output folder and ensure it exists
    output_folder="$output_base/$base_name"
    mkdir -p "$output_folder"

    # Run Kraken2
    echo "Processing $file with Kraken2..."
    kraken2 --db "$db_path" \
            --report "$output_folder/kraken2.kreport" \
            --output "$output_folder/kraken2.kraken" \
            --use-names "$file" >> "$output_folder/kraken2.log" 2>&1

    # Wait for Kraken2 to finish before proceeding to the next file
    wait
done

echo "Kraken2 processing complete for all files."

# Step 2: Run Bracken on all generated kraken2.kreport files
for report in "$output_base"/*/kraken2.kreport; do
    # Extract folder name where the report is located
    output_folder=$(dirname "$report")

    # Run Bracken
    echo "Processing $report with Bracken..."
    bracken -d "$db_path" \
            -i "$report" \
            -o "$output_folder/kraken2.bracken" \
            -l S

    # Wait for Bracken to finish before proceeding to the next report
    wait
done

echo "Bracken processing complete for all files."
