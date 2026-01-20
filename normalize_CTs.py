import pandas as pd
import os
from pathlib import Path

# Define your reference genes
reference_genes = [
    'UBC21', 'PP2A-1', 'YLS8', 'PP2AA3', 'RHIP1', 
    'His3.3', 'AT5G12240', 'MON1', 'F-Box', 'UPL7'
]

# Get input and output folders
input_folder = input("Enter path to input folder containing TPM tables: ").strip()
output_folder = input("Enter path to output folder for normalized tables: ").strip()

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Process each species table
for file_path in Path(input_folder).glob("*.txt"):  
    print(f"Processing: {file_path.name}")
    
    # Read the count table
    df = pd.read_csv(file_path, index_col=0, sep="\t")
    
    # Extract reference gene rows (including duplicated ones if multiple homologs exist)
    reference_df = df.loc[df.index.intersection(reference_genes)].copy()
    
    if reference_df.empty:
        print(f"Warning: No reference genes found in {file_path.name}")
        continue

    # For each sample (column), compute the mean TPM across reference genes (taking duplicates naturally)
    reference_means = reference_df.sum(axis=0) / len(reference_df)

    # Normalize the entire table: divide each value by the mean reference expression of its column
    normalized_df = df.div(reference_means, axis=1)

    # Save the normalized table
    output_file = Path(output_folder) / file_path.name
    normalized_df.to_csv(output_file, sep="\t")
    
    print(f"Saved normalized table: {output_file.name}")

print("Normalization complete for all species.")
