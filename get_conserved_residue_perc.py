import os
import argparse

def process_kipes_results(kipes_dir, input_file, output_file, target_gene):
    # Read the input file and store species-gene mappings
    species_gene_map = {}
    with open(input_file, 'r') as infile:
        for line in infile:
            line = line.strip()
            if not line:
                continue
            parts = line.split('_')
            species_name = '_'.join(parts[:2])  # First two words as species name
            gene_id = '_'.join(parts[2:])  # The rest as gene ID
            species_gene_map[line] = (species_name, gene_id)

    # Process KIPEs result directory
    with open(output_file, 'w') as outfile:
        for original_id, (species_name, gene_id) in species_gene_map.items():
            species_dir = os.path.join(kipes_dir, species_name)
            summary_file = os.path.join(species_dir, "summary.txt")

            if not os.path.exists(summary_file):
                print(f"Warning: {summary_file} not found. Skipping {species_name}.")
                continue

            with open(summary_file, 'r') as summary:
                for line in summary:
                    columns = line.strip().split('\t')
                    if len(columns) < 5:
                        continue  # Skip malformed lines
                    
                    seq_id, gene, _, conserved_residues, _ = columns
                    
                    # Extract the base gene name before underscore
                    base_gene_name = gene.split('_')[0]

                    if base_gene_name == target_gene and seq_id == gene_id:
                        conserved_residues = float(conserved_residues)  
                        diff_value = 100 - conserved_residues
                        outfile.write(f"{original_id},{conserved_residues}\n")
                        break  

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract conserved residue information from KIPEs summary files.")
    parser.add_argument("-k", "--kipes_dir", required=True, help="Path to KIPEs result directory.")
    parser.add_argument("-i", "--input", required=True, help="Input file containing FASTA headers.")
    parser.add_argument("-o", "--output", required=True, help="Output file to store results.")
    parser.add_argument("-g", "--gene", required=True, help="Target gene name to search for in summary.txt.")

    args = parser.parse_args()
    process_kipes_results(args.kipes_dir, args.input, args.output, args.gene)
