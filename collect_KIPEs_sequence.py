import os
import argparse
from pathlib import Path

def parse_summary(summary_file, enzymes, threshold):
    """
    Parses summary.txt and returns a dictionary of valid sequence IDs.
    Ensures that enzyme names match correctly and each sequence appears only once.
    """
    valid_sequences = {}

    with open(summary_file, 'r') as f:
        for line in f:
            columns = line.strip().split('\t')
            if len(columns) < 5:
                continue  # Skip malformed lines
            
            seq_id, enzyme, _, conserved_residues, _ = columns

            try:
                conserved_residues = float(conserved_residues)
            except ValueError:
                continue  # Skip invalid numerical values

            # Extract base enzyme name (FLS_1 → FLS)
            enzyme_base = enzyme.split('_')[0]

            if enzyme_base in enzymes and conserved_residues >= threshold:
                valid_sequences[seq_id] = enzyme_base  # Store only the first valid enzyme

    print(f"Found {len(valid_sequences)} valid sequences in {summary_file}")
    return valid_sequences

def process_kipes_results(main_dir, enzymes, output_file, threshold):
    """
    Processes KIPEs result folder, filters sequences by conserved residues, and concatenates valid sequences.
    """
    main_dir = Path(main_dir)
    if not main_dir.is_dir():
        print(f"Error: {main_dir} is not a valid directory.")
        return
    
    written_sequences = set()  # Ensure unique sequence IDs in the output
    skipped_sequences = {}  # To track how many times a sequence ID is skipped

    with open(output_file, 'w') as out_fh:
        for sub_dir in main_dir.iterdir():
            if sub_dir.is_dir():
                summary_file = sub_dir / "summary.txt"
                final_pep_files_dir = sub_dir / "final_pep_files"
                
                if not summary_file.exists() or not final_pep_files_dir.is_dir():
                    continue  # Skip if required files/directories are missing
                
                # Get valid sequence IDs based on summary.txt filtering
                valid_sequences = parse_summary(summary_file, enzymes, threshold)

                if not valid_sequences:
                    print(f"Warning: No valid sequences found in {summary_file}")
                    continue

                for fasta_file in final_pep_files_dir.iterdir():
                    if fasta_file.suffix == ".fasta":
                        with open(fasta_file, 'r') as enzyme_fh:
                            write_sequence = False
                            for line in enzyme_fh:
                                if line.startswith(">"):
                                    seq_id = line.strip().split()[0][1:]  # Extract sequence ID
                                    
                                    # Write only if this ID meets the threshold and hasn't been written yet
                                    if seq_id in valid_sequences and seq_id not in written_sequences:
                                        write_sequence = True
                                        written_sequences.add(seq_id)  # Mark it as written
                                        species_name = sub_dir.name
                                        out_fh.write(f">{species_name}_{seq_id}\n")
                                    else:
                                        write_sequence = False
                                        # Track and print duplicate sequence IDs
                                        if seq_id in written_sequences:
                                            if seq_id not in skipped_sequences:
                                                skipped_sequences[seq_id] = []
                                            skipped_sequences[seq_id].append(sub_dir.name)

                                elif write_sequence:
                                    out_fh.write(line)

    print(f"Concatenation completed. Results saved to {output_file}")
    # Print skipped sequences (duplicates) for debugging
    if skipped_sequences:
        print("\nDuplicate sequences found in the following subdirectories:")
        for seq_id, sub_dirs in skipped_sequences.items():
            print(f"Sequence ID: {seq_id} appears in directories: {', '.join(sub_dirs)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate enzyme files from KIPEs result folders after filtering by conserved residues.")
    parser.add_argument("-d", "--directory", required=True, help="Path to the KIPEs result folder.")
    parser.add_argument("-e", "--enzymes", required=True, nargs='+', help="List of enzyme names to concatenate.")
    parser.add_argument("-o", "--output", required=True, help="Output file to save concatenated results.")
    parser.add_argument("-t", "--threshold", type=float, default=70.0, help="Minimum percentage of conserved residues (default: 70.0)")
    
    args = parser.parse_args()
    process_kipes_results(args.directory, args.enzymes, args.output, args.threshold)
