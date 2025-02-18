import os
import json
from itertools import product
import subprocess


def load_config_json(config_file_path):
    """
    Load configuration from a JSON file.
    The JSON file should be in the following format:
    {
       "contig_name": ["disallowed_contig1", "disallowed_contig2", ...],
       "another_contig": ["disallowed_contigX", ...]
    }
    """
    with open(config_file_path, "r") as f:
        config_data = json.load(f)

    # Convert all lists to sets for faster membership checks
    config_dict = {k: set(v) for k, v in config_data.items()}
    return config_dict


def index_fasta(fasta_path):
    """
    Index the generated FASTA file using samtools.
    This function requires that samtools is installed and accessible in the system path.
    """
    try:
        subprocess.run(["samtools", "faidx", fasta_path], check=True)
        print(f"Index file '{fasta_path}.fai' created successfully.")
    except FileNotFoundError:
        print(
            "Error: samtools not found. Please ensure that samtools is installed and in your PATH."
        )
    except subprocess.CalledProcessError:
        print("Error: Failed to index the FASTA file with samtools.")


def merge_all_contigs(input_path, output_path, config_path):
    # Load configuration from JSON
    config = load_config_json(config_path)

    with open(input_path, "r") as infile, open(output_path, "w") as outfile:
        lines = infile.readlines()

        # Separate motif IDs and sequences
        headers = [line.strip() for line in lines if line.startswith(">")]
        sequences = [line.strip() for line in lines if not line.startswith(">")]

        # Create a dictionary for easier access
        contigs = {header[1:]: seq for header, seq in zip(headers, sequences)}

        # Generate all combinations, including self-combinations
        for header1, header2 in product(contigs.keys(), repeat=2):
            # Check if there are any disallowed combinations from config
            if header1 in config and header2 in config[header1]:
                # This combination is disallowed based on config
                continue

            # Create the new header
            new_header = f">{header1}-{header2}"

            # Construct the new sequence as per original logic
            new_sequence = contigs[header2] + contigs[header1]

            outfile.write(f"{new_header}\n{new_sequence}\n")

    # After writing the FASTA file, index it using samtools faidx
    index_fasta(output_path)


# Example usage:
if __name__ == "__main__":
    input_file_path = "MUC1_motifs_Rev_com.fa"
    output_file_path = "All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa"
    config_file_path = "filter_config.json"

    merge_all_contigs(input_file_path, output_file_path, config_file_path)
