# Reference folder for vntyper

## 1. VNTR Reference Generator

This repository contains a Python script `generate_vntr_reference.py` which generates pairwise combinations of VNTR motifs from a given input FASTA file and writes them out to a new FASTA file. The script also uses a JSON configuration file to filter out certain disallowed motif combinations.

## Overview

The script:
1. Reads an input FASTA file that contains multiple contigs (motifs).
2. Generates all pairwise combinations of these contigs, including self-combinations.
3. Concatenates the sequences in a specified order to form the new reference entries.
4. Uses a JSON configuration file to skip (filter out) certain disallowed combinations.

## Files

- **MUC1_motifs_Rev_com.fa**: The input FASTA file containing the original VNTR motifs.  
  Each contig is represented as:
  ```
  >contigName
  ACTGACTG...
  ```

- **filter_config.json**: The JSON configuration file that specifies disallowed combinations.  
  Example format:
  ```json
  {
    "1": ["1", "3", "4", ...],
    "2": ["1", "2", "4", ...],
    ...
  }
  ```
  This means if the first contig is `"1"`, it should not be paired with `"1"`, `"3"`, `"4"`, etc.

- **generate_vntr_reference.py**: The Python script that:
  - Loads the configuration JSON file.
  - Iterates over all contigs from the input FASTA file.
  - Generates all pairwise combinations.
  - Skips any combinations found in the JSON configuration file.
  - Writes the allowed combinations to the output FASTA file.

## Prerequisites

- Python 3.x
- A JSON configuration file with the appropriate format.
- A FASTA file with your source motifs.

## Usage

1. **Prepare Your Input Files**:
   - Ensure `MUC1_motifs_Rev_com.fa` (or your chosen input FASTA) is in the same directory as the script.
   - Place `filter_config.json` in the same directory and confirm it has the correct format.

2. **Run the Script**:
   From the terminal or command prompt:
   ```bash
   python3 generate_vntr_reference.py
   ```
   
   By default, the script is set to:
   - Read from `MUC1_motifs_Rev_com.fa`.
   - Write to `All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa`.
   - Use `filter_config.json` for filtering.

   If you want to customize paths, open `generate_vntr_reference.py` and modify:
   ```python
   input_file_path = 'MUC1_motifs_Rev_com.fa'
   output_file_path = 'All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa'
   config_file_path = 'filter_config.json'
   ```
   Then run the script again.

3. **Check the Output**:
   After running the script, you will have a new FASTA file (`All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa`) containing only the allowed motif combinations.

## Configuration Details

- The `filter_config.json` file maps each "first contig" to a list of "second contigs" that cannot follow it.
- If a combination `first-second` is disallowed, it will be skipped entirely.
- Modify `filter_config.json` as needed to add or remove disallowed combinations. Any combination where the first contig is not listed or the second contig is not in the corresponding array/set will be allowed.

## Example

If you have:
```json
{
  "1": ["1", "3"]
}
```

This means:
- When the first contig is `1`, do not allow `1-1` or `1-3`.
- All other combinations starting with `1` (except `1-1` and `1-3`) are allowed.

## 2. pseudonymize.py

## 3. Adapted adVNTR references

The referneces in vntr_db_advntr.zip have been created by removing all non-Muc1 entries from the hg9 database and adding the MUC1 reference with the same ID to the hgg38 database where it was missing and changing the start position to the corresponding start position after finding this with UCSC BLAT usng the left sequence.