#!/usr/bin/env python3

import pandas as pd
import os
import re
import argparse

"""
conda activate Sangopore_Prep_and_Run

sango-create-table-onedir-allBCs \
    --parent_input_dir "/home/wordenp/AGTC_EMAI_NAS/SRA/GridION/Biotechnology/SP25_015_Sangopore_28thOCT25/20251028_SP25_015/no_sample_id/20251028_1314_X2_BAJ991_de9556df/fastq_hac_pass" \
    --output_dir "/home/wordenp/projects/sangopore_project/analyses/Sangopore_analysis/SP26-016_Sangopore_31stJUL2026" \
    --output_tsv_name "AllInOneDir_table.tsv"
"""

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description='Process .fastq files to create a TSV file with barcode information.')
    parser.add_argument('--parent_input_dir', type=str, required=True, help='Path to the parent input directory containing .fastq files, one fastq file per barcode.')
    parser.add_argument('--output_dir', type=str, required=True, help='Path to the output directory where the TSV file will be saved.')
    parser.add_argument('--output_tsv_name', type=str, required=True, help='User name of the output TSV file.')

    args = parser.parse_args()

    # User-defined inputs from command-line arguments
    parentInputDir = args.parent_input_dir
    outputDir = args.output_dir
    outTableName = args.output_tsv_name
    outputTablePath = os.path.join(outputDir, outTableName)

    # Create the output directory if it doesn't exist
    os.makedirs(outputDir, exist_ok=True)

    # Initialize a list to hold the rows of data
    data = []

    print(f"Searching for .fastq.gz files in {parentInputDir}...")

    # Traverse the directory to find all .fastq files
    for root, _, files in os.walk(parentInputDir):
        for file in files:
            if file.endswith(".fastq.gz"):
                # Extract barcode number from filename
                match = re.search(r'barcode([0-9]{2})', file)
                if match:
                    barcode_number = f"barcode{match.group(1)}"
                    
                    # Extract folder path from the file path (full path to the directory containing the .fastq.gz file)
                    folder_path = os.path.dirname(os.path.join(root, file))
                    
                    # Append the data
                    data.append([barcode_number, folder_path, "", "", ""])

    print(f"Found {len(data)} entries.")

    # Create a DataFrame from the collected data
    df = pd.DataFrame(data, columns=["barcode_number", "input_barcode_file", "subgroup_folder", "min_range_bp", "max_range_bp"])

    # Check if the DataFrame is empty
    if df.empty:
        print("No data to process.")
    else:
        # Remove duplicates based on both columns 1 and 2
        df = df.drop_duplicates(subset=["barcode_number", "input_barcode_file"])
        # Extract the numeric part of the barcode_number for sorting
        df['barcode_number_numeric'] = df['barcode_number'].str.extract(r'barcode([0-9]{2})').astype(int)
        # Sort the DataFrame based on the numeric part of barcode_number
        df = df.sort_values(by='barcode_number_numeric')
        # Drop the temporary numeric column used for sorting
        df = df.drop(columns=['barcode_number_numeric'])
        # Save the DataFrame to a TSV file
        df.to_csv(outputTablePath, sep='\t', index=False)
        print(f"TSV file saved to {outputTablePath}.")

if __name__ == "__main__":
    main()
