#!/usr/bin/env python3

"""
conda activate Sangopore_Prep_and_Run

sango-create-table-multidir-multifiles \
    --parentInputDir "Sangopore_Prep_and_Run/Test_Barcoded_Fastq_Inputs" \
    --outputDir "Sangopore_Prep_and_Run/Test_Outputs" \
    --outTableName "MultiDirMultiFiles_table.tsv"
"""

#!/usr/bin/env python3

import os
import re
import pandas as pd
import argparse

class FASTQProcessor:
    def __init__(self, parentInputDir, outputDir, outTableName):
        self.parentInputDir = parentInputDir
        self.outputDir = outputDir
        self.outTableName = outTableName
        self.outputTablePath = os.path.join(outputDir, outTableName)
        self.data = []

    def traverse_directory(self):
        barcode_regex = re.compile(r'barcode(\d{2})')

        for root, dirs, files in os.walk(self.parentInputDir):
            for file in files:
                if file.endswith(".fastq.gz"):
                    self.process_file(file, root, barcode_regex)

    def process_file(self, file, root, barcode_regex):
        file_path = os.path.join(root, file)
        match = barcode_regex.search(file)
        if match:
            barcode_number = f"barcode{match.group(1)}"
            directory = os.path.dirname(file_path)
            dir_match = barcode_regex.search(os.path.basename(directory))
            if dir_match:
                dir_barcode = f"barcode{dir_match.group(1)}"
                if barcode_number == dir_barcode:
                    self.data.append([barcode_number, directory, "", "", ""])

    def create_dataframe(self):
        df = pd.DataFrame(self.data, columns=["barcode_number", "input_barcode_file", "subgroup_folder", "min_range_bp", "max_range_bp"])
        df = df.drop_duplicates()
        df['barcode_number_numeric'] = df['barcode_number'].str.extract(r'(\d+)').astype(float)
        df_sorted = df.sort_values(by='barcode_number_numeric').drop(columns='barcode_number_numeric')
        return df_sorted

    def save_table(self, df):
        os.makedirs(self.outputDir, exist_ok=True)
        df.to_csv(self.outputTablePath, sep='\t', index=False, header=True)
        print(f"Table has been written to {self.outputTablePath}")

    def run_pipeline(self):
        self.traverse_directory()
        df = self.create_dataframe()
        self.save_table(df)

# This is the main function that will be linked in pyproject.toml
def main():
    parser = argparse.ArgumentParser(description='Process FASTQ files and generate a TSV file with barcode information.')
    parser.add_argument('--parentInputDir', required=True, help='Path to the parent input directory containing FASTQ files.')
    parser.add_argument('--outputDir', required=True, help='Path to the output directory where the TSV file will be saved.')
    parser.add_argument('--outTableName', required=True, help='Name of the output TSV file.')
    args = parser.parse_args()

    processor = FASTQProcessor(args.parentInputDir, args.outputDir, args.outTableName)
    processor.run_pipeline()

if __name__ == "__main__":
    main()
