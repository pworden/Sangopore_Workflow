#!/usr/bin/env python3

"""
# This script inputs a CSV file manually captured from the Excel sample booking sheet (Sangopore) that includes both barcodes and sample identifiers.
# The ONT barcodes and sample identifiers in the input CSV are concatenated and used to rename files and directories under a user defined parent directory.
# A CSV of identifiers and barcodes, as well as their concatenated products (strings) is also output. 

conda activate Sangopore_Prep_and_Run_Dev

sangopore-expand-bc-names \
    --input-barcode-and-id-csv "Sangopore_Prep_and_Run_Dev/Test_Barcoded_Fastq_Inputs/tutorial.csv" \
    --data-parent-dir "Sangopore_Prep_and_Run_Dev/Test_Outputs/sangoAnalysis_out/Group1/results_per_barcode"
    
    # Optionally add for 'cleaned_sample_names.csv' output: --output-csv "/some/other/path/cleaned_sample_names.csv"


conda activate Sangopore_Prep_and_Run_Dev

python "/home/wordenp/Scripts/General_Scripts/Sangopore_Workflow/Sangopore_Prep_and_Run_Dev/Sangopore_Run/expand_barcodes_with_sampleIDs.py" \
    --input-barcode-and-id-csv "/home/wordenp/projects/sangopore_project/analyses/Sangopore_analysis/SP25-006_Sangopore_17JUN25_rerun/backup_17-6-2025_Sangopore.csv" \
    --data-parent-dir "/home/wordenp/projects/sangopore_project/analyses/Sangopore_analysis/SP25-006_Sangopore_17JUN25_rerun/sangoAnalysis_out/Biotech/results_per_barcode2"
    
    # Optionally add for 'cleaned_sample_names.csv' output: --output-csv "/some/other/path/cleaned_sample_names.csv"
"""

import os
import re
import csv
import argparse
from typing import List, Dict, Union

def parse_args():
    parser = argparse.ArgumentParser(description="Sangopore workflow: extract and expand barcode info.")
    parser.add_argument('--input-barcode-and-id-csv', required=True, help='Input CSV file from Excel (for extraction step)')
    parser.add_argument('--data-parent-dir', required=True, help='Parent directory to rename files/folders in')
    parser.add_argument('--output-csv', help='Output CSV file for cleaned sample names (default: cleaned_sample_names.csv in data-parent-dir)')
    return parser.parse_args()

def clean_sample_name(name: str) -> str:
    name = re.sub(r'[^A-Za-z0-9]+', '_', name)
    name = re.sub(r'_+', '_', name)
    return name.strip('_')

def strip_metadata_line(lines: List[str]) -> List[str]:
    for i, line in enumerate(lines):
        for delimiter in [',', '\t']:
            fields = [field.strip() for field in line.split(delimiter)]
            headers_lower = [f.lower() for f in fields]
            if 'sample name' in headers_lower and 'ont barcode' in headers_lower:
                return lines[i:]
    raise ValueError("Could not find valid header row with 'Sample Name' and 'ONT Barcode'")

def format_barcode(value: str) -> str:
    try:
        num = int(value)
        return f"barcode{num:02d}"
    except (ValueError, TypeError):
        return "barcode??"

class SampleBarcodeExtractor:
    def __init__(self, input_barcode_and_id_csv: str, output_csv: str):
        self.input_barcode_and_id_csv = input_barcode_and_id_csv
        self.output_csv = output_csv

    def extract(self) -> List[Dict[str, str]]:
        with open(self.input_barcode_and_id_csv, "r", encoding="utf-8") as f:
            lines = f.read().strip().splitlines()
        lines = strip_metadata_line(lines)
        text = "\n".join(lines)
        sniffer = csv.Sniffer()
        try:
            dialect = sniffer.sniff(text, delimiters=[',', '\t'])
        except csv.Error:
            dialect = csv.get_dialect('excel')
        reader = csv.DictReader(lines, dialect=dialect)
        results = []
        for row in reader:
            original_name = row.get("Sample Name", "").strip()
            barcode_val = row.get("ONT Barcode", "").strip()
            if original_name:
                cleaned_name = clean_sample_name(original_name)
                formatted_barcode = format_barcode(barcode_val)
                id_and_barcode = f"{formatted_barcode}_{cleaned_name}"
                results.append({
                    "Cleaned_Sample_Name": cleaned_name,
                    "ONT_Barcode": formatted_barcode,
                    "id_and_barcode": id_and_barcode
                })
        return results

    def write_output(self, rows: List[Dict[str, str]]):
        if not rows:
            print("No valid rows to write.")
            return
        fieldnames = ["Cleaned_Sample_Name", "ONT_Barcode", "id_and_barcode"]
        with open(self.output_csv, "w", newline='', encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print(f"Wrote output to: {self.output_csv}")

    def run(self):
        results = self.extract()
        for row in results:
            print(f"{row['Cleaned_Sample_Name']}\t{row['ONT_Barcode']}\t{row['id_and_barcode']}")
        self.write_output(results)

class BarcodeExpander:
    def __init__(self, mapping_csv: str, parent_dir: str):
        self.mapping_csv = mapping_csv
        self.parent_dir = parent_dir

    def read_mapping(self) -> Dict[str, str]:
        barcode_to_sample = {}
        with open(self.mapping_csv, mode='r', encoding='utf-8') as file:
            reader = csv.DictReader(file, delimiter=',')
            for row in reader:
                barcode_to_sample[row['ONT_Barcode'].strip()] = row['id_and_barcode'].strip()
        return barcode_to_sample

    def rename(self, barcode_to_sample: Dict[str, str]):
        for root, dirs, files in os.walk(self.parent_dir, topdown=False):
            for name in files:
                for barcode, sample in barcode_to_sample.items():
                    if barcode in name:
                        old_path = os.path.join(root, name)
                        new_name = name.replace(barcode, sample)
                        os.rename(old_path, os.path.join(root, new_name))
            for name in dirs:
                for barcode, sample in barcode_to_sample.items():
                    if barcode in name:
                        old_path = os.path.join(root, name)
                        new_name = name.replace(barcode, sample)
                        os.rename(old_path, os.path.join(root, new_name))

    def run(self):
        barcode_to_sample = self.read_mapping()
        self.rename(barcode_to_sample)
        print("Renaming completed successfully.")

def pipeline():
    args = parse_args()

    # Determine output CSV path
    if args.output_csv:
        output_csv = args.output_csv
    else:
        output_csv = os.path.join(args.data_parent_dir, "cleaned_sample_names.csv")

    # Extraction step
    extractor = SampleBarcodeExtractor(args.input_barcode_and_id_csv, output_csv)
    extractor.run()

    # Expansion step (always uses the output CSV from extraction)
    expander = BarcodeExpander(output_csv, args.data_parent_dir)
    expander.run()

def main():
    pipeline()

if __name__ == "__main__":
    main()
