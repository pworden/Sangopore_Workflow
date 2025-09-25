#!/usr/bin/env python3

"""
conda activate Sangopore_Prep_and_Run

# Make "Sangopore_Workflow" (Containing folder for github package) the working directory
cd /path/to/Sangopore_Workflow

sangopore-prep \
    --input_tsv "Sangopore_Prep_and_Run/Test_Barcoded_Fastq_Inputs/Test_sango_prep_table.tsv" \
    --parent_output_folder "Sangopore_Prep_and_Run/Test_Outputs"
"""

import pandas as pd
import os
import argparse
import gzip
import subprocess
import glob
import shutil

class TSVSubsetter:
    def __init__(self):
        self.input_tsv = None
        self.parent_output_folder = None
        self.df = None
        self.unique_values = None
        self.new_folder_paths = []
        self.chopped_file_paths = {}
        self.fastp_file_paths = {}

    def argument_parser(self):
        parser = argparse.ArgumentParser(description='Subset input TSV and create output folders with options.')
        parser.add_argument('--input_tsv', required=True, help='Path to the input TSV file.')
        parser.add_argument('--parent_output_folder', required=True, help='Path to the parent output folder.')
        parser.add_argument('--r_script_path', required=False, help='Path to the R script.')
        args = parser.parse_args()
        
        self.input_tsv = args.input_tsv
        self.parent_output_folder = args.parent_output_folder
        if args.r_script_path:
            self.r_script_path = args.r_script_path
        else:
            # Get directory where the current python script resides
            script_dir = os.path.dirname(os.path.realpath(__file__))
            self.r_script_path = os.path.join(script_dir, "tidy_multiqc_script.R")

    def read_input_tsv(self):
        self.df = pd.read_csv(self.input_tsv, sep='\t', encoding='utf-8-sig')
        self.df.columns = self.df.columns.str.strip()

    def get_unique_values(self):
        self.unique_values = self.df['subgroup_folder'].unique()

    def create_parent_output_folder(self):
        os.makedirs(self.parent_output_folder, exist_ok=True)

    def create_subdirectories(self):
        for value in self.unique_values:
            folder_path = os.path.join(self.parent_output_folder, value)
            os.makedirs(folder_path, exist_ok=True)
            self.new_folder_paths.append(folder_path)
            self.create_subset_tsv(folder_path, value)

    def create_subset_tsv(self, folder_path, value):
        subset_df = self.df[self.df['subgroup_folder'] == value].copy()
        subset_tsv_filename = os.path.join(folder_path, f"{value}_subset.tsv")
        subset_df.to_csv(subset_tsv_filename, sep='\t', index=False)

    def write_new_folder_paths(self):
        list_file_path = os.path.join(self.parent_output_folder, "new_folder_paths.txt")
        with open(list_file_path, 'w') as list_file:
            for path in self.new_folder_paths:
                list_file.write(f"{path}\n")

    def process_subset_tsv_files(self):
        """
        Process only original input TSVs in the provided parent folder.
        Ignores any generated subset TSVs in the output folder.
        """
        # If you supplied the original input TSV, process it directly
        if os.path.exists(self.input_tsv):
            print(f"Processing original TSV: {self.input_tsv}")
            # Each subgroup folder will be created as usual
            subset_folder_name = os.path.splitext(os.path.basename(self.input_tsv))[0]
            folder_path = os.path.join(self.parent_output_folder, subset_folder_name)
            os.makedirs(folder_path, exist_ok=True)
            self.process_tsv(self.input_tsv, folder_path)
        else:
            print(f"Warning: Original TSV {self.input_tsv} not found, skipping.")


    def process_tsv(self, tsv_path, output_base_folder):
        print(f"Reading TSV: {tsv_path}")
        df = pd.read_csv(tsv_path, sep='\t', encoding='utf-8-sig')
        df.columns = df.columns.str.strip()
        print("Columns found:", df.columns.tolist())
        fastp_dirs = set()
        all_files_in_one_folder = df['input_barcode_file'].nunique() == 1
        common_folder = df['input_barcode_file'].iloc[0] if all_files_in_one_folder else None
        
        for index, row in df.iterrows():
            barcode = row['barcode_number']
            subgroup_folder = row['subgroup_folder']
            search_path = row['input_barcode_file']
            output_subfolder = os.path.join(output_base_folder, subgroup_folder)
            os.makedirs(output_subfolder, exist_ok=True)

            cat_file = os.path.join(output_subfolder, f"{barcode}_{subgroup_folder}_cat.fastq")
            chopped_file = os.path.join(output_subfolder, f"{barcode}_{subgroup_folder}_chopped.fastq")
            fastp_json_file = os.path.join(output_subfolder, f"{barcode}_{subgroup_folder}_fastp.json")
            fastp_html_file = os.path.join(output_subfolder, f"{barcode}_{subgroup_folder}_fastp.html")
            fastp_output_file = os.path.join(output_subfolder, f"{barcode}_{subgroup_folder}_fastp.fastq.gz")
            
            if all_files_in_one_folder:
                matching_files = [
                    file for file in os.listdir(common_folder)
                    if barcode in file and file.endswith('.fastq.gz')
                ]
                matching_files = [os.path.join(common_folder, file) for file in matching_files]
            else:
                matching_files = [
                    os.path.join(root, file)
                    for root, _, files in os.walk(search_path)
                    for file in files
                    if barcode in file and file.endswith('.fastq.gz')
                ]

            if not matching_files:
                continue
            
            with open(cat_file, 'wb') as cat_f:
                for file in matching_files:
                    with gzip.open(file, 'rb') as f:
                        cat_f.write(f.read())
            
            try:
                chopper_command = f"chopper -q 10 -i {cat_file} > {chopped_file}"
                subprocess.run(chopper_command, shell=True, check=True, universal_newlines=True)
                self.chopped_file_paths[(barcode, subgroup_folder)] = chopped_file
            except subprocess.CalledProcessError as e:
                print(f"Error running chopper for {barcode}: {e}")
            
            try:
                fastp_command = (
                    f"fastp -i {chopped_file} "
                    f"-j {fastp_json_file} "
                    f"-h {fastp_html_file} "
                    f"-o {fastp_output_file} "
                )
                subprocess.run(fastp_command, shell=True, check=True, universal_newlines=True)
                fastp_dirs.add(output_subfolder)
                self.fastp_file_paths[(barcode, subgroup_folder)] = fastp_output_file
            except subprocess.CalledProcessError as e:
                print(f"Error running fastp for {barcode}: {e}")

            os.remove(cat_file)
            os.remove(chopped_file)
        
        try:
            multiqc_command = f"multiqc {' '.join(fastp_dirs)} -o {output_base_folder}"
            subprocess.run(multiqc_command, shell=True, check=True, universal_newlines=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running multiqc: {e}")

    def run_multiqc_on_parent(self):
        try:
            multiqc_command = f"multiqc {self.parent_output_folder} -o {self.parent_output_folder}"
            subprocess.run(multiqc_command, shell=True, check=True, universal_newlines=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running multiqc on parent folder: {e}")

    def run_r_script_on_final_json(self):
        final_json_path = os.path.join(self.parent_output_folder, 'multiqc_data', 'multiqc_data.json')
        if not self.r_script_path:
            raise ValueError("R script path is not set.")
        command = [
            'Rscript',
            self.r_script_path,
            final_json_path,
            os.path.join(self.parent_output_folder, 'summary_multiqc_table.csv')
        ]
        try:
            subprocess.run(command, check=True, universal_newlines=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running R script: {e}")

    def create_final_tsv(self):
        df = self.df.copy()
        df.rename(columns={"input_barcode_file": "fastq_file_paths"}, inplace=True)
        for index, row in df.iterrows():
            key = (row['barcode_number'], row['subgroup_folder'])
            if key in self.fastp_file_paths:
                df.at[index, "fastq_file_paths"] = self.fastp_file_paths[key]
            else:
                df.at[index, "fastq_file_paths"] = "Path not found"
        
        # Get the parent directory name
        parent_dir_name = os.path.basename(self.parent_output_folder)
        
        # Add the parent directory name as a prefix to the TSV filename
        final_tsv_filename = f"{parent_dir_name}_sangopore_prep_qc.tsv"
        final_tsv_path = os.path.join(self.parent_output_folder, final_tsv_filename)
        
        df.to_csv(final_tsv_path, sep='\t', index=False)
        
        for value in df['subgroup_folder'].unique():
            subset_df = df[df['subgroup_folder'] == value].copy()
            folder_path = os.path.join(self.parent_output_folder, value)
            os.makedirs(folder_path, exist_ok=True)
            direct_path_tsv_filename = os.path.join(folder_path, f"{value}_directPath_subset.tsv")
            subset_df.to_csv(direct_path_tsv_filename, sep='\t', index=False)

    def run_pipeline(self):
        self.argument_parser()
        self.read_input_tsv()
        self.get_unique_values()
        self.create_parent_output_folder()
        self.create_subdirectories()
        self.write_new_folder_paths()
        self.process_subset_tsv_files()
        self.run_multiqc_on_parent()
        self.create_final_tsv()
        self.run_r_script_on_final_json()
        print("Pipeline completed.")

def main():
    subsetter = TSVSubsetter()
    subsetter.run_pipeline()

if __name__ == "__main__":
    main()
