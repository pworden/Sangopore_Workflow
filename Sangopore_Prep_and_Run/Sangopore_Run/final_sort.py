#!/usr/bin/env python3

"""
GitHub Copilot
# The current script organizes and collates analysis results from a Sangopore run by creating barcode-specific subfolders under a results_per_barcode directory.
# It then copies relevant FASTQ, PDF, FASTA, Excel, and BLASTN output files and folders into their corresponding barcode folders.
# This script can be used alone or it is called by the Sangopore_ampliconSorter.py script through "import".
    # For example, in the Sangopore_ampliconSorter.py script: "from Sangopore_Prep_and_Run_Dev.Sangopore_Run import final_sort"


conda activate Sangopore_Prep_and_Run_Dev
    
python sangopore-results-sort \
    --parent-input-output-dir "/home/wordenp/Test_Sango_Out_19-6-25/sangoAnalysis_out/All_M24-13575"
"""

import re
from pathlib import Path
import shutil

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Collate compressed FASTQ files into results_per_barcode structure")
    parser.add_argument("--parent-input-output-dir", required=True, help="Path to parent directory of Sangopore results, containing many subfolders of output")
    return parser.parse_args()

class ResultsCollation:
    def __init__(self, parent_input_output_dir: Path):
        self.parent_input_output_dir = parent_input_output_dir.resolve()
        self.compressed_path = self.parent_input_output_dir / "compressed"
        self.barcode_ids = []  # Populated after get_created_barcode_folders()
        self.results_barcode_path = self.parent_input_output_dir / "results_per_barcode"

    def create_barcode_subfolders(self):
        pattern = "*barcode[0-9][0-9]*.fastq.gz"
        for fastq_file in self.compressed_path.glob(pattern):
            barcode_match = re.search(r"barcode\d{2}", fastq_file.name)
            if barcode_match:
                barcode_id = barcode_match.group(0)
                target_dir = self.results_barcode_path / barcode_id
                target_dir.mkdir(parents=True, exist_ok=True)
                print(f"[DEBUG] Created directory: {target_dir}")
            else:
                print(f"[WARNING] No barcode found in filename: {fastq_file.name}")

    def copy_compressed_fastq_files(self):
        pattern = "*barcode[0-9][0-9]*.fastq.gz"
        for fastq_file in self.compressed_path.glob(pattern):
            barcode_match = re.search(r"barcode\d{2}", fastq_file.name)
            if barcode_match:
                barcode_id = barcode_match.group(0)
                target_path = self.results_barcode_path / barcode_id / fastq_file.name
                shutil.copy2(fastq_file, target_path)
                print(f"[INFO] Copied: {fastq_file.name} → {target_path}")
            else:
                print(f"[WARNING] Skipping file (no barcode found): {fastq_file.name}")

    def get_created_barcode_folders(self):
        """
        Populate self.barcode_ids based on folders under results_per_barcode.
        """
        if not self.results_barcode_path.exists():
            print(f"[WARNING] results_per_barcode path does not exist: {self.results_barcode_path}")
            return []

        self.barcode_ids = [
            folder.name
            for folder in self.results_barcode_path.iterdir()
            if folder.is_dir() and re.match(r"barcode\d{2}", folder.name)
        ]

        print(f"[INFO] Found barcode folders: {self.barcode_ids}")
        return self.barcode_ids

    def copy_qc_pdfs_to_barcode_folders(self):
        """
        For each barcode folder, find a matching PDF in Amplicon_Sorter_QC and copy it into the corresponding
        results_per_barcode/barcodeXX/ folder.
        """
        qc_path = self.parent_input_output_dir / "Amplicon_Sorter_QC"
        if not qc_path.exists():
            print(f"[WARNING] No Amplicon_Sorter_QC directory found in: {self.parent_input_output_dir}")
            return

        for barcode in self.barcode_ids:
            # Match any subfolder starting with barcodeXX, e.g. barcode01, barcode01_fp, barcode01_xyz
            subfolder_pattern = f"{barcode}*"
            subfolders = list(qc_path.glob(subfolder_pattern))

            if not subfolders:
                print(f"[WARNING] No QC subfolder found for {barcode}")
                continue

            for subfolder in subfolders:
                # Look for a PDF file inside the subfolder like '*barcode01*.pdf'
                pdf_pattern = f"*{barcode}*.pdf"
                matching_pdfs = list(subfolder.glob(pdf_pattern))

                if not matching_pdfs:
                    print(f"[WARNING] No PDF found in {subfolder} matching {pdf_pattern}")
                    continue

                for pdf_file in matching_pdfs:
                    dest = self.results_barcode_path / barcode / pdf_file.name
                    shutil.copy2(pdf_file, dest)
                    print(f"[INFO] Copied QC PDF for {barcode}: {pdf_file.name} → {dest}")

    def copy_consensus_fastas_to_barcode_folders(self):
        """
        For each barcode folder, find matching FASTA files in amplicon_sorter_consensus_multiFasta
        and copy them into the corresponding results_per_barcode/barcodeXX/ folder.
        """
        consensus_path = self.parent_input_output_dir / "amplicon_sorter_consensus_multiFasta"
        if not consensus_path.exists():
            print(f"[WARNING] No amplicon_sorter_consensus_multiFasta directory found in: {self.parent_input_output_dir}")
            return

        for barcode in self.barcode_ids:
            # Find all FASTA files matching *barcodeXX*.fasta
            fasta_pattern = f"*{barcode}*.fasta"
            matching_fastas = list(consensus_path.glob(fasta_pattern))

            if not matching_fastas:
                print(f"[WARNING] No consensus FASTA files found for {barcode} matching {fasta_pattern}")
                continue

            for fasta_file in matching_fastas:
                dest = self.results_barcode_path / barcode / fasta_file.name
                shutil.copy2(fasta_file, dest)
                print(f"[INFO] Copied consensus FASTA for {barcode}: {fasta_file.name} → {dest}")

    def copy_blastn_excel_to_barcode_folders(self):
        """
        Copy BLASTN Excel files matching each barcode into results_per_barcode/barcodeXX/ folders.
        """
        blast_excel_path = self.parent_input_output_dir / "amplicon_sorter_blastn_20_ordered__excel"
        if not blast_excel_path.exists():
            print(f"[WARNING] No amplicon_sorter_blastn_20_ordered__excel directory found in: {self.parent_input_output_dir}")
            return

        for barcode in self.barcode_ids:
            # Pattern example: blastn_collation_barcode01_ordered20.xlsx
            excel_pattern = f"blastn_collation_{barcode}_*.xlsx"
            matching_excels = list(blast_excel_path.glob(excel_pattern))

            if not matching_excels:
                print(f"[WARNING] No BLASTN Excel files found for {barcode} matching {excel_pattern}")
                continue

            for excel_file in matching_excels:
                dest = self.results_barcode_path / barcode / excel_file.name
                shutil.copy2(excel_file, dest)
                print(f"[INFO] Copied BLASTN Excel for {barcode}: {excel_file.name} → {dest}")

    def copy_consensus_singleSeq_fastas_to_barcode_folders(self):
        """
        For each barcode, copy all matching single consensus FASTA files into
        a dedicated subfolder 'singleSeq_fasta_barcodeXX' inside results_per_barcode/barcodeXX/.
        """
        single_fasta_path = self.parent_input_output_dir / "amplicon_sorter_consensus_singleFasta"
        if not single_fasta_path.exists():
            print(f"[WARNING] No amplicon_sorter_consensus_singleFasta directory found in: {self.parent_input_output_dir}")
            return

        for barcode in self.barcode_ids:
            # Pattern example: barcode01_fp_0_0__352_consensus.fasta
            fasta_pattern = f"{barcode}*.fasta"
            matching_fastas = list(single_fasta_path.glob(fasta_pattern))

            if not matching_fastas:
                print(f"[WARNING] No single consensus FASTA files found for {barcode} matching {fasta_pattern}")
                continue

            dest_dir = self.results_barcode_path / barcode / f"singleSeq_fasta_{barcode}"
            dest_dir.mkdir(parents=True, exist_ok=True)

            for fasta_file in matching_fastas:
                dest = dest_dir / fasta_file.name
                shutil.copy2(fasta_file, dest)
                print(f"[INFO] Copied single consensus FASTA for {barcode}: {fasta_file.name} → {dest}")

    def copy_blastn_out_folders_to_barcode_folders(self):
        """
        Copy entire blastn result folders (e.g. blastn_barcode01) from
        amplicon_sorter_blastn_out into results_per_barcode/barcodeXX/
        """
        blast_out_path = self.parent_input_output_dir / "amplicon_sorter_blastn_out"
        if not blast_out_path.exists():
            print(f"[WARNING] No amplicon_sorter_blastn_out directory found in: {self.parent_input_output_dir}")
            return

        for barcode in self.barcode_ids:
            folder_name = f"blastn_{barcode}"
            source_folder = blast_out_path / folder_name

            if not source_folder.exists() or not source_folder.is_dir():
                print(f"[WARNING] No blastn output folder found for {barcode} at {source_folder}")
                continue

            dest_folder = self.results_barcode_path / barcode / folder_name

            if dest_folder.exists():
                print(f"[INFO] Destination folder already exists, skipping copy: {dest_folder}")
                continue

            # Use shutil.copytree to copy entire folder
            shutil.copytree(source_folder, dest_folder)
            print(f"[INFO] Copied BLASTN folder for {barcode}: {source_folder} → {dest_folder}")

    def run_pipeline(self):
        if not self.compressed_path.exists():
            print(f"[ERROR] 'compressed/' folder not found in: {self.parent_input_output_dir}")
            return

        self.results_barcode_path.mkdir(parents=True, exist_ok=True)
        print(f"[INFO] Created results root: {self.results_barcode_path}")

        self.create_barcode_subfolders()
        self.copy_compressed_fastq_files()
        self.get_created_barcode_folders()
        self.copy_qc_pdfs_to_barcode_folders()
        self.copy_consensus_fastas_to_barcode_folders()
        # self.copy_consensus_singleSeq_fastas_to_barcode_folders()
        self.copy_blastn_excel_to_barcode_folders()
        self.copy_blastn_out_folders_to_barcode_folders()

def main():
    args = parse_args()
    collator = ResultsCollation(Path(args.parent_input_output_dir))
    collator.run_pipeline()

if __name__ == "__main__":
    main()
