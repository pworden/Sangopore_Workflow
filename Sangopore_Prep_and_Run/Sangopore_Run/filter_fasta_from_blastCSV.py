#!/usr/bin/env python3

"""
This script filters multi-amplicon FASTA files based on sequence abundance from corresponding BLASTN result CSVs.
For each barcode, it retains amplicon sequences that meet or exceed a specified abundance threshold and writes the filtered FASTA to the output directory.

conda activate Sangopore_Prep_and_Run_Dev

python "/home/wordenp/Scripts/General_Scripts/Sangopore_Workflow/Sangopore_Prep_and_Run_Dev/Sangopore_Run/filter_fasta_from_blastCSV.py" \
    --csv_dir "/home/wordenp/projects/sangopore_project/analyses/Sangopore_analysis/SP25-006_Sangopore_17JUN25_max10000/sangoAnalysis_out/Biotech/amplicon_sorter_blastn_20_ordered_collation" \
    --fasta_dir "/home/wordenp/projects/sangopore_project/analyses/Sangopore_analysis/SP25-006_Sangopore_17JUN25_max10000/sangoAnalysis_out/Biotech/amplicon_sorter_consensus_multiFasta" \
    --output_dir "/home/wordenp/projects/sangopore_project/analyses/Sangopore_analysis/SP25-006_Sangopore_17JUN25_max10000/sangoAnalysis_out/Biotech/amplicon_sorter_consensus_multiFasta_abund_filtered" \
    --threshold 1.0

# ----------

conda activate Sangopore_Prep_and_Run_Dev

sango-filter-fasta-with-blastCSV \
    --csv_dir "/home/wordenp/projects/sangopore_project/analyses/Sangopore_analysis/SP25-006_Sangopore_17JUN25_max10000/sangoAnalysis_out/Biotech/amplicon_sorter_blastn_20_ordered_collation" \
    --fasta_dir "/home/wordenp/projects/sangopore_project/analyses/Sangopore_analysis/SP25-006_Sangopore_17JUN25_max10000/sangoAnalysis_out/Biotech/amplicon_sorter_consensus_multiFasta" \
    --output_dir "/home/wordenp/projects/sangopore_project/analyses/Sangopore_analysis/SP25-006_Sangopore_17JUN25_max10000/sangoAnalysis_out/Biotech/amplicon_sorter_consensus_multiFasta_abund_filtered" \
    --threshold 1.0
"""

import os, re
import argparse
from pathlib import Path
import pandas as pd
from Bio import SeqIO

def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter FASTA sequences based on CSV 'Percent_Amp_reads_of_analysis_inputs_per_BC' thresholds.")
    parser.add_argument("--csv_dir", required=True, help="Path to directory containing CSV files like 'blastn_collation_barcode*_ordered20.csv'")
    parser.add_argument("--fasta_dir", required=True, help="Path to directory containing multi-amplicon FASTA files per barcode")
    parser.add_argument("--output_dir", required=True, help="Directory to write filtered FASTA files to (e.g., 'filtered_amp_sorter_consensus_multiFasta')")
    parser.add_argument("--threshold", type=float, default=1.0, help="Minimum Percent_Amp_reads_of_analysis_inputs_per_BC required to keep an amplicon (default: 1.0)")
    return parser.parse_args()

def get_passing_query_ids(csv_path, threshold):
    df = pd.read_csv(csv_path)
    if "Query_Seq_id" not in df.columns or "Percent_Amp_reads_of_analysis_inputs_per_BC" not in df.columns:
        raise ValueError(f"Missing required columns in {csv_path}")
    return set(df.loc[df["Percent_Amp_reads_of_analysis_inputs_per_BC"] >= threshold, "Query_Seq_id"])

def filter_fasta(input_fasta, output_fasta, allowed_ids):
    sequences = list(SeqIO.parse(input_fasta, "fasta"))

    # DEBUG: show a few original and cleaned FASTA IDs
    original_ids = [rec.id for rec in sequences[:3]]
    cleaned_ids = [re.sub(r"\(\d+\)$", "", rec.id) for rec in sequences[:3]]

    print(f"[DEBUG] Sample FASTA headers: {original_ids}")
    print(f"[DEBUG] Cleaned FASTA IDs: {cleaned_ids}")
    print(f"[DEBUG] Allowed IDs (from CSV): {list(allowed_ids)[:3]}")

    # Filter
    kept = [rec for rec in sequences if re.sub(r"\(\d+\)$", "", rec.id) in allowed_ids]
    SeqIO.write(kept, output_fasta, "fasta")
    return len(kept), len(sequences)

def pipeline(csv_dir, fasta_dir, output_dir, threshold):
    csv_dir = Path(csv_dir)
    fasta_dir = Path(fasta_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    csv_files = sorted(csv_dir.glob("blastn_collation_barcode*_ordered20.csv"))

    for csv_file in csv_files:
        try:
            barcode = csv_file.stem.split("barcode")[1].split("_")[0]
        except IndexError:
            print(f"[WARNING] Skipping malformed CSV name: {csv_file.name}")
            continue

        matching_fastas = [f for f in fasta_dir.glob("*.fasta") if f"barcode{barcode}" in f.name]
        if not matching_fastas:
            print(f"[WARNING] No FASTA file found containing 'barcode{barcode}' in name")
            continue
        elif len(matching_fastas) > 1:
            print(f"[WARNING] Multiple FASTA files found for barcode{barcode}, using first: {[f.name for f in matching_fastas]}")
        fasta_file = matching_fastas[0]

        try:
            allowed_ids = get_passing_query_ids(csv_file, threshold)
            output_fasta = output_dir / f"barcode{barcode}_amplicons_filtered.fasta"
            kept, total = filter_fasta(fasta_file, output_fasta, allowed_ids)
            print(f"[INFO] Barcode {barcode}: retained {kept}/{total} sequences")
        except Exception as e:
            print(f"[ERROR] Failed processing barcode {barcode}: {e}")

def main():
    args = parse_arguments()
    pipeline(
        csv_dir=args.csv_dir,
        fasta_dir=args.fasta_dir,
        output_dir=args.output_dir,
        threshold=args.threshold
    )

if __name__ == "__main__":
    main()
