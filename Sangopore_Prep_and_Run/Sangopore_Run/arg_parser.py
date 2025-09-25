#!/usr/bin/env python3

import argparse
import os
from pathlib import Path
from .version import __version__

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run Sangopore Amplicon Pipeline")
    parser.add_argument(
        "-s", "--scriptPath",
        required=False,
        help="The path to the 'amplicon_sorter.py' script (default: look in the same directory as this script)"
    )
    parser.add_argument("--input_tsv", type=str, required=True, help='Path to the TSV file with multi-inputs')
    parser.add_argument("-o", "--output_folder", required=True, help="Final output directory")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads")
    parser.add_argument("--blast_db", type=str, required=True, help="Path to the BLAST nt database")
    parser.add_argument("--blast-files-prefix", type=str, default="", help="Prefix for BLAST database files (e.g., 'nt' if your files are nt.*)")
    parser.add_argument("--max-reads", type=int, default=10000, help='Set max reads for amplicon_sorter (default: 10000)')
    parser.add_argument("--use-ram-disk", action="store_true",help="Use a RAM disk to copy the BLAST database for faster processing (default: False)")
    parser.add_argument("--set-ram-size", type=float, default=None, help="Optional fixed RAM disk size in GB. If not set, automatically computes 1.1x the database folder size.")
    parser.add_argument("--run_fasta_filtering", action="store_true", help="Run postprocessing step to filter FASTA files based on CSV abundance threshold")
    parser.add_argument("--fasta_filter_threshold", type=float, default=1.0, help="Threshold for filtering FASTA based on Percent_Amp_reads_of_analysis_inputs_per_BC (default: 1.0)")
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args()

    # If scriptPath not provided, set default to amplicon_sorter.py in this script's folder
    if not args.scriptPath:
        args.scriptPath = str(Path(__file__).parent / "amplicon_sorter.py")

    # Validate it exists
    if not os.path.isfile(args.scriptPath):
        parser.error(f"amplicon_sorter.py script not found at '{args.scriptPath}'")

    return args

if __name__ == "__main__":
    args = parse_arguments()
    print(args)
