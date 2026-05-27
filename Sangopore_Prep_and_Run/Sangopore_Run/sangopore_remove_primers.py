#!/usr/bin/env python3

"""
Regex-based primer removal from consensus FASTA sequences.
Supports degenerate IUPAC primer codes and removes primers anywhere in sequences.

EXAMPLES:

1. Per-barcode consensus sequences (default settings):
   conda activate Sangopore_Prep_and_Run
   python sangopore_remove_primers.py \
       --primers-fasta-file-path "primers.fasta" \
       --output-dir "path/to/sangopore/output"
   Output: amplicon_sorter_run_pr_trim/
   (Uses defaults: amplicon_sorter_run, _fp_consensussequences.fasta)

2. Multi-FASTA consensus files:
   python sangopore_remove_primers.py \
       --primers-fasta-file-path "primers.fasta" \
       --output-dir "path/to/sangopore/output" \
       --input-amplicon-multifasta-folder "amplicon_sorter_consensus_multiFasta" \
       --input-search-suffix "_multiFasta.fasta"
   Output: amplicon_sorter_consensus_multiFasta_pr_trim/

3. Custom output filename suffix:
   python sangopore_remove_primers.py \
       --primers-fasta-file-path "primers.fasta" \
       --output-dir "path/to/sangopore/output" \
       --output-file-base-suffix "_pr_trim.fasta"

PRIMER FILE FORMAT (supports degenerate IUPAC codes):
   >barcode60_forward
   CCIGAYATRGCITTYCCICG
   >barcode60_reverse
   TCDGGRTGICCRAARAAYCA

DEGENERATE CODES:
   R=A/G, Y=C/T, W=A/T, S=C/G, K=G/T, M=A/C, N=any
   D=A/G/T, H=A/C/T, B=C/G/T, V=A/C/G, I=A/C/G
"""

import re
import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq

# IUPAC degenerate base codes
IUPAC_CODES = {
    'R': '[AG]',      # Purine (A or G)
    'Y': '[CT]',      # Pyrimidine (C or T)
    'W': '[AT]',      # Weak (A or T)
    'S': '[CG]',      # Strong (C or G)
    'K': '[GT]',      # Keto (G or T)
    'M': '[AC]',      # Amino (A or C)
    'N': '[ACGT]',    # aNy base
    'D': '[AGT]',     # not C
    'H': '[ACT]',     # not G
    'B': '[CGT]',     # not A
    'V': '[ACG]',     # not T
    'I': '[ACG]',     # Inosine (A or C or G)
}

def parse_args():
    parser = argparse.ArgumentParser(
        description="Remove primer sequences from consensus FASTA files. "
        "Supports degenerate IUPAC codes and finds primers anywhere in sequences."
    )
    parser.add_argument("--primers-fasta-file-path", required=True,
                       help="Path to primers FASTA file (supports degenerate IUPAC codes: R,Y,W,S,K,M,N,D,H,B,V)")
    parser.add_argument("--output-dir", required=True,
                       help="Parent output directory containing the input amplicon folder")
    parser.add_argument("--input-amplicon-multifasta-folder", default="amplicon_sorter_run",
                       help="Input folder name containing consensus FASTA files (default: 'amplicon_sorter_run')")
    parser.add_argument("--input-search-suffix", default="_fp_consensussequences.fasta",
                       help="Suffix of FASTA files to search and trim (default: '_fp_consensussequences.fasta')")
    parser.add_argument("--output-file-base-suffix", default="_fp_consensussequences_pr_trim.fasta",
                       help="Suffix for output files (default: '_fp_consensussequences_pr_trim.fasta')")
    return parser.parse_args()

def iupac_to_regex(primer_seq):
    """Convert IUPAC degenerate primer sequence to regex pattern."""
    pattern = ""
    for base in primer_seq.upper():
        pattern += IUPAC_CODES.get(base, base)
    return pattern

def load_primers(primers_file):
    """Load primer sequences from FASTA file and convert to regex patterns."""
    primers = []
    for record in SeqIO.parse(primers_file, "fasta"):
        seq = str(record.seq)
        regex_pattern = iupac_to_regex(seq)
        primers.append({
            'name': record.id,
            'sequence': seq,
            'pattern': regex_pattern
        })
    return primers

def remove_primers_from_sequence(sequence, primers):
    """Remove all primer matches from sequence."""
    trimmed_seq = sequence
    matches_found = {}

    for primer in primers:
        pattern = primer['pattern']
        matches = list(re.finditer(pattern, trimmed_seq))
        if matches:
            matches_found[primer['name']] = len(matches)
        trimmed_seq = re.sub(pattern, "", trimmed_seq)

    return trimmed_seq, matches_found

def process_fasta_file(input_file, output_file, log_file, primers):
    """Process a single FASTA file, removing primers."""
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    total_seqs = 0
    seqs_with_primers = 0
    total_bases_removed = 0

    with open(output_file, "w") as out_fasta, open(log_file, "a") as log:
        log.write(f"\n==== Primer removal log for {input_file} ====\n")

        for record in SeqIO.parse(input_file, "fasta"):
            total_seqs += 1
            original_length = len(record.seq)

            trimmed_seq, matches = remove_primers_from_sequence(str(record.seq), primers)
            trimmed_length = len(trimmed_seq)
            bases_removed = original_length - trimmed_length

            if matches:
                seqs_with_primers += 1
                total_bases_removed += bases_removed
                log.write(f"{record.id}\t{original_length}bp\t{bases_removed}bp removed\t{matches}\n")

            record.seq = Seq(trimmed_seq)
            SeqIO.write(record, out_fasta, "fasta")

        log.write(f"\nSummary: {seqs_with_primers}/{total_seqs} sequences had primers\n")
        log.write(f"Total bases removed: {total_bases_removed}bp\n")

    print(f"Processed {input_file}: {seqs_with_primers}/{total_seqs} sequences trimmed ({total_bases_removed}bp removed)")

def find_files_with_suffix(parent_dir, suffix):
    """Recursively find all files under parent_dir ending with suffix."""
    matches = []
    for root, dirs, files in os.walk(parent_dir):
        for fname in files:
            if fname.endswith(suffix):
                matches.append(os.path.join(root, fname))
    return matches

def pipeline(primers_file, output_dir, input_amplicon_multifasta_folder, input_search_suffix, output_file_base_suffix):
    output_dir = os.path.abspath(output_dir)
    primers_file = os.path.abspath(primers_file)

    input_dir = os.path.join(output_dir, input_amplicon_multifasta_folder)
    output_base = os.path.join(output_dir, input_amplicon_multifasta_folder + "_pr_trim")

    if not os.path.isdir(input_dir):
        print(f"Error: Input directory {input_dir} does not exist.")
        return

    # Load primers once
    print(f"Loading primers from {primers_file}...")
    primers = load_primers(primers_file)
    print(f"Loaded {len(primers)} primer(s): {[p['name'] for p in primers]}")

    # Find files to process
    files_to_process = find_files_with_suffix(input_dir, input_search_suffix)
    if not files_to_process:
        print(f"Error: No files ending with '{input_search_suffix}' found under {input_dir}")
        return

    print(f"Processing {len(files_to_process)} file(s)...\n")

    # Process each file
    for input_path in files_to_process:
        # Build output paths maintaining directory structure
        rel_path = os.path.relpath(input_path, input_dir)

        # Replace input suffix with output suffix in filename
        if rel_path.endswith(input_search_suffix):
            base_path = rel_path[:-len(input_search_suffix)]
            output_filename = base_path + output_file_base_suffix
        else:
            output_filename = rel_path

        output_path = os.path.join(output_base, output_filename)

        barcode_folder = os.path.dirname(input_path)
        log_path = os.path.join(output_base, os.path.relpath(barcode_folder, input_dir), "primer_removal.log")

        process_fasta_file(input_path, output_path, log_path, primers)

    print(f"\nDone. Trimmed sequences saved to: {output_base}")

def main():
    args = parse_args()
    pipeline(args.primers_fasta_file_path, args.output_dir, args.input_amplicon_multifasta_folder,
             args.input_search_suffix, args.output_file_base_suffix)

if __name__ == "__main__":
    main()
