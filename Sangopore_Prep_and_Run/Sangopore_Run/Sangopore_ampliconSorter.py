#!/usr/bin/env python3

"""
# Github test command:
conda activate Sangopore_Prep_and_Run

# Can run this first to delete any existing ram disk (although)
rm -rf /mnt/ramfs2/nt
# Create the nt folder ram disk
mkdir /mnt/ramfs2/nt

sangopore-run \
    --input_tsv "Sangopore_Prep_and_Run/Test_Outputs/Test_Outputs_sangopore_prep_qc.tsv" \
    -o "Sangopore_Prep_and_Run/Test_Outputs" \
    -t "40" \
    --blast_db "$HOME/databases/blast_dbs/nt" \
    --blast-files-prefix "nt" \
    --max-reads 10000 \
    --use-ram-disk \
    --run_fasta_filtering \
    --fasta_filter_threshold 1
"""

from .arg_parser import parse_arguments
from pathlib import Path
import os, subprocess, csv, sys, re
import pandas as pd
from Bio import Phylo
from Bio import SeqIO
import time
import math
import matplotlib.pyplot as plt
from openpyxl.styles import PatternFill
from openpyxl.utils import get_column_letter
from openpyxl import load_workbook
import shutil
import threading
from Sangopore_Prep_and_Run.Sangopore_Run import final_sort


def create_parent_output_folder():
    os.makedirs(output_folder, exist_ok=True)

def read_tsv(tsv_path):
    data = []
    with open(tsv_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        headers = next(reader)
        for row in reader:
            data.append(row)
    return headers, data

def group_by_subfolder(data):
    grouped_data = {}
    for row in data:
        subgroup = row[2]
        if subgroup not in grouped_data:
            grouped_data[subgroup] = []
        grouped_data[subgroup].append(row)
    return grouped_data

def create_output_folders(output_folder, subgroups):
    for subgroup in subgroups:
        subgroup_path = os.path.join(output_folder, subgroup)
        if not os.path.exists(subgroup_path):
            os.makedirs(subgroup_path, exist_ok=True)

class Sangopore:
    def __init__(self,
                 script_path: str,
                 output_folder: str,
                 threads: int,
                 subgroup_data: dict,
                 blast_db_path: str,
                 max_reads: int = 500,
                 blast_files_prefix: str = "",
                 run_fasta_filtering: bool = False,
                 fasta_filter_threshold: float = 1.0,
                 use_ram_disk: bool = True,
                 fixed_ram_size_gb: float = None
                 ):
        self.script_path = script_path
        self.output_folder = output_folder
        self.threads = threads
        self.subgroup_data = subgroup_data
        self.blast_db_path = Path(blast_db_path)
        self.blast_files_prefix = blast_files_prefix
        self.max_reads = max_reads
        self.run_fasta_filtering = run_fasta_filtering
        self.fasta_filter_threshold = fasta_filter_threshold
        self.use_ram_disk = use_ram_disk

        # For ram disk usage
        self.ramdisk_path = Path("/mnt/ramfs2")
        self.ramdisk_copy_proc = None
        self.ram_db_path = None

        # Determine the actual path to BLAST files
        if self.use_ram_disk:
            # Use only the database path; the prefix will be appended when calling BLAST
            self.blast_path = Path(self.blast_db_path)
        else:
            # Use database path + prefix
            self.blast_path = Path(self.blast_db_path) / self.blast_files_prefix


    def copy_compressed_fastq_files(self, subgroup, files):
        compressed_path = os.path.join(self.output_folder, subgroup, "compressed")
        os.makedirs(compressed_path, exist_ok=True)

        for row in files:
            input_file = row[1]
            prefix = row[0]
            dest_file = os.path.join(compressed_path, os.path.basename(input_file))

            # Check if the file already exists to avoid overwriting
            if not os.path.exists(dest_file):
                shutil.copy2(input_file, dest_file)
                print(f"Copied compressed file to: {dest_file}")
            else:
                print(f"Compressed file already exists: {dest_file}")


    def uncompress_fastq_files(self, subgroup, files):
        uncompressed_path = os.path.join(self.output_folder, subgroup, "uncompressed")
        if not os.path.exists(uncompressed_path):
            os.makedirs(uncompressed_path, exist_ok=True)
        
        for row in files:
            input_file = row[1]
            prefix = row[0]
            output_file = os.path.join(uncompressed_path, prefix + "_fp.fastq")
            
            # Check if the file already exists to avoid reprocessing
            if not os.path.exists(output_file):
                with open(output_file, 'w') as output_file_handle:
                    gunzip_process = subprocess.Popen(['gunzip', '-c', input_file], stdout=subprocess.PIPE, text=True)
                    subprocess.run(['cat'], stdin=gunzip_process.stdout, stdout=output_file_handle, text=True)
                print(f"Uncompressed file created: {output_file}")
            else:
                print(f"File already exists: {output_file}")


    def amplicon_sorter_QC(self, subgroup, files):
        qc_path = os.path.join(self.output_folder, subgroup, "Amplicon_Sorter_QC")
        if not os.path.exists(qc_path):
            os.mkdir(qc_path)
        
        trimmed_path = os.path.join(self.output_folder, subgroup, "uncompressed")
        for row in files:
            prefix = row[0]
            trimmed_file = prefix + "_fp.fastq"
            if os.path.exists(os.path.join(trimmed_path, trimmed_file)):
                InputTrimmedFile = os.path.join(trimmed_path, trimmed_file)
                subprocess.run(['python', self.script_path, '-i', InputTrimmedFile, '-o', qc_path, '-ho', '-maxr', str(self.max_reads), '-np', str(self.threads)])


    def amplicon_sorter_run(self, subgroup, files):
        import math

        ampSort_path = os.path.join(self.output_folder, subgroup, "amplicon_sorter_run")
        os.makedirs(ampSort_path, exist_ok=True)

        trimmed_path = os.path.join(self.output_folder, subgroup, "uncompressed")
        cmds = []

        # Compute cores per file (at least 1 core per file)
        threads_per_file = max(1, self.threads // len(files))
        print(f"[INFO] Allocating {threads_per_file} cores per FASTQ file for {len(files)} files.")

        for row in files:
            prefix = row[0]
            trimmed_file = prefix + "_fp.fastq"
            InputTrimmedFile = os.path.join(trimmed_path, trimmed_file)

            if os.path.exists(InputTrimmedFile):
                minLength, maxLength = row[3], row[4]
                cmd = [
                    sys.executable, self.script_path, 
                    "-i", InputTrimmedFile,
                    "-o", ampSort_path, 
                    "-ldc", "0.5", 
                    "-min", str(minLength), 
                    "-max", str(maxLength), 
                    "-ar", 
                    "-maxr", str(self.max_reads), 
                    "-np", str(threads_per_file),  # cores allocated per file
                    "--save_fastq"
                ]
                cmds.append(" ".join(cmd))

        tmpfile = os.path.join(ampSort_path, ".amp_sort_cmds.tmp")
        with open(tmpfile, "w") as pcmds:
            pcmds.write("\n".join(cmds))

        # Run jobs in parallel; parallel will queue jobs if more than available cores
        subprocess.run(
            ["parallel", "--jobs", str(self.threads), "-a", tmpfile],
            check=True
        )


    def amplicon_sorter_copy_fastqs(self, subgroup):
        """
        - Copies individual FASTQ files listed in consensus FASTA
        - Concatenates all reads per barcode into *_concat.fastq (deduplicated)
        - Writes a read count summary per barcode
        """
        run_dir = os.path.join(self.output_folder, subgroup, "amplicon_sorter_run")
        out_dir = os.path.join(self.output_folder, subgroup, "amplicon_sorter_fastq_out")
        os.makedirs(out_dir, exist_ok=True)

        for barcode_folder in os.listdir(run_dir):
            barcode_path = os.path.join(run_dir, barcode_folder)
            if not os.path.isdir(barcode_path):
                continue

            fasta_file = os.path.join(barcode_path, f"{barcode_folder}_consensussequences.fasta")
            if not os.path.isfile(fasta_file):
                print(f"[WARNING] Consensus FASTA missing: {fasta_file}")
                continue

            dest_dir = os.path.join(out_dir, barcode_folder)
            os.makedirs(dest_dir, exist_ok=True)

            seen_sequences = set()
            concat_path = os.path.join(dest_dir, f"{barcode_folder}_concat.fastq")
            read_counts = {}

            with open(fasta_file, 'r') as f, open(concat_path, 'w') as concat_out:
                for line in f:
                    if not line.startswith(">"):
                        continue
                    clean_id = line.strip().lstrip(">").split("(")[0].replace("consensus_", "")
                    fastq_name = f"{clean_id}.fastq"
                    fastq_path = os.path.join(barcode_path, fastq_name)
                    copied_path = os.path.join(dest_dir, fastq_name)

                    if os.path.isfile(fastq_path):
                        # Copy to output directory
                        shutil.copy2(fastq_path, copied_path)
                        print(f"[COPIED] {fastq_name} → {dest_dir}")

                        # Count reads in this file
                        read_count = 0
                        with open(fastq_path, 'r') as fq:
                            while True:
                                lines = [fq.readline() for _ in range(4)]
                                if not lines[0]:
                                    break
                                read_count += 1
                                read_block = "".join(lines)
                                if read_block not in seen_sequences:
                                    concat_out.write(read_block)
                                    seen_sequences.add(read_block)

                        read_counts[fastq_name] = read_count
                    else:
                        print(f"[MISSING] {fastq_name} not found in {barcode_path}")
                        read_counts[fastq_name] = 0

            # Add total unique reads to summary
            read_counts[f"{barcode_folder}_concat.fastq"] = len(seen_sequences)

            # Write summary file
            summary_file = os.path.join(dest_dir, f"{barcode_folder}_read_counts.txt")
            with open(summary_file, 'w') as sf:
                sf.write("fastq_file\tread_count\n")
                for fname, count in sorted(read_counts.items()):
                    sf.write(f"{fname}\t{count}\n")

            print(f"[SUMMARY] {summary_file} created.")

        print(f"\nFASTQ copy, concatenation, and count summaries complete for subgroup: {subgroup}")


    def collect_and_rename_consensus_sequences(self, subgroup, search_term="_fp_consensussequences.fasta", output_subdir="amplicon_sorter_consensus_multiFasta"):
        """
        Finds all consensus FASTA files matching `search_term` in the amplicon_sorter_run subfolder for a subgroup,
        copies them to a shared folder (e.g. 'amplicon_sorter_consensus_multiFasta'), and prepends the prefix.
        """
        input_dir = os.path.join(self.output_folder, subgroup, "amplicon_sorter_run")
        output_dir = os.path.join(self.output_folder, subgroup, output_subdir)
        prefix = subgroup  # or any other prefix you prefer

        if not os.path.isdir(input_dir):
            print(f"[ERROR] Input directory does not exist: {input_dir}")
            return

        os.makedirs(output_dir, exist_ok=True)

        count = 0
        for root, _, files in os.walk(input_dir):
            for file in files:
                if search_term in file:
                    src = os.path.join(root, file)
                    dst = os.path.join(output_dir, file)
                    shutil.copy2(src, dst)
                    count += 1

        if count == 0:
            print(f"[INFO] No files matching '{search_term}' found in {input_dir}")
        else:
            print(f"[INFO] {count} consensus FASTA file(s) copied to {output_dir}")


    def collect_and_rename_single_consensus_sequences(
        self,
        subgroup,
        search_term="_consensussequences.fasta",
        output_subdir="amplicon_sorter_consensus_singleFasta",
        prefix=""
    ):
        """
        Finds all single-sequence consensus FASTA files matching `search_term` in the
        amplicon_sorter_run subfolder for a subgroup, processes and copies them to a
        shared folder (e.g. 'amplicon_sorter_consensus_singleFasta'), and prepends the prefix.

        Args:
            subgroup (str): Subfolder name inside output_folder.
            search_term (str): Pattern to match files.
            output_subdir (str): Subfolder name for output within subgroup.
            prefix (str): Prefix to prepend to output filenames.
        """
        input_dir = os.path.join(self.output_folder, subgroup, "amplicon_sorter_run")
        output_dir = os.path.join(self.output_folder, subgroup, output_subdir)

        if not os.path.isdir(input_dir):
            print(f"[ERROR] Input directory does not exist: {input_dir}")
            return

        os.makedirs(output_dir, exist_ok=True)

        import re
        from Bio import SeqIO

        pattern = re.compile(r"barcode\d{2}_fp_\d+_.*consensussequences\.fasta$")

        matched_files = [
            os.path.join(root, f)
            for root, _, files in os.walk(input_dir)
            for f in files
            if search_term in f and pattern.match(f)
        ]

        if not matched_files:
            print(f"[INFO] No matching consensus FASTA files found in: {input_dir}")
            return

        print(f"[INFO] Processing {len(matched_files)} consensus FASTA files from: {input_dir}")
        processed_count = 0

        for filepath in matched_files:
            records = list(SeqIO.parse(filepath, "fasta"))
            if not records:
                print(f"[WARNING] Skipping {os.path.basename(filepath)} (no sequences)")
                continue

            first_record = records[0]
            header_match = re.search(r'fp_(\d+_\d+)\((\d+)\)', first_record.description)
            if not header_match:
                print(f"[WARNING] Skipping {os.path.basename(filepath)} (bad header): {first_record.description}")
                continue

            fp_id = header_match.group(1)
            number = header_match.group(2)

            filename_match = re.match(r"(.*fp)_\d+_consensussequences\.fasta", os.path.basename(filepath))
            if not filename_match:
                print(f"[WARNING] Skipping {os.path.basename(filepath)} (unexpected filename structure)")
                continue

            base = filename_match.group(1)
            new_filename = f"{prefix}_{base}_{fp_id}__{number}_consensus.fasta" if prefix else f"{base}_{fp_id}__{number}_consensus.fasta"
            out_path = os.path.join(output_dir, new_filename)

            with open(out_path, "w") as out_f:
                SeqIO.write(first_record, out_f, "fasta")

            processed_count += 1
            print(f"[INFO] Saved: {new_filename}")

        print(f"[SUCCESS] {processed_count} single-sequence consensus FASTA files written to: {output_dir}")


    def blastn_run(self, subgroup, use_ramdisk=False):
        import os
        import subprocess
        from pathlib import Path
        import sys

        parent_input_folder = Path(self.output_folder) / subgroup / "amplicon_sorter_run"
        parent_output_folder = Path(self.output_folder) / subgroup / "amplicon_sorter_blastn_out"

        # Determine BLAST DB path
        if use_ramdisk:
            # Use RAM disk copy
            blast_db_for_run = self.ram_db_path / self.blast_files_prefix
        else:
            # Normal DB path: append prefix so BLAST sees ../nt/nt
            blast_db_for_run = Path(self.blast_db_path) / self.blast_files_prefix  # e.g., ../nt/nt

        # ---- Check that BLAST DB files exist ----
        blast_files = list(blast_db_for_run.parent.glob(blast_db_for_run.name + ".*"))
        if not blast_files:
            print(f"[ERROR] BLAST DB path does not exist or no files found: {blast_db_for_run}")
            return

        # Find input FASTA files
        common_string = "_fp_consensussequences"
        input_files = list(parent_input_folder.rglob(f"*{common_string}.fasta"))
        if not input_files:
            print(f"[WARNING] No files found with pattern {common_string} in {parent_input_folder}")
            return

        # Run BLAST for each input file
        for input_file in input_files:
            if not input_file.is_file():
                print(f"[WARNING] File not found: {input_file}")
                continue

            output_file_base = input_file.stem.replace(common_string, "")
            blast_out_path = parent_output_folder / ("blastn_" + output_file_base)
            blast_out_path.mkdir(parents=True, exist_ok=True)

            # ---- Run blastn (ASN output) ----
            blast_asn_file = blast_out_path / f"blastn_{output_file_base}.asn"
            blastn_command_asn = [
                "blastn",
                "-query", str(input_file.resolve()),
                "-db", str(blast_db_for_run.resolve()),
                "-outfmt", "11",
                "-max_target_seqs", "20",
                "-max_hsps", "2",
                "-evalue", "1e-05",
                "-num_threads", str(self.threads),
                "-out", str(blast_asn_file.resolve())
            ]

            print(f"[DEBUG] Running BLAST ASN:\n{' '.join(blastn_command_asn)}")
            try:
                subprocess.run(blastn_command_asn, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                print(f"[INFO] BLAST ASN output written to {blast_asn_file}")
            except subprocess.CalledProcessError as e:
                print(f"[ERROR] BLAST failed for {input_file.name}")
                print(f"Return code: {e.returncode}")
                print(f"stdout: {e.stdout.decode('utf-8')}")
                print(f"stderr: {e.stderr.decode('utf-8')}")
                continue

            # ---- Convert ASN → TSV ----
            blast_tsv_file = blast_out_path / f"blastn_{output_file_base}.tsv"
            subprocess.run([
                "blast_formatter",
                "-archive", str(blast_asn_file),
                "-outfmt", "6 qseqid sseqid qcovs qlen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids stitle sscinames",
                "-out", str(blast_tsv_file)
            ], check=True)

            # ---- Convert ASN → XML ----
            blast_xml_file = blast_out_path / f"blastn_{output_file_base}.xml"
            subprocess.run([
                "blast_formatter",
                "-archive", str(blast_asn_file),
                "-outfmt", "5",
                "-out", str(blast_xml_file)
            ], check=True)

            print(f"[INFO] TSV and XML outputs written to {blast_out_path}")


# ------
    def load_passed_filter_reads_from_multiqc(self):
        """
        Loads passed filter read counts and maps barcode numbers (e.g. 78) to read counts.
        """

        summary_path = os.path.join(os.path.dirname(self.output_folder), "summary_multiqc_table.csv")

        if not os.path.exists(summary_path):
            print(f"[WARNING] summary_multiqc_table.csv not found at expected path: {summary_path}")
            return {}

        try:
            df = pd.read_csv(summary_path)

            # Extract barcode number as integer (e.g. 78 from 'barcode78_...')
            df["barcode_num"] = df["metadata.sample_id"].str.extract(r"barcode(\d+)")[0].astype(int)

            # Create a mapping: 78 → 21173
            passed_reads_lookup = dict(zip(df["barcode_num"], df["general.filtering_result_passed_filter_reads"]))

            print(f"[INFO] Loaded passed filter reads for {len(passed_reads_lookup)} barcodes.")
            return passed_reads_lookup

        except Exception as e:
            print(f"[ERROR] Could not parse summary_multiqc_table.csv: {e}")
            return {}


# -----
    def process_blastn_files(self, subgroup):
        parent_collation_input_dir = Path(os.path.join(self.output_folder, subgroup, "amplicon_sorter_blastn_out"))
        print(f"[DEBUG] Running BLASTN collation for subgroup: {subgroup}")
        print(f"[DEBUG] Looking in folder: {parent_collation_input_dir}")

        passed_reads_lookup = self.load_passed_filter_reads_from_multiqc()

        tsv_files = list(parent_collation_input_dir.rglob("*.tsv"))
        print(f"[DEBUG] Found TSV files: {[str(f) for f in tsv_files]}")

        if not tsv_files:
            print("No TSV files found to process.")
            return

        parent_collation_output_dir = Path(os.path.join(self.output_folder, subgroup, "amplicon_sorter_blastn_collation"))
        parent_collation_output_dir.mkdir(parents=True, exist_ok=True)

        def process_file(tsv_path, output_dir):
            filename = os.path.basename(tsv_path)
            outfile = output_dir / filename.replace("blastn_", "blastn_collation_").replace(".tsv", ".csv")
            print(f"Processing: {filename}")

            headers = [
                "Query_Seq_id",
                "Subject_Seq_id",
                "Query_Coverage",
                "Query_Length",
                "Percent_identical_matches",
                "Alignment_length",
                "Number_mismatches",
                "Number_gap_openings",
                "Start_alignment_in_query",
                "End_alignment_in_query",
                "Start_alignment_in_subject",
                "End_alignment_in_subject",
                "Expect_value",
                "Bit_score",
                "Subject_Taxonomy_ID",
                "Subject_Title",
                "Subject_scientific_name"
            ]

            df = pd.read_csv(tsv_path, sep="\t", header=None)
            if df.shape[1] != len(headers):
                raise ValueError(f"Expected {len(headers)} columns, got {df.shape[1]} in file {filename}")
            df.columns = headers

            # Extract base Query_Seq_id and read counts
            def extract_parts(query_id):
                match = re.match(r"(.+?)\((\d+)\)", query_id)
                if match:
                    return match.group(1), int(match.group(2))
                return query_id, 1  # fallback

            df[["Query_Seq_id", "Read_count"]] = df["Query_Seq_id"].apply(lambda x: pd.Series(extract_parts(x)))

            # Compute total reads per base ID
            abundance_lookup = (
                df.groupby("Query_Seq_id")["Read_count"]
                .first()  # Each base ID should appear with a consistent count
                .to_dict()
            )
            total_reads = sum(abundance_lookup.values())

            # Map percent abundance to each row based on Query_Seq_id
            df["Percent_abundance"] = df["Query_Seq_id"].map(
                lambda x: (abundance_lookup.get(x, 0) / total_reads * 100) if total_reads > 0 else 0
            )

            # Add  -----
            barcode_num_match = re.search(r"barcode(\d+)", os.path.basename(tsv_path))
            barcode_num = int(barcode_num_match.group(1)) if barcode_num_match else None
            print(f"[DEBUG] Extracted barcode_num: {barcode_num}")

            print(f"[DEBUG] subgroup: {subgroup}")
            print(f"[DEBUG] Extracted barcode_num: {barcode_num}")
            print(f"[DEBUG] Keys in passed_reads_lookup: {list(passed_reads_lookup.keys())[:10]}")
            print(f"[DEBUG] Passed_Filter_Reads found: {passed_reads_lookup.get(barcode_num, 'NOT FOUND')}")

            # Apply to all rows in the DataFrame
            df["set_input_maxreads"] = self.max_reads
            df["barcode_number"] = barcode_num
            df["Passed_Filter_Reads"] = passed_reads_lookup.get(barcode_num, "NA")

            # Reorder columns to match requested output
            col_order = [
                "Query_Seq_id",
                "Subject_scientific_name",
                "Subject_Taxonomy_ID",
                "Subject_Title",
                "Read_count",
                "Passed_Filter_Reads",
                "Percent_abundance",
                "set_input_maxreads",
                "Subject_Seq_id",
                "Query_Coverage",
                "Query_Length",
                "Percent_identical_matches",
                "Alignment_length",
                "Number_mismatches",
                "Number_gap_openings",
                "Start_alignment_in_query",
                "End_alignment_in_query",
                "Start_alignment_in_subject",
                "End_alignment_in_subject",
                "Expect_value",
                "Bit_score"
            ]

            df = df[col_order]
            df.to_csv(outfile, index=False)
            return outfile

        for tsv_file in tsv_files:
            try:
                csv_path = process_file(tsv_file, parent_collation_output_dir)
                print(f"Written CSV: {csv_path}")
            except Exception as e:
                print(f"Error processing {tsv_file}: {e}")


    def sort_blastn_collated_csvs(self, subgroup):
        from pathlib import Path
        import pandas as pd

        input_dir = Path(self.output_folder) / subgroup / "amplicon_sorter_blastn_collation"
        output_dir = Path(self.output_folder) / subgroup / "amplicon_sorter_blastn_20_ordered_collation"
        output_dir.mkdir(parents=True, exist_ok=True)

        drop_cols = [
            "Number_mismatches",
            "Number_gap_openings",
            "Start_alignment_in_query",
            "End_alignment_in_query",
            "Start_alignment_in_subject",
            "End_alignment_in_subject"
        ]

        csv_files = list(input_dir.glob("*.csv"))
        print(f"[DEBUG] Sorting and trimming {len(csv_files)} CSV(s). Output to: {output_dir}")

        for csv_file in csv_files:
            try:
                df = pd.read_csv(csv_file)

                # Drop unwanted columns if present
                df = df.drop(columns=[col for col in drop_cols if col in df.columns])

                # Ensure numeric columns are correct type
                for col in ["set_input_maxreads", "Passed_Filter_Reads", "Read_count"]:
                    if col in df.columns:
                        df[col] = pd.to_numeric(df[col], errors="coerce")

                grouped = []
                for query_id, group in df.groupby("Query_Seq_id"):
                    group_sorted = group.sort_values(
                        by=["Percent_abundance", "Query_Coverage", "Bit_score"],
                        ascending=[False, False, False]
                    ).head(20).copy()
                    # Rename Percent_abundance for clarity
                    group_sorted = group_sorted.rename(columns={
                        "Percent_abundance": "Percent_abundance_of_Amplicon_reads"
                    })

                    read_count = group["Read_count"].unique()
                    if len(read_count) == 1:
                        read_count = read_count[0]
                    else:
                        raise ValueError(f"Inconsistent Read_count for {query_id}: {read_count}")
                    passed_filter = group["Passed_Filter_Reads"].iloc[0]
                    set_input = group["set_input_maxreads"].iloc[0]

                    if pd.notna(read_count) and pd.notna(passed_filter) and pd.notna(set_input):
                        denom = min(passed_filter, set_input)
                        if denom > 0:
                            percent_amplicon_reads = (read_count / denom) * 100
                        else:
                            percent_amplicon_reads = None
                    else:
                        percent_amplicon_reads = None

                    group_sorted["Percent_Amp_reads_of_analysis_inputs_per_BC"] = percent_amplicon_reads

                    print(f"[DEBUG] Query_ID: {query_id}, Read_count: {read_count}, "
                        f"set_input: {set_input}, passed_filter: {passed_filter}, "
                        f"percent_amplicon_reads: {percent_amplicon_reads}")

                    max_abundance = group_sorted["Percent_abundance_of_Amplicon_reads"].max()
                    grouped.append((max_abundance, query_id, group_sorted))

                # Sort and combine all top 20 rows across amplicons
                grouped_sorted = sorted(grouped, key=lambda x: (-x[0], x[1]))
                df_sorted = pd.concat([g[2] for g in grouped_sorted], ignore_index=True)

                # Write to output file
                output_path = output_dir / f"{csv_file.stem}_ordered20.csv"
                df_sorted.to_csv(output_path, index=False)
                print(f"[DEBUG] Wrote: {output_path.name}")

            except Exception as e:
                print(f"[ERROR] Could not process {csv_file.name}: {e}")


    def write_grouped_excel_with_shading(self, subgroup):
        from pathlib import Path
        import pandas as pd
        from openpyxl import load_workbook
        from openpyxl.styles import PatternFill

        # Define input and output directories
        input_dir = Path(self.output_folder) / subgroup / "amplicon_sorter_blastn_20_ordered_collation"
        output_dir = Path(self.output_folder) / subgroup / "amplicon_sorter_blastn_20_ordered__excel"
        output_dir.mkdir(parents=True, exist_ok=True)

        drop_cols = [
            "Number_mismatches",
            "Number_gap_openings",
            "Start_alignment_in_query",
            "End_alignment_in_query",
            "Start_alignment_in_subject",
            "End_alignment_in_subject"
        ]

        csv_files = list(input_dir.glob("*.csv"))
        print(f"[DEBUG] Writing shaded Excel files from {len(csv_files)} CSV(s). Output to: {output_dir}")

        # Define fill colors
        grey_fill = PatternFill(start_color="D9D9D9", end_color="D9D9D9", fill_type="solid")  # light grey header
        blue_fill = PatternFill(start_color="DDEBF7", end_color="DDEBF7", fill_type="solid")  # light blue
        yellow_fill = PatternFill(start_color="FFF2CC", end_color="FFF2CC", fill_type="solid")  # light yellow

        for csv_file in csv_files:
            try:
                df = pd.read_csv(csv_file)

                # Drop unwanted columns if present
                df = df.drop(columns=[col for col in drop_cols if col in df.columns])

                # Rename abundance column
                if "Percent_abundance" in df.columns:
                    df = df.rename(columns={"Percent_abundance": "Percent_abundance_of_Amplicon_reads"})

                # Group and sort each group individually
                grouped = []
                for query_id, group in df.groupby("Query_Seq_id"):
                    group_sorted = group.sort_values(
                        by=["Percent_abundance_of_Amplicon_reads", "Query_Coverage", "Bit_score"],
                        ascending=[False, False, False]
                    ).head(20)
                    max_abundance = group_sorted["Percent_abundance_of_Amplicon_reads"].max()
                    grouped.append((max_abundance, query_id, group_sorted))

                # Sort the groups by max Percent_abundance
                grouped_sorted = sorted(grouped, key=lambda x: (-x[0], x[1]))

                # Concatenate all groups into one DataFrame
                df_sorted = pd.concat([g[2] for g in grouped_sorted], ignore_index=True)

                output_path = output_dir / f"{csv_file.stem}.xlsx"
                df_sorted.to_excel(output_path, index=False)

                # Load workbook and sheet to apply styling
                wb = load_workbook(output_path)
                ws = wb.active

                # Apply grey fill to header
                for cell in ws[1]:
                    cell.fill = grey_fill

                # Apply alternating fill by Query_Seq_id
                current_fill = blue_fill
                previous_id = None
                for row in ws.iter_rows(min_row=2, max_row=ws.max_row):
                    query_id = row[0].value  # assuming Query_Seq_id is first column
                    if query_id != previous_id:
                        current_fill = yellow_fill if current_fill == blue_fill else blue_fill
                        previous_id = query_id
                    for cell in row:
                        cell.fill = current_fill

                # Auto-adjust column widths
                for col in ws.columns:
                    max_length = 0
                    column_letter = col[0].column_letter
                    for cell in col:
                        if cell.value:
                            max_length = max(max_length, len(str(cell.value)))
                    ws.column_dimensions[column_letter].width = max_length + 2

                wb.save(output_path)
                print(f"[DEBUG] Excel file written: {output_path.name}")
            except Exception as e:
                print(f"[ERROR] Failed to process {csv_file.name}: {e}")


    # ----------------------
    # RAM-disk helper methods
    # ----------------------
    def mount_ramdisk_user(self, ramdisk_path: Path):
        """
        Mount RAM disk at the given path as a normal user.
        This blocks until the mount completes.
        """
        ramdisk_path.mkdir(parents=True, exist_ok=True)
        print(f"[DEBUG] Mounting RAM disk at {ramdisk_path} (user)")
        subprocess.run(["mount", str(ramdisk_path)], check=True)  # no sudo required
        print(f"[INFO] RAM disk mounted at {ramdisk_path}")
        
        
    def copy_db_to_ramdisk(self, db_folder: Path, ramdisk_path: Path):
        """
        Copy BLAST DB to RAM disk in background thread.
        """
        def copy():
            self.ram_db_path = ramdisk_path  # <- use the path as-is
            print(f"[INFO] Copying BLAST DB to RAM disk: {self.ram_db_path}")
            shutil.copytree(db_folder, self.ram_db_path, dirs_exist_ok=True)
            print("[INFO] BLAST DB copy complete")

        self.ramdisk_copy_thread = threading.Thread(target=copy)
        self.ramdisk_copy_thread.start()


    def wait_for_ramdisk_copy(self):
        """Wait for background DB copy to complete before starting BLAST."""
        if hasattr(self, "ramdisk_copy_thread") and self.ramdisk_copy_thread.is_alive():
            print("[INFO] Waiting for RAM disk copy to finish...")
            self.ramdisk_copy_thread.join()
            print("[INFO] RAM disk ready")


    def cleanup_ramdisk(self):
        """Unmount and clean up the RAM disk."""
        if not self.use_ram_disk:
            return  # Skip if RAM disk was not used

        if not self.ramdisk_path.exists():
            print("[DEBUG] RAM disk path does not exist, skipping cleanup.")
            return

        print("[DEBUG] Cleaning up RAM disk...")
        try:
            subprocess.run(["sudo", "umount", str(self.ramdisk_path)], check=True)
        except subprocess.CalledProcessError as e:
            print(f"[WARNING] Could not unmount RAM disk: {e}")

        try:
            shutil.rmtree(self.ramdisk_path)
        except Exception as e:
            print(f"[WARNING] Could not remove RAM disk folder: {e}")


    def process_fasta_files(self, subgroup):
        parent_folder = os.path.join(self.output_folder, subgroup, "amplicon_sorter_run")
        outfolder = os.path.join(self.output_folder, subgroup, "amplicon_sorter_alnAndTree")
        parent_folder = Path(parent_folder)
        outfolder = Path(outfolder)

        # Check if the parent folder exists
        if not os.path.exists(parent_folder):
            print(f"Error: The parent folder '{parent_folder}' does not exist.")
            return


        def adjust_fasta_headers(file_path):
            with open(file_path, 'r') as infile:
                lines = infile.readlines()
            with open(file_path, 'w') as outfile:
                for line in lines:
                    if line.startswith('>'):
                        line = line.replace('(', '__').replace(')', '')
                    outfile.write(line)


        def align_and_create_tree(fasta_path, output_folder):
            # Ensure the output folder exists
            os.makedirs(output_folder, exist_ok=True)
            # Define the output file base name
            base_name = os.path.basename(fasta_path).replace("fp_consensussequences.fasta", "")
            # Define the paths for the alignment and tree output
            alignment_path = os.path.join(output_folder, f'{base_name}_alignment.fasta')
            tree_path = os.path.join(output_folder, f'{base_name}_tree.newick')
            pdf_tree_path = os.path.join(output_folder, f'{base_name}_tree.pdf')
            
            # Check if the alignment and tree files already exist
            if os.path.exists(alignment_path) or os.path.exists(tree_path):
                print(f"Output files already exist for {fasta_path}. Skipping processing.")
                return

            try:
                # Run MAFFT to create the alignment
                mafft_command = f"mafft --auto {fasta_path} > {alignment_path}"
                print(f"Running MAFFT: {mafft_command}")
                result = subprocess.run(mafft_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                print(f"MAFFT stdout: {result.stdout.decode('utf-8')}")
                print(f"MAFFT stderr: {result.stderr.decode('utf-8')}")
                if result.returncode != 0:
                    print(f"MAFFT failed with exit code {result.returncode}")
                    return
                print(f"MAFFT alignment created: {alignment_path}")

                # Adjust headers in the alignment file
                adjust_fasta_headers(alignment_path)
                print(f"Headers adjusted in: {alignment_path}")

                # Run FastTree to create the Newick tree from the alignment
                fasttree_command = f"FastTree -nt {alignment_path} > {tree_path}"
                print(f"Running FastTree: {fasttree_command}")
                result = subprocess.run(fasttree_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                print(f"FastTree stdout: {result.stdout.decode('utf-8')}")
                print(f"FastTree stderr: {result.stderr.decode('utf-8')}")
                if result.returncode != 0:
                    print(f"FastTree failed with exit code {result.returncode}")
                    return
                print(f"FastTree output created: {tree_path}")

                # Load the tree and create a PDF
                tree = Phylo.read(tree_path, "newick")
                fig = plt.figure(figsize=(11.7, 8.3))  # Set the figure size to A4 landscape
                ax = fig.add_subplot(1, 1, 1)
                Phylo.draw(tree, do_show=False, axes=ax)
                # Adjust the text size and layout
                for text in ax.texts:
                    text.set_fontsize(8)  # Set text size to make it readable
                # Save the figure to a PDF file
                plt.savefig(pdf_tree_path, format='pdf', bbox_inches='tight')
                print(f"Tree PDF saved at: {pdf_tree_path}")

            except subprocess.CalledProcessError as e:
                print(f"Error during subprocess execution: {e}")
                print(f"Command output: {e.output}")

        # Ensure the output directory exists
        os.makedirs(outfolder, exist_ok=True)
        found_any_fasta_files = False  # Flag to check if any relevant files are found
        for root, dirs, files in os.walk(parent_folder):
            print(f"Traversing directory: {root}")
            for file in files:
                print(f"Found file: {file}")
                if file.endswith("fp_consensussequences.fasta"):
                    found_any_fasta_files = True
                    fasta_path = os.path.join(root, file)
                    # Create a corresponding output directory without nested subdirectories
                    base_name = file.replace("fp_consensussequences.fasta", "")
                    output_folder = os.path.join(outfolder, base_name + "_tree")
                    print(f"Processing file: {fasta_path}, output folder: {output_folder}")
                    align_and_create_tree(fasta_path, output_folder)

        if not found_any_fasta_files:
            print("No files ending with 'fp_consensussequences.fasta' were found.")


    def run_pipeline(self):
        if self.use_ram_disk:
            ramdisk_path = Path("/mnt/ramfs2")

            # 1. Mount the RAM disk (blocks until ready)
            self.mount_ramdisk_user(ramdisk_path)
            
            # 1a. Remove old database on RAM disk if present
            ramdisk_nt_path = ramdisk_path / "nt"
            if ramdisk_nt_path.exists():
                print(f"[INFO] Removing old RAM disk database at {ramdisk_nt_path}")
                shutil.rmtree(ramdisk_nt_path)

            # 2. Copy the database in background
            self.copy_db_to_ramdisk(Path(self.blast_db_path), ramdisk_nt_path)
            
            # 3. Tell the object where the RAM disk copy is
            self.ram_db_path = ramdisk_nt_path

        try:
            # ---------- Phase 1: Preprocessing ----------
            for subgroup, files in self.subgroup_data.items():   ### CHANGED
                self.copy_compressed_fastq_files(subgroup, files)
                self.uncompress_fastq_files(subgroup, files)
                self.amplicon_sorter_QC(subgroup, files)
                self.amplicon_sorter_run(subgroup, files)
                self.amplicon_sorter_copy_fastqs(subgroup)
                self.collect_and_rename_consensus_sequences(subgroup)
                self.collect_and_rename_single_consensus_sequences(subgroup)

            if self.use_ram_disk:   ### CHANGED: wait once after preprocessing
                self.wait_for_ramdisk_copy()

            # ---------- Phase 2: BLAST + Post-BLAST ----------
            for subgroup, files in self.subgroup_data.items():   ### CHANGED
                self.blastn_run(subgroup, use_ramdisk=self.use_ram_disk)

                self.process_blastn_files(subgroup)
                self.sort_blastn_collated_csvs(subgroup)
                self.write_grouped_excel_with_shading(subgroup)
                final_sort.ResultsCollation(
                    Path(self.output_folder) / subgroup
                ).run_pipeline()
                # self.process_fasta_files(subgroup)

                # Optional postprocessing: filter FASTA by abundance
                if self.run_fasta_filtering:
                    try:
                        filter_fasta_script_dir = Path(__file__).parent
                        filter_fasta_script_path = filter_fasta_script_dir / "filter_fasta_from_blastCSV.py"                 
                        csv_dir = os.path.join(self.output_folder, subgroup, "amplicon_sorter_blastn_20_ordered_collation")
                        fasta_dir = os.path.join(self.output_folder, subgroup, "amplicon_sorter_consensus_multiFasta")
                        output_dir = os.path.join(self.output_folder, subgroup, "amplicon_sorter_consensus_multiFasta_abund_filtered")

                        command = [
                            "python3", str(filter_fasta_script_path),
                            "--csv_dir", str(csv_dir),
                            "--fasta_dir", str(fasta_dir),
                            "--output_dir", str(output_dir),
                            "--threshold", str(self.fasta_filter_threshold)
                        ]

                        print(f"[DEBUG] Running FASTA filtering for subgroup '{subgroup}': {' '.join(command)}")
                        subprocess.run(command, check=True)

                    except subprocess.CalledProcessError as e:
                        print(f"[ERROR] FASTA filtering script failed for subgroup {subgroup}: {e}")

        finally:
            if self.use_ram_disk:
                print(f"[INFO] Cleaning up RAM disk {ramdisk_path}")
                subprocess.run(["umount", str(ramdisk_path)], check=False)


class TeeOutput:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, message):
        for stream in self.streams:
            stream.write(message)
            stream.flush()

    def flush(self):
        for stream in self.streams:
            stream.flush()


def main():
    start_time = time.time()  # Record start

    args = parse_arguments()
    # Create a subdirectory called "sangoAnalysis_out" inside the user-specified output folder
    analysis_output_folder = os.path.join(args.output_folder, "sangoAnalysis_out")
    
    # Create log file inside sangoAnalysis_out
    log_file_path = Path(analysis_output_folder) / "sangopore_stdout.log"
    log_file_path.parent.mkdir(parents=True, exist_ok=True)  # Ensure folder exists
    
    # Redirect stdout and stderr
    sys.stdout = TeeOutput(sys.stdout, open(log_file_path, "w"))
    sys.stderr = sys.stdout
    
    # Read and process the TSV file
    headers, data = read_tsv(args.input_tsv)
    # Group data by subgroup folder
    subgroup_data = group_by_subfolder(data)
    # Create parent output folder and subfolders (including "sangoAnalysis_out")
    create_output_folders(analysis_output_folder, subgroup_data.keys())
    # Initialize and run the Sangopore pipeline
    sangopore_run = Sangopore(
        args.scriptPath, analysis_output_folder, args.threads, subgroup_data, args.blast_db, args.max_reads,
        use_ram_disk=args.use_ram_disk,
        fixed_ram_size_gb=args.set_ram_size,
        run_fasta_filtering=args.run_fasta_filtering,
        fasta_filter_threshold=args.fasta_filter_threshold,
        blast_files_prefix=args.blast_files_prefix
    )
    print("[DEBUG] Subgroup data keys:")
    print(list(subgroup_data.keys()))
    sangopore_run.run_pipeline()

    # Record end time
    end_time = time.time()
    elapsed_seconds = end_time - start_time
    print(f"[INFO] Total runtime: {elapsed_seconds:.2f} seconds "
          f"({elapsed_seconds/60:.2f} minutes)")

if __name__ == "__main__":
    main()
