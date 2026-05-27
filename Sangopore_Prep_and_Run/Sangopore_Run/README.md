# Sangopore Workflow Part 2 - Sangopore Analysis

## Description
**Sangopore Analysis** is the <ins>second</ins> part of the **Sangopore Workflow**. It uses Amplicon Sorter (https://github.com/avierstr/amplicon_sorter) to discover one or more consensus sequences per barcoded set of reads, followed by taxonomic identification by BLASTN. Summary files are also produced that describe the  qualities and relative abundances of ampicon detected per barcode.

#### Sangopore Analysis `(sangopore-run)` command options:
   - `--input_tsv` (essentially the same as the TSV for the `sangopore-prep` command but with direct paths to the fastq files filtered by `sangopore-prep`)

      <details>
      <summary>Click here to see table example</summary>

      | barcode_number | fastq_file_paths                                                                                           | subgroup_folder | min_range_bp | max_range_bp |
      |----------------|------------------------------------------------------------------------------------------------------------|-----------------|--------------|--------------|
      | barcode01      | Sangopore_Prep_and_Run_Dev/Test_Outputs/All_M24-13575/All_M24-13575/barcode01_All_M24-13575_fastp.fastq.gz     | All_M24-13575   | 100          | 550          |
      | barcode03      | Sangopore_Prep_and_Run_Dev/Test_Outputs/All_M24-13575/All_M24-13575/barcode03_All_M24-13575_fastp.fastq.gz     | All_M24-13575   | 100          | 550          |
      | barcode09      | Sangopore_Prep_and_Run_Dev/Test_Outputs/All_M24-13575/All_M24-13575/barcode09_All_M24-13575_fastp.fastq.gz     | All_M24-13575   | 500          | 1000         |
      
      </details>

   - `-o` or `--output_folder` (The path to the output folder.)
        - Best to use the same output folder as for the previous `sangopore-prep` command (`--parent_output_folder`).
   - `-t` or `--threads` (Set the number of computer threads/cores)
   - `--blast_db` (Path to the BLAST nt database)
   - `--max-reads` (How many reads will be used in calculating all amplicons per barcode)
   - `-h` or `--help` (To show script command options)
   - `--version` (Current script version)

**Example command using test data:**

```
sangopore-run \
    --input_tsv "Sangopore_Prep_and_Run_Dev/Test_Outputs/Test_Outputs_sangopore_prep_qc.tsv" \
    -o "Sangopore_Prep_and_Run_Dev/Test_Outputs" \
    -t "40" \
    --blast_db "/path/to/local/blast_nt/database/nt/nt" \
    --max-reads 1000
```

<small>**Note:** If running the above test commands the user must be in the parent package directory (default = Sangopore_Prep_and_Run_Dev).</small>

### Workflow summary
**The input table (TSV) determines:**
1. The barcodes to be analysed.
1. The input file for each barcode and which subdirectory (analysis) to save the results

**With these data the script continues and will:**
1. Run [amplicon_sorter](https://github.com/avierstr/amplicon_sorter) on each barcoded fastq file, where:
   - A Histogram (per fastq) is created to examine the input reads.
   - Consensus sequence or sequences (amplicon/s) are determined, polished, and saved.
1. Perform a local BLASTN search to taxonomically identify each detected consensus sequence/s
1. Produce **BLASTN** outputs (each barcode)
   - TSV file of blastn results
   - XML file of blastn results for visualisation with Geneious
1. Output a Blastn summary
   - Concatenation of all blastn result tables (TSVs) into one summary TSV

---

## Post-Analysis: Remove Primer Sequences from Consensus FASTA Files

After running `sangopore-run`, you can optionally remove primer sequences from the consensus amplicon FASTA files.

### `sangopore_remove_primers.py`

This standalone script removes primer sequences from multi-FASTA files produced by Sangopore using regex-based matching. It is useful when:
- Your amplicons contain PCR primers that should be removed before downstream analysis
- You have primers at sequence ends or anywhere within the sequence
- You need to handle degenerate (ambiguous) primer sequences

#### Script location:
```
Sangopore_Prep_and_Run/Sangopore_Run/sangopore_remove_primers.py
```

#### Usage:

```bash
# Activate conda environment:
conda activate Sangopore_Prep_and_Run

# Run the script (minimal - uses all defaults):
python sangopore_remove_primers.py \
    --primers-fasta-file-path "/path/to/primers.fasta" \
    --output-dir "/path/to/sangopore/output"
```

#### Arguments:

- `--primers-fasta-file-path` (required) — Path to FASTA file containing primer sequences (forward and/or reverse)
  - **Fully supports degenerate IUPAC codes:** R, Y, W, S, K, M, N, D, H, B, V, I
  - Example: `GCGCGCAGATAG` or `RTTAATTCCGGAA` (R = A or G) or `CCIGAYATRGCITTYCCICG` (I = A/C/G)
- `--output-dir` (required) — Parent output directory containing the input amplicon folder
- `--input-amplicon-multifasta-folder` (optional, default: `amplicon_sorter_run`) — Input folder name containing consensus FASTA files
  - For per-barcode sequences: `amplicon_sorter_run`
  - For multi-FASTA files: `amplicon_sorter_consensus_multiFasta`
  - Or any other folder name with consensus sequences
- `--input-search-suffix` (optional, default: `_fp_consensussequences.fasta`) — File suffix to search for
  - Default works for standard per-barcode consensus files
  - Use `_multiFasta.fasta` for multi-FASTA consensus files
- `--output-file-base-suffix` (optional, default: `_fp_consensussequences_pr_trim.fasta`) — Suffix for output files
  - Controls output filename (replaces input suffix)
  - Example: `_pr_trim.fasta` or `_trimmed.fasta`

#### What it does:

1. Recursively finds all FASTA files matching the suffix in `amplicon_sorter_run/`
2. Searches for primer sequences **anywhere** in each consensus sequence using regex patterns
3. Removes found primers while **preserving flanking sequences** (including multiple primers per sequence)
4. Outputs trimmed sequences to: `amplicon_sorter_run_trimmed/`
5. Creates a `primer_removal.log` file per barcode with detailed trimming statistics

#### Key Features:

- ✅ Finds primers anywhere in sequence (not just at 5'/3' ends)
- ✅ Removes only the primer, preserves flanking sequences
- ✅ Handles multiple primers in the same sequence correctly
- ✅ Natively supports degenerate IUPAC codes
- ✅ Pure Python implementation (no external tool dependencies)
- ✅ Detailed per-sequence logging

#### Output structure:

```
parent-dir/
├── amplicon_sorter_run/          (original files)
│   ├── barcode01/
│   │   └── barcode01_fp_consensussequences.fasta
│   └── barcode02/
│       └── barcode02_fp_consensussequences.fasta
└── amplicon_sorter_run_trimmed/  (primer-removed output)
    ├── barcode01/
    │   ├── barcode01_fp_consensussequences.fasta
    │   └── primer_removal.log
    └── barcode02/
        ├── barcode02_fp_consensussequences.fasta
        └── primer_removal.log
```

#### Examples:

**Basic usage (minimal - all defaults):**
```bash
conda activate Sangopore_Prep_and_Run

python Sangopore_Prep_and_Run/Sangopore_Run/sangopore_remove_primers.py \
    --primers-fasta-file-path "path/to/my_primers.fasta" \
    --output-dir "path/to/sangopore/output"
```
Output: `amplicon_sorter_run_pr_trim/`

**Multi-FASTA consensus files:**
```bash
conda activate Sangopore_Prep_and_Run

python Sangopore_Prep_and_Run/Sangopore_Run/sangopore_remove_primers.py \
    --primers-fasta-file-path "path/to/my_primers.fasta" \
    --output-dir "path/to/sangopore/output" \
    --input-amplicon-multifasta-folder "amplicon_sorter_consensus_multiFasta" \
    --input-search-suffix "_multiFasta.fasta"
```
Output: `amplicon_sorter_consensus_multiFasta_pr_trim/`

**With custom output file suffix (per-barcode):**
```bash
python Sangopore_Prep_and_Run/Sangopore_Run/sangopore_remove_primers.py \
    --primers-fasta-file-path "path/to/my_primers.fasta" \
    --output-dir "path/to/sangopore/output" \
    --output-file-base-suffix "_pr_trim.fasta"
```

**With degenerate primers (primers.fasta):**
```fasta
>primer_forward_exact
GCGCGCAGATAG
>primer_reverse_degenerate
RTTAATTCCGGAA
>primer_internal
CCNNGGAA
```

IUPAC codes automatically interpreted:
- `R` = A or G (puRine)
- `Y` = C or T (pYrimidine)
- `N` = any base
- `W` = A or T (Weak)
- `S` = C or G (Strong)
- `K` = G or T (Keto)
- `M` = A or C (aMino)
- `D` = A, G, or T (not C)
- `H` = A, C, or T (not G)
- `B` = C, G, or T (not A)
- `V` = A, C, or G (not T)

#### Example Results:

| Scenario | Original | Trimmed | Removed |
|----------|----------|---------|---------|
| End primers | `GCGCGCAGATAGCGATACGATACGA` | `CGATACGATACGA` | 12bp primer |
| Middle primer | `AAAAAAAGCGCGCAGATAGAAAAAA` | `AAAAAAAAAAAAAA` | 12bp (keeps both sides) |
| Multiple primers | `GCGCGCAGATAGCCCCGCGCGCAGATAG` | `CCCC` | 24bp (2 primers) |

---

## Workflow: Trim Primers and Continue BLASTN Analysis

To trim primers from amplicons and then run BLASTN analysis:

### Step 1: Run primer removal
```bash
conda activate Sangopore_Prep_and_Run

python Sangopore_Prep_and_Run/Sangopore_Run/sangopore_remove_primers.py \
    --primers-fasta-file-path "path/to/primers.fasta" \
    --output-dir "path/to/sangopore/output" \
    --input-amplicon-multifasta-folder "amplicon_sorter_consensus_multiFasta_untrimmed"
```

### Step 2: Rename output folder
The script creates `amplicon_sorter_consensus_multiFasta_untrimmed_pr_trim/`. Rename this to `amplicon_sorter_consensus_multiFasta/` to replace the untrimmed sequences:
```bash
mv amplicon_sorter_consensus_multiFasta_untrimmed_pr_trim amplicon_sorter_consensus_multiFasta
```

### Step 3: Rename files (remove _pr_trim suffix)
```bash
cd amplicon_sorter_consensus_multiFasta
for f in *_pr_trim.fasta; do mv "$f" "${f%_pr_trim.fasta}.fasta"; done
```

### Step 4: Comment out Phase 1 in Sangopore_ampliconSorter.py
In `Sangopore_Prep_and_Run/Sangopore_Run/Sangopore_ampliconSorter.py`, comment out Phase 1 (lines 987-995) to skip amplicon_sorter and go directly to BLASTN:

```python
# ---------- Phase 1: Preprocessing ----------
# for subgroup, files in self.subgroup_data.items():
#     self.copy_compressed_fastq_files(subgroup, files)
#     self.uncompress_fastq_files(subgroup, files)
#     self.amplicon_sorter_QC(subgroup, files)
#     self.amplicon_sorter_run(subgroup, files)
#     self.amplicon_sorter_copy_fastqs(subgroup)
#     self.collect_and_rename_consensus_sequences(subgroup)
#     self.collect_and_rename_single_consensus_sequences(subgroup)
```

### Step 5: Run sangopore-run
Now run `sangopore-run` normally. It will skip Phase 1 and proceed directly to BLASTN analysis using your trimmed consensus sequences:

```bash
sangopore-run \
    --input_tsv "path/to/input.tsv" \
    -o "path/to/output" \
    -t "60" \
    --blast_db "$HOME/databases/blast_dbs/nt" \
    --blast-files-prefix "nt"
```
