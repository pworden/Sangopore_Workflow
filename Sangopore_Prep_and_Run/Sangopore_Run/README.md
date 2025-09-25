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
