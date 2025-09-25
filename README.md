# Sangopore

## Overview:
Sangopore is a workflow that takes Oxford Nanopore reads, and for reads belonging to each barcode, finds one or more consensus DNA sequences and taxonomically identifies each consensus sequence.

**The Sangopore workflow consists of two parts**

1. **Sangopore preparation**
   - Quality control and basic filtering of reads
   - Automated table creation for later input into <ins>Sangopore Analysis</ins>
1. **Sangopore Analysis**
   - For each set of barcoded reads [Amplicon Sorter](https://github.com/avierstr/amplicon_sorter) finds consensus DNA sequences
   - The taxonomic identity of each discovered consensus sequence is determined using BLASTN
   - Summary documents are also output

   ---

## Installation instructions 
1. Git clone this workflow
   - `git clone https://github.com/pworden/Sangopore_Workflow.git`
1. Use the `environment.yml` file to generate a conda environment called *Sangopore_Prep_and_Run*
   - `conda env create -f /path/to/Sangopore_Workflow/environment.yml`
1. Activate the newly created conda environment
   - `conda activate Sangopore_Prep_and_Run`
1. Move into the newly (github) cloned 'Sangopore_Workflow' directory on your computer
   - Check this directory contains a `pyproject.toml` file
1. Install the Sangopore Workflow package
   - `pip install -e .`
     - Or (if you want to give the full path)
   - `pip install -e /Path/to/Sangopore_Workflow`
1. Download and install blastn database
   - Check how many gz files there are for the nt database, at the following site
   - https://ftp.ncbi.nlm.nih.gov/blast/db/
      - (For example at writing the largest number is: nt.243.tar.gz)
   - Copy the following onto the linux command line, once you know the maximum number of nt database files, for example if "max_number=243":
      - wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt.{000..243}.tar.gz -P $HOME/blast_nt
      - Where "-P" directs the download to a directory called "blast_nt" (or whichever path is chosen)
        - The directory will be created if it doesn't exist
   - Move into the directory of newly downloaded nt files (.tar.gz) such as "blast_nt"
      - `cd /path/to/blast_nt`
   - Run the following command in the linux terminal to decompress the database:
      ```
      for f in nt.*.tar.gz; do
         tar -xvzf "$f"
      done
      ```
1. Download and install taxonomy database into the same directory as the blastn database (see above)
     - `wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P $HOME/blast_nt`
   - Move into the "blast_nt" directory, if not already there
     - `cd $HOME/blast_nt`
   - Expand the files so they can be accessed
     - `tar -xvzf taxdump.tar.gz`
1. With the following commands check that the install has been successful for the following scripts:
   - `sango-create-table-onedir-allBCs -h`
   (Partially create input table if all bar-coded fastq files are in one directory)
   - `sango-create-table-multidir-multifiles -h`
   (Partially create input table for multiple subdirectories, each representing one barcode)
   - `sangopore-prep -h`
   (Prepare, organise, and check quality of input data)
   - `sangopore-run -h`
   (For each set of bar-coded reads, find one or more amplicons, and taxonomically identify each. Finally, collate results for each barcode as well as grouping barcode results by study.)
   - `sangopore-results-sort`
   (Sort results based on barcode)
   - `sangopore-expand-bc-names -h`
   (Using the study input metadata, as exampled by "turorial.csv" within the "Test_Barcoded_Fastq_Inputs" folder. This script will substitute all barcodeXX (XX = 01 to 99) for a concatenation of barcode number and sample name defined by the metadata)
   - `sango-filter-fasta-with-blastCSV`
   (Will subset the FASTA amplicon file based on amplicon sequences that meet or exceed a specified abundance threshold in the blastn-CSV output)

---

## Run Scripts

### Pre-Sangopore preparation - create initial input table
Before the Sangopore preparation workflow can be run a table (TSV file) must be created that directs the inputs, script logic, and outputs of the Sangopore Preparation script.

For a detailed explanation of the initial **input table** see the [create initial input table](Sangopore_Prep_and_Run/Sangopore_Preparation/Sangopore_create_input_table/README.md) README.

- The table can be created either:
   1. <ins>Manually by following the TSV example here:</ins>
   `Sangopore_Prep_and_Run/Test_Barcoded_Fastq_Inputs/Test_sango_prep_table.tsv`

      <details>
      <summary>Click to see table example</summary>

      | barcode_number | input_barcode_file                                                 | subgroup_folder      | min_range_bp | max_range_bp |
      |----------------|---------------------------------------------------------------------|----------------------|--------------|--------------|
      | barcode01      | Sangopore_Prep_and_Run/Test_Barcoded_Fastq_Inputs/barcode01        | Group1     | 300          | 500          |
      | barcode03      | Sangopore_Prep_and_Run/Test_Barcoded_Fastq_Inputs/barcode03        | Group1     | 300          | 500          |
      | barcode09      | Sangopore_Prep_and_Run/Test_Barcoded_Fastq_Inputs/barcode09        | Group2     | 700          | 900         |

      </details>

   1. <ins>Or, by using one of two scripts to partially auto-create the initial input TSV</ins>
      - If all fastq files are in one directory with one fastq file per barcode:
         ```
         # Activate conda environment:
         conda activate Sangopore_Prep_and_Run

         # TEST SCRIPT EXAMPLE
         sango-create-table-onedir-allBCs \
            --parent_input_dir "Sangopore_Prep_and_Run/Test_OneDir_BC_per_fastq" \
            --output_dir "Sangopore_Prep_and_Run/Test_Outputs" \
            --output_tsv_name "AllInOneDir_table.tsv"
         # Use full paths for "--parent_input_dir" and "--output_dir" for user data
         ```

      - If there are multiple folders, each folder representing one barcode, and possibly many fastq files per barcode:
         ```
         # Activate conda environment:
         conda activate Sangopore_Prep_and_Run

         # TEST SCRIPT EXAMPLE
         sango-create-table-multidir-multifiles \
            --parentInputDir "Sangopore_Prep_and_Run/Test_Barcoded_Fastq_Inputs" \
            --outputDir "Sangopore_Prep_and_Run/Test_Outputs" \
            --outTableName "MultiDirMultiFiles_table.tsv"
         # Use full paths for "--parentInputDir" and "--outputDir" for user data
         ```
<sup>ℹ️ **Note:** If running the above test commands the user must be in the parent package directory (default = Sangopore_Prep_and_Run).</sup>
<sup>❗ **Note:** Activate conda environment </sup>

   - If using either of the two auto table-creation scripts, the remaining three columns of data still need to be filled in ("subgroup_folder", "min_range_bp", "max_range_bp"). Different character strings (using underscores for spaces) with the "subgroup_folder" column are used to divide output into one or more sub-folders.
     - For example if the "subgroup_folder" column has either "Run1_Group1" or "Run1_Group2" in each row, script outputs for each barcode will be into either a subdirectory named "Run1-Group1" or "Run1_Group2".

---

### Sangopore Preparation

The `sangopore-prep` command invokes a python script that takes both a TSV file (described above) and a parent output path as inputs. The data in the TSV file allows the script to direct input file locations, expected read sizes, and organise QC and other outputs into one or more subdirectories.

For a detailed explanation see the [Sangopore Preparation](Sangopore_Prep_and_Run/Sangopore_Preparation/README.md) README.

#### sangopore-prep command options:
   - `--input_tsv` (The path to the TSV that will be used as input fo the `sangopore-prep` command)
   - `--parent_output_folder` (The path to the parent directory under which all output folders and files will be placed)
   - `-h` or `--help` (To show script command options)
   - `--version` (Current script version)

```
# Activate conda environment:
conda activate Sangopore_Prep_and_Run

sangopore-prep \
    --input_tsv "Sangopore_Prep_and_Run/Test_Barcoded_Fastq_Inputs/Test_sango_prep_table.tsv" \
    --parent_output_folder "Sangopore_Prep_and_Run/Test_Outputs"
```

<sup>ℹ️ **Note:** If running the above test commands the user must be in the parent package directory (default = Sangopore_Prep_and_Run).</sup>

---

### Sangopore analysis using Amplicon Sorter

The `sangopore-run` command takes a TSV file output from the `sangopore-prep` command, and uses it as input to direct Amplicon Sorter to find consensus DNA sequences, for each barcode, from oxford nanopore reads. The TSV file is in the `sangopore-prep` parent output directory, with the suffix `*_sangopore_prep_qc.tsv`.

For a detailed explanation see the [Sangopore Analysis](Sangopore_Prep_and_Run/Sangopore_Run/README.md) README.

#### sangopore-run command options:
   - `--input_tsv` (essentially the same as the TSV for the `sangopore-prep` command but with direct paths to the fastq files filtered by `sangopore-prep`)
   - `-o` or `--output_folder` (The path to the output folder.)
        - Best to use the same output folder as for the previous `sangopore-prep` command (`--parent_output_folder`).
   - `-t` or `--threads` (Set the number of computer threads/cores)
   - `--blast_db` (Path to the BLAST nt database)
   - `--max-reads` (How many reads will be used in calculating all amplicons per barcode)
   - `-h` or `--help` (To show script command options)
   - `--version` (Current script version)

**Example command using test data:**
```
# If you want to just run the workflow:

# Activate conda environment:
conda activate Sangopore_Prep_and_Run

sangopore-run \
    --input_tsv "Sangopore_Prep_and_Run/Test_Outputs/Test_Outputs_sangopore_prep_qc.tsv" \
    -o "Sangopore_Prep_and_Run/Test_Outputs" \
    -t "40" \
    --blast_db "/path/to/local/blast_nt/database/nt/nt" \
    --max-reads 1000 \
    --run_fasta_filtering \
    --fasta_filter_threshold 1

# ----------

# Activate conda environment:
conda activate Sangopore_Prep_and_Run

# Per barcode, if you also want to filter the amplicons based on percent abundance relative to read input for the analysis (default = 1):
sangopore-run \
    --input_tsv "Sangopore_Prep_and_Run/Test_Outputs/Test_Outputs_sangopore_prep_qc.tsv" \
    -o "Sangopore_Prep_and_Run/Test_Outputs" \
    -t "40" \
    --blast_db "/path/to/local/blast_nt/database/nt/nt" \
    --max-reads 1000 \
    --run_fasta_filtering \
    --fasta_filter_threshold 1
```

<sup>ℹ️ **Note:** If running the above test commands the user must be in the parent package directory (default = Sangopore_Prep_and_Run).</sup>

---
## Metadata gathering scripts

### `sangopore-results-sort`

This is a module that is automatically invoked by the "sangopore-run" command. It takes results from the Sangopore run, and for each barcode gathers specific outputs into a folder, copying relevant FASTQ, PDF, FASTA, Excel, and BLASTN output files and folders into their corresponding barcode folders. However it can also be run as a standalone script, which if required, is invoked by the command below.

```
# To run as a standalone script:
conda activate Sangopore_Prep_and_Run

sangopore-results-sort \
    --parent-input-output-dir "Sangopore_Prep_and_Run/Test_Outputs/sangoAnalysis_out/Group1"
```

<sup>ℹ️ **Note:** If running the above test commands the user must be in the parent package directory (default = Sangopore_Prep_and_Run).</sup>

<p align="center">❖❖❖❖❖❖❖❖❖❖</p>

### `sangopore-expand-bc-names`

This is a stand alone command that takes as input a CSV file manually captured from the Excel Sangopore booking sheet (copy sheet into new document, save as CSV and copy to server - note the CSV path and use it as input for `--input-barcode-and-id-csv`). It then extracts barcodes and sample identifier, concatenates them and uses this data to add sample data to barcode names in folder and filenames. A CSV of identifiers, barcodes, and concatenated strings (identifier + barcode) is also output.

```
# Activate conda environment:
conda activate Sangopore_Prep_and_Run

sangopore-expand-bc-names \
    --input-barcode-and-id-csv "Sangopore_Prep_and_Run/Test_Barcoded_Fastq_Inputs/tutorial.csv" \
    --data-parent-dir "Sangopore_Prep_and_Run/Test_Outputs/sangoAnalysis_out/Group1/results_per_barcode"
    
    # To direct the output to a specific path add the following command:
      --output-csv "/some/other/path/cleaned_sample_names.csv"
```

<sup>ℹ️ **Note:** If running the above test commands the user must be in the parent package directory (default = Sangopore_Prep_and_Run).</sup>

<p align="center">❖❖❖❖❖❖❖❖❖❖</p>

### `filter_fasta_from_blastCSV`

This script filters multi-amplicon FASTA files based on amplicon sequence abundance from corresponding BLASTN result CSVs. For each barcode, it retains amplicon sequences that meet or exceed a specified abundance threshold and writes the filtered FASTA to the output directory.

```
# Activate conda environment:
conda activate Sangopore_Prep_and_Run

sango-filter-fasta-with-blastCSV \
    --csv_dir "/home/wordenp/projects/sangopore_project/analyses/Sangopore_analysis/SP25-006_Sangopore_17JUN25_max10000/sangoAnalysis_out/Biotech/amplicon_sorter_blastn_20_ordered_collation" \
    --fasta_dir "/home/wordenp/projects/sangopore_project/analyses/Sangopore_analysis/SP25-006_Sangopore_17JUN25_max10000/sangoAnalysis_out/Biotech/amplicon_sorter_consensus_multiFasta" \
    --output_dir "/home/wordenp/projects/sangopore_project/analyses/Sangopore_analysis/SP25-006_Sangopore_17JUN25_max10000/sangoAnalysis_out/Biotech/amplicon_sorter_consensus_multiFasta_abund_filtered" \
    --threshold 1.0
```

---

### Output folder and file structure

**If the above test scripts were run, the parent output folder should contain the following:**

- **Group1** (Directory - of outputs for the group named in the `Sangopore_Prep` input TSV. In this instance only one group was named in the test data TSV)
  - If the `Sangopore_Prep` input [TSV defined multiple groups](#barcoded-table), multiple subfolders for each group would also be present.
- **multiqc_data** (Directory - output of quality control data from MultiQC data)
- **multiqc_report.html** (File - The MultiQC report in HTML format)
- **new_folder_paths.txt** (File - Text file of subdirectory paths, with each representing a path to an output group as defined within the input TSV for `Sangopore_Prep`.)
- **sangoAnalysis_out** (Directory - The parent folder containing all output for the `Sangopore_Run` command, as defined by the `-o` or `--output_folder` flag)
- **summary_multiqc_table.csv** (CSV version of `multiqc_report.html`)
- **Test_Outputs_sangopore_prep_qc.tsv** (This file with the "_sangopore_prep_qc.tsv" suffix will be used as the TSV input for the `Sangopore_Run` command.)

### Full folder and file structure of test output:

```
Parent_Output_Folder     # Path named by the user
├── Group1
│   ├── Group1
│   │   ├── barcode01_Group1_fastp.fastq.gz
│   │   ├── barcode01_Group1_fastp.html
│   │   ├── barcode01_Group1_fastp.json
│   │   ├── barcode03_Group1_fastp.fastq.gz
│   │   ├── barcode03_Group1_fastp.html
│   │   ├── barcode03_Group1_fastp.json
│   │   ├── barcode09_Group1_fastp.fastq.gz
│   │   ├── barcode09_Group1_fastp.html
│   │   └── barcode09_Group1_fastp.json
│   ├── Group1_directPath_subset.tsv
│   ├── Group1_subset.tsv
│   ├── multiqc_data
│   │   ├── fastp_filtered_reads_plot.txt
│   │   ├── fastp-seq-content-gc-plot_Read_1_After_filtering.txt
│   │   ├── fastp-seq-content-gc-plot_Read_1_Before_filtering.txt
│   │   ├── fastp-seq-content-n-plot_Read_1_After_filtering.txt
│   │   ├── fastp-seq-content-n-plot_Read_1_Before_filtering.txt
│   │   ├── fastp-seq-quality-plot_Read_1_After_filtering.txt
│   │   ├── fastp-seq-quality-plot_Read_1_Before_filtering.txt
│   │   ├── multiqc_citations.txt
│   │   ├── multiqc_data.json
│   │   ├── multiqc_fastp.txt
│   │   ├── multiqc_general_stats.txt
│   │   ├── multiqc.log
│   │   ├── multiqc_software_versions.txt
│   │   └── multiqc_sources.txt
│   └── multiqc_report.html
├── multiqc_data
│   ├── fastp_filtered_reads_plot.txt
│   ├── fastp-seq-content-gc-plot_Read_1_After_filtering.txt
│   ├── fastp-seq-content-gc-plot_Read_1_Before_filtering.txt
│   ├── fastp-seq-content-n-plot_Read_1_After_filtering.txt
│   ├── fastp-seq-content-n-plot_Read_1_Before_filtering.txt
│   ├── fastp-seq-quality-plot_Read_1_After_filtering.txt
│   ├── fastp-seq-quality-plot_Read_1_Before_filtering.txt
│   ├── multiqc_citations.txt
│   ├── multiqc_data.json
│   ├── multiqc_fastp.txt
│   ├── multiqc_general_stats.txt
│   ├── multiqc.log
│   ├── multiqc_software_versions.txt
│   └── multiqc_sources.txt
├── multiqc_report.html
├── new_folder_paths.txt
├── sangoAnalysis_out
│   └── Group1
│       ├── amplicon_sorter_blastn_20_ordered_collation
│       ├── amplicon_sorter_blastn_20_ordered__excel
│       ├── amplicon_sorter_blastn_collation
│       ├── amplicon_sorter_blastn_out
│       ├── amplicon_sorter_consensus_multiFasta
│       ├── amplicon_sorter_consensus_singleFasta
│       ├── Amplicon_Sorter_QC
│       ├── amplicon_sorter_run
│       ├── compressed
│       ├── results_per_barcode # <----- Folder of final per-barcode results
│       └── uncompressed
├── summary_multiqc_table.csv
├── Test_Outputs_README.md
└── Test_Outputs_sangopore_prep_qc.tsv
```

<sup>ℹ️ **Note:** In the above folder structure the multiqc outputs are duplicates, with multiqc results on level 1 (upper) and level 2 (lower) levels of the folder structure. However, if the analysis consisted of multiple groups then level 1 would contain all multiqc results, while level 2 would contain multiqc results specific for each group defined in the TSV input.</sup>

<p align="center">❖❖❖❖❖❖❖❖❖❖</p>

### The per-barcode collated data in the  `results_per_barcode` folder
- The sample information in file and folder names is a result of barcode information gathered from `sangopore-run` accessing `sangopore-results-sort`, as well as sample information from `sangopore-expand-bc-names`.

```
path/to/folder/results_per_barcode
├── barcode01_sample1
│   ├── barcode01_sample1_fp_consensussequences.fasta
│   ├── barcode01_sample1_fp_total_outputfig.pdf
│   ├── barcode01_sample1_Group1_fastp.fastq.gz
│   ├── blastn_barcode01_sample1
│   │   ├── blastn_barcode01_sample1.tsv
│   │   └── blastn_barcode01_sample1.xml
│   └── blastn_collation_barcode01_sample1_ordered20.xlsx
├── barcode03_sample2
│   ├── barcode03_sample2_fp_consensussequences.fasta
│   ├── barcode03_sample2_fp_total_outputfig.pdf
│   ├── barcode03_sample2_Group1_fastp.fastq.gz
│   ├── blastn_barcode03_sample2
│   │   ├── blastn_barcode03_sample2.tsv
│   │   └── blastn_barcode03_sample2.xml
│   └── blastn_collation_barcode03_sample2_ordered20.xlsx
├── barcode09_sample3
│   ├── barcode09_sample3_fp_consensussequences.fasta
│   ├── barcode09_sample3_fp_total_outputfig.pdf
│   ├── barcode09_sample3_Group1_fastp.fastq.gz
│   ├── blastn_barcode09_sample3
│   │   ├── blastn_barcode09_sample3.tsv
│   │   └── blastn_barcode09_sample3.xml
│   └── blastn_collation_barcode09_sample3_ordered20.xlsx
│       └── barcode09_sample3_fp_0_0__57_consensus.fasta
└── cleaned_sample_names.csv
```

### Multiple Output subfolders - example

The `Sangopore_Prep` input TSV can be set to give multiple subfolder outputs (groups). Simply assign required names to each row within the `subgroup_folder` column as shown below (Group1 and Group2).

<a name="barcoded-table"></a>
<details>

***<summary>Click to expand table</summary>***

| barcode_number | input_barcode_file                                                 | subgroup_folder      | min_range_bp | max_range_bp |
|----------------|---------------------------------------------------------------------|----------------------|--------------|--------------|
| barcode01      | Sangopore_Prep_and_Run/Test_Barcoded_Fastq_Inputs/barcode01        | Group1     | 100          | 550          |
| barcode03      | Sangopore_Prep_and_Run/Test_Barcoded_Fastq_Inputs/barcode03        | Group2     | 100          | 550          |
| barcode09      | Sangopore_Prep_and_Run/Test_Barcoded_Fastq_Inputs/barcode09        | Group2     | 500          | 1000         |

</details>

Within the `Sangopore_Prep` input TSV, in the `subgroup_folder` column are two groups defined by the `subgroup_folder` column (Group1 and Group2). This would result in a folder structure as follows:

```
Parent_Output_Folder     # Path named by the user
├── Group1
│   ├──Group1_QC_Files_and_Folders
├── Group2
│   ├──Group2_QC_Files_and_Folders
├── multiqc_data
├── multiqc_report.html
├── new_folder_paths.txt
├── sangoAnalysis_out # All output from the Sangopore_Run command 
├── summary_multiqc_table.csv
└── Test_Outputs_sangopore_prep_qc.tsv
```
