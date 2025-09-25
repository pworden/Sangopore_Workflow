# Sangopore Workflow Part 1 - <span style="color:#2e86c1">Sangopore Preparation</span>

## Description
**Sangopore Preparation** is the <u>first</u> part of the **Sangopore Workflow**. It requires a tab-separated value (TSV) file as input, and through a user determined output path, creates a number of subdirectories containing quality filtered read files using [chopper](https://github.com/wdecoster/chopper) and [fastp](https://github.com/OpenGene/fastp), quality control summaries through [MultiQC](https://github.com/MultiQC/MultiQC) and [TidyMultiqc](https://cran.r-project.org/web/packages/TidyMultiqc/vignettes/TidyMultiqc.html), as well as auto-generating a TSV file that will be used as input for part 2 of the Sangopore Workflow.

### Sangopore-Prep command options:
   - `--input_tsv` (The path to the TSV that will be used as input fo the `sangopore-prep` command)
   - `--parent_output_folder` (The path to the parent directory under which all output folders and files will be placed)
   - `-h` or `--help` (To show script command options)
   - `--version` (Current script version)

```
sangopore-prep \
    --input_tsv "Sangopore_Prep_and_Run_Dev/Test_Barcoded_Fastq_Inputs/Test_sango_prep_table.tsv" \
    --parent_output_folder "Sangopore_Prep_and_Run_Dev/Test_Outputs"
```

#### Input Table (TSV)
| barcode_number | input_barcode_file                                                 | subgroup_folder      | min_range_bp | max_range_bp |
|----------------|---------------------------------------------------------------------|----------------------|--------------|--------------|
| barcode01      | Sangopore_Prep_and_Run_Dev/Test_Barcoded_Fastq_Inputs/barcode01        | Group1     | 100          | 550          |
| barcode03      | Sangopore_Prep_and_Run_Dev/Test_Barcoded_Fastq_Inputs/barcode03        | Group1     | 100          | 550          |
| barcode09      | Sangopore_Prep_and_Run_Dev/Test_Barcoded_Fastq_Inputs/barcode09        | Group1     | 500          | 1000         |

<small>**Note:** When constructing an input table the column names must be the same as above or the script will fail.</small>

#### Table description
The above table sets up three different analyses for 2 barcodes, with each analysis using different input folder/file structures (see **"Column-2"** heading below). Once the table is constructed save it in tab-separated format (TSV) and use this path for table input (`--input_tsv`).

- **Column-1:** Barcode number in the format "barcode" + a number in the format "01", "02", etc. ... "99"
   - So for barcodes 1 and 3 the rows would contain:
      - barcode01
      - barcode03

- **Column-2:** Path to the parent directory containing files for a specific barcode.

- **Column-3:** The name of the subgroup folder or folders. The script will put the data into one or more subfolders based on how many unique names are in this column. If there is only to be one subfolder keep all names the same for this columm as in the above table.

- **Column-4:** The expected minimum length of the target sequence.

- **Column-5:** The expected maximum length of the target sequence.

---

### Workflow
**The input table (TSV) determines:**
1. The barcodes to be analysed.
1. The parent folder containing either:
   - A set of barcode subfolders.
   - Or, one subfolder with one fastq (read) file for each barcode.
1. Locates the reads for each barcode and allows the script to concatenate reads for each barcode, if needed.
1. Determines the number of analyses (3 in the example table above)
1. Sets the size range for later analyses.
   - *To be acted upon in part-2 of the Sangopore workflow*.

**With these data the script continues:**

1. Runs quality and size filtering through [chopper](https://github.com/wdecoster/chopper) and [fastp](https://github.com/OpenGene/fastp) on reads for each barcode.
1. For each analyses, quality filtered reads and quality control summaries are saved.
   - Subfolder (analysis) names are defined in the `--input_tsv` (input table) in the column `subgroup_folder`.
      - In the TSV example above the subfolder names would be:
         - 1_marteilia_multiDir_multiFile
         - 2_marteilia_multiDir_oneFile
         - 3_marteilia_oneDir_oneFile
1. MultiQC summaries (HTML) are saved for each analyses as well as an overall summary (HTML)
   - In the TSV example above there would be three subfolders each with a Multiqc summary, as well as a single Multiqc summary in the parent output folder
1. A summary table (CSV) is also created in the parent output folder from the MultiQC summaries.
   - Created through `tidy_multiqc_script.R` invoking  [TidyMultiqc](https://cran.r-project.org/web/packages/TidyMultiqc/vignettes/TidyMultiqc.html).
1. TSV tables for running part-2 of the sangopore workflow are saved.
   - One TSV table for each analysis
   - An overall table for all separate analyses (as determine by table input)
      - In the obove example this would include all three analyses
         - 1_marteilia_multiDir_multiFile
         - 2_marteilia_multiDir_oneFile
         - 3_marteilia_oneDir_oneFile

## Example output

**Parent_Directory**:
├── multiqc_data (Dir)
├── multiqc_report.html
├── new_folder_paths.txt
├── Sangopore_test_output_sangopore_prep_qc.tsv
└── summary_multiqc_table.csv

#### File Descriptions
- **multiqc_report.html**: Contains quality control information for all analyses
- **new_folder_paths.txt**: Paths to folders containing data for each analysis.
- **Sangopore_test_output_sangopore_prep_qc.tsv**: Input TSV used in the next stage of the Sangopore workflow (with direct paths to fastq files).
- **summary_multiqc_table.csv**: Summary of quality control information for all analyses

---
---

## Script Breakdown

1. **Initialization**
   - The `TSVSubsetter` class is initialized with attributes to store paths, dataframes, and processed file paths.

2. **Argument Parsing**
   - The `argument_parser()` method uses `argparse` to parse command-line arguments for the script. These arguments are stored as instance variables.

3. **Reading the Input TSV**
   - The `read_input_tsv()` method reads the input TSV file into a pandas DataFrame.

4. **Getting Unique Values**
   - The `get_unique_values()` method extracts unique values from the `subgroup_folder` column of the DataFrame.

5. **Creating the Parent Output Folder**
   - The `create_parent_output_folder()` method creates the parent output directory if it does not already exist.

6. **Creating Subdirectories**
   - The `create_subdirectories()` method creates subdirectories for each unique value in the `subgroup_folder` column. It also creates subset TSV files for each subgroup.

7. **Writing New Folder Paths**
   - The `write_new_folder_paths()` method writes paths to the newly created subdirectories into a text file.

8. **Processing Subset TSV Files**
   - The `process_subset_tsv_files()` method reads the folder paths from the text file, processes each subset TSV file by concatenating FASTQ files, running `chopper` and `fastp` QC tools, and generating summary outputs with `multiqc`.

9. **Running the R Script**
   - The `run_r_script()` method executes the provided R script to process the final QC summary data.

10. **Final Output**
    - The final output includes processed TSV files, QC results, and a summary generated by the R script.

