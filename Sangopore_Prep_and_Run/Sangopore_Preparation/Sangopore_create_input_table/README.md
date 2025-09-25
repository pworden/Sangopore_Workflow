
# Manual or partial autocreation of input tables for Sangopore

Before running the `sangopore-prep` command the user needs to create an input table with information on the barcode number, path to the barcode directory holding one or more fastq read files, one or more group names, and a minimum and maximum size range in base pairs.

**Input Table structure example:**
| barcode_number | input_barcode_file                                                 | subgroup_folder      | min_range_bp | max_range_bp |
|----------------|---------------------------------------------------------------------|----------------------|--------------|--------------|
| barcode01      | Sangopore_Prep_and_Run_Dev/Test_Barcoded_Fastq_Inputs/barcode01        | Group1     | 100          | 550          |
| barcode03      | Sangopore_Prep_and_Run_Dev/Test_Barcoded_Fastq_Inputs/barcode03        | Group1     | 100          | 550          |
| barcode09      | Sangopore_Prep_and_Run_Dev/Test_Barcoded_Fastq_Inputs/barcode09        | Group1     | 500          | 1000         |

<small>Note: In the `input_barcode_file` column absolute addresses should be used. The above example is for test code only.</small> 

## Table creation, either:
**Manual** - Use the above table example as a template (column headers must be exactly the same)
**Partially auto-filled** - Use one of the two examples below (only the headers and first two columns will be auto-filled).

1. If all fastq files are in one directory with one fastq file per barcode, as below:
   ```bash
   Parent_Directory
   ├── barcode01.fastq.gz
   ├── barcode03.fastq.gz
   └── barcode09.fastq.gz
   ```
   Run the `sango-create-table-onedir-allBCs` command
   ```
   # TEST EXAMPLE
   sango-create-table-onedir-allBCs \
      --parent_input_dir "Sangopore_Prep_and_Run_Dev/Test_OneDir_BC_per_fastq" \
      --output_dir "Sangopore_Prep_and_Run_Dev/Test_Outputs" \
      --output_tsv_name "AllInOneDir_table.tsv"
   # Use full paths for "--parent_input_dir" and "--output_dir" for user data
   ```

1. If there are multiple folders, each folder representing one barcode, and possibly many fastq files per barcode, as below:

   ```bash
   Parent_Directory
   ├── barcode01
   │   ├── barcode01_0.fastq.gz
   │   ├── barcode01_1.fastq.gz
   │   └── barcode01_2.fastq.gz
   ├── barcode03
   │   ├── barcode03_0.fastq.gz
   │   └── barcode03_1.fastq.gz
   └── barcode09
       ├── barcode03_0.fastq.gz
       └── barcode03_1.fastq.gz
   ```

   ```
   # TEST EXAMPLE
   sango-create-table-multidir-multifiles \
      --parentInputDir "Sangopore_Prep_and_Run_Dev/Test_Barcoded_Fastq_Inputs" \
      --outputDir "Sangopore_Prep_and_Run_Dev/Test_Outputs" \
      --outTableName "MultiDirMultiFiles_table.tsv"
   # Use full paths for "--parentInputDir" and "--outputDir" for user data
   ```

   <small>**Note:** If running the above test commands the user must be in the parent directory (default = Sangopore_Prep_and_Run_Dev).</small>

### Summary notes
   - If using either of the table-creation scripts, the remaining three columns still need to be filled in ("subgroup_folder", "min_range_bp", "max_range_bp").
   - To group outputs into specific subdirectories name them within the "subgroup_folder" column.
     - For example if the "subgroup_folder" column has either "Run1_Group1" or "Run1_Group2" in each row, script outputs for each barcode will be into either a subdirectory named "Run1-Group1" or "Run1_Group2".
- The script assumes that barcode numbers in filenames are always in the format `barcodeXX`, where `XX` is a two-digit number (such as barcode01).
