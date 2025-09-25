## Detailed Script Workflow

1. **Output Directory Creation:**
    - The script checks if the specified output directory exists. If not, it creates the directory.

1. **Searching for `.fastq.gz` Files:**
    - The script recursively traverses the parent input directory to find either all barcode sub-folders and their `.fastq.gz` or the single read folder and their each `.fastq.gz` files representing a barcode.

1. **Barcode Extraction:**
   - For each barcode folder or `.fastq.gz` file, the script uses a regular expression to extract the barcode number from the barode folder or barcode filename (e.g., `barcode12`).

1. **Table Creation (5 columns)**
   - Auto filled in by script:
     - `barcode number`
     - `input directory path`
   - Filled in by user, after script run:
     - `subgroup_folder`
     - `min_range_bp`
     - `max_range_bp`

1. **Data Collection:**
   - The script collects the `barcode number`, as well as the `input directory path` containing the `.fastq.gz` file.
   - Placeholder columns for `subgroup_folder`, `min_range_bp`, and `max_range_bp` are included.

1. **DataFrame Creation and Processing:**
   - The script creates a pandas DataFrame from the collected data.
   - It removes duplicate entries based on the `barcode_number` and `input_barcode_file` columns.
   - The DataFrame is sorted based on the numeric part of the barcode number.

1. **Saving the TSV File:**
   - The processed DataFrame is saved as a TSV file in the specified output directory with the given file name.

