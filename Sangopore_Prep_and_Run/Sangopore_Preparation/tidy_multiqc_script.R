#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments
if (length(args) != 2) {
    stop("Usage: Rscript tidy_multiqc_script.R <json_file_path> <output_csv_path>")
}

json_file_path <- args[1]
output_csv_path <- args[2]

# Check if the JSON file exists
if (!file.exists(json_file_path)) {
    stop(paste("Error: JSON file not found at", json_file_path))
}

# Load required libraries
if (!requireNamespace("TidyMultiqc", quietly = TRUE)) {
    install.packages("TidyMultiqc")
}
library(TidyMultiqc)

# Load the JSON file
tryCatch({
    df <- TidyMultiqc::load_multiqc(json_file_path)
    df_both <- TidyMultiqc::load_multiqc(json_file_path, sections = c('raw'))
    
    # Write the CSV file
    write.csv(df, output_csv_path, row.names = FALSE)
    cat("TidyMultiqc processing completed successfully.\n")
}, error = function(e) {
    cat("Error:", e$message, "\n")
    quit(status = 1)
})
