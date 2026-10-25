# Copy this file to config/project_config.R and edit the path.
# The copied file is ignored by Git so local paths are not published.

project_folder <- Sys.getenv("TCGA_BRCA_DIR")

if (!nzchar(project_folder)) {
    stop(
        paste(
            "Set the TCGA_BRCA_DIR environment variable or edit",
            "config/project_config.R with your local project directory."
        )
    )
}

rna_download_folder <- file.path(project_folder, "RNA_data")
maf_download_folder <- file.path(project_folder, "MAF_data")
results_folder <- file.path(project_folder, "results")

dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)

