# Query and audit TCGA-BRCA RNA and mutation metadata ----------------------

library(TCGAbiolinks)
library(dplyr)
library(stringr)

if (!exists("project_folder")) {
    source("config/project_config.R")
}

query_rna <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
)

rna_metadata <- getResults(query_rna)

rna_primary <- rna_metadata %>%
    filter(sample_type == "Primary Tumor") %>%
    mutate(
        Sample_ID = sample.submitter_id,
        Case_ID = cases.submitter_id
    )

rna_summary <- rna_primary %>%
    summarise(
        RNA_file_rows = n(),
        unique_samples = n_distinct(Sample_ID),
        unique_cases = n_distinct(Case_ID)
    )

query_maf <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open"
)

maf_metadata <- getResults(query_maf)

maf_metadata_cases <- maf_metadata %>%
    transmute(
        Case_ID = str_sub(cases, 1, 12)
    ) %>%
    distinct()

print(rna_summary)
cat("RNA query rows:", nrow(rna_metadata), "\n")
cat("MAF query rows:", nrow(maf_metadata), "\n")
cat("MAF metadata cases:", nrow(maf_metadata_cases), "\n")

saveRDS(query_rna, file.path(project_folder, "query_RNA_all_STAR_counts.rds"))
saveRDS(rna_primary, file.path(project_folder, "RNA_primary_metadata.rds"))
saveRDS(query_maf, file.path(project_folder, "query_MAF_masked_somatic.rds"))
saveRDS(maf_metadata, file.path(project_folder, "MAF_metadata.rds"))

