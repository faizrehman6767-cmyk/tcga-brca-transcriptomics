# Prepare one RNA count column per case ------------------------------------

library(SummarizedExperiment)
library(dplyr)
library(stringr)

if (!exists("project_folder")) {
    source("config/project_config.R")
}

# Expected input: rna_se returned by GDCprepare(..., summarizedExperiment=TRUE)
stopifnot(exists("rna_se"))

raw_counts <- assay(rna_se, "unstranded")
rna_coldata <- as.data.frame(colData(rna_se))

rna_coldata <- rna_coldata %>%
    mutate(
        RNA_column = colnames(raw_counts),
        Sample_ID = str_sub(sample_submitter_id, 1, 16),
        Case_ID = str_sub(sample_submitter_id, 1, 12),
        column_number = seq_len(n())
    )

# Sum raw counts from technical aliquots belonging to one biological sample.
counts_by_sample <- t(rowsum(
    t(raw_counts),
    group = rna_coldata$Sample_ID,
    reorder = FALSE
))

sample_metadata <- rna_coldata %>%
    group_by(Sample_ID) %>%
    summarise(
        across(everything(), first),
        number_of_aliquots = n(),
        .groups = "drop"
    )

sample_metadata$collapsed_total_counts <- colSums(counts_by_sample)
sample_metadata$collapsed_detected_genes <- colSums(counts_by_sample > 0)

# Select one primary-tumor sample per case. Non-FFPE material is prioritized.
# Total counts are only a final tie-breaker after biological/technical criteria.
selected_sample_metadata <- sample_metadata %>%
    mutate(
        ffpe_priority = if_else(
            preservation_method == "FFPE",
            1L,
            0L,
            missing = 0L
        )
    ) %>%
    arrange(
        Case_ID,
        ffpe_priority,
        days_to_collection,
        desc(collapsed_total_counts),
        Sample_ID
    ) %>%
    group_by(Case_ID) %>%
    slice(1L) %>%
    ungroup()

selected_sample_ids <- selected_sample_metadata$Sample_ID
counts_by_case <- counts_by_sample[, selected_sample_ids, drop = FALSE]
colnames(counts_by_case) <- selected_sample_metadata$Case_ID
rownames(selected_sample_metadata) <- selected_sample_metadata$Case_ID

stopifnot(
    ncol(counts_by_case) == nrow(selected_sample_metadata),
    identical(colnames(counts_by_case), rownames(selected_sample_metadata))
)

saveRDS(counts_by_case, file.path(project_folder, "RNA_counts_by_case.rds"), compress = FALSE)
saveRDS(selected_sample_metadata, file.path(project_folder, "RNA_selected_sample_metadata.rds"))

