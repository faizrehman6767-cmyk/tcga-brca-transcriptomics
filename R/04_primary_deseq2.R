# Primary DESeq2 analysis ---------------------------------------------------

library(DESeq2)
library(dplyr)

if (!exists("project_folder")) {
    source("config/project_config.R")
}

counts_by_case <- readRDS(file.path(project_folder, "RNA_counts_by_case.rds"))
analysis_metadata <- readRDS(file.path(project_folder, "DESeq2_analysis_metadata_966_cases.rds"))

analysis_metadata <- as.data.frame(analysis_metadata)
analysis_metadata <- analysis_metadata[colnames(counts_by_case), , drop = FALSE]

stopifnot(identical(colnames(counts_by_case), rownames(analysis_metadata)))

keep_gene <- rowSums(counts_by_case >= 10) >= 8
counts_filtered <- counts_by_case[keep_gene, , drop = FALSE]

analysis_metadata$Mutation_Group_model <- factor(
    analysis_metadata$Mutation_Group_model,
    levels = c("WT", "BRCA1_mutant", "BRCA2_mutant")
)

dds_primary <- DESeqDataSetFromMatrix(
    countData = round(counts_filtered),
    colData = analysis_metadata,
    design = ~ Mutation_Group_model
)

dds_primary <- DESeq(dds_primary)

res_BRCA1_vs_WT <- results(
    dds_primary,
    contrast = c("Mutation_Group_model", "BRCA1_mutant", "WT"),
    alpha = 0.05
)

res_BRCA2_vs_WT <- results(
    dds_primary,
    contrast = c("Mutation_Group_model", "BRCA2_mutant", "WT"),
    alpha = 0.05
)

res_BRCA1_vs_BRCA2 <- results(
    dds_primary,
    contrast = c("Mutation_Group_model", "BRCA1_mutant", "BRCA2_mutant"),
    alpha = 0.05
)

saveRDS(
    dds_primary,
    file.path(project_folder, "DESeq2_primary_model_966_cases.rds"),
    compress = FALSE
)

saveRDS(res_BRCA1_vs_WT, file.path(project_folder, "DESeq2_BRCA1_vs_WT_raw_results.rds"))
saveRDS(res_BRCA2_vs_WT, file.path(project_folder, "DESeq2_BRCA2_vs_WT_raw_results.rds"))
saveRDS(res_BRCA1_vs_BRCA2, file.path(project_folder, "DESeq2_BRCA1_vs_BRCA2_raw_results.rds"))

