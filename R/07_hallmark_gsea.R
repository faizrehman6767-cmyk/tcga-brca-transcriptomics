# Protein-coding Hallmark ranked GSEA --------------------------------------

library(DESeq2)
library(dplyr)
library(fgsea)
library(msigdbr)

if (!exists("project_folder")) {
    source("config/project_config.R")
}

gene_annotation <- readRDS(
    file.path(project_folder, "RNA_gene_annotations_filtered.rds")
)

hallmark_df <- msigdbr(
    species = "Homo sapiens",
    collection = "H"
)

hallmark_pathways <- split(
    hallmark_df$gene_symbol,
    hallmark_df$gs_name
)

make_protein_coding_rank <- function(result_object, annotation) {
    result_df <- as.data.frame(result_object) %>%
        mutate(Ensembl_ID = rownames(.)) %>%
        left_join(
            annotation %>%
                as.data.frame() %>%
                mutate(Ensembl_ID = rownames(.)),
            by = "Ensembl_ID"
        ) %>%
        filter(
            gene_type == "protein_coding",
            !is.na(gene_name),
            gene_name != "",
            !is.na(stat)
        ) %>%
        arrange(desc(abs(stat))) %>%
        distinct(gene_name, .keep_all = TRUE)

    ranks <- result_df$stat
    names(ranks) <- result_df$gene_name
    sort(ranks, decreasing = TRUE)
}

run_hallmark_gsea <- function(ranks, comparison) {
    fgseaMultilevel(
        pathways = hallmark_pathways,
        stats = ranks,
        minSize = 15,
        maxSize = 500,
        eps = 0
    ) %>%
        as.data.frame() %>%
        mutate(
            Comparison = comparison,
            Enriched_toward = if_else(NES > 0, "numerator", "denominator")
        )
}

# Apply make_protein_coding_rank() and run_hallmark_gsea() to each DESeq2
# contrast from the adjusted and matched unadjusted common-support models.

