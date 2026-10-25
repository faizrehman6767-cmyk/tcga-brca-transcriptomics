# Variance-stabilizing transformation and PCA ------------------------------

library(DESeq2)
library(dplyr)
library(ggplot2)

if (!exists("project_folder")) {
    source("config/project_config.R")
}

dds_primary <- readRDS(
    file.path(project_folder, "DESeq2_primary_model_966_cases.rds")
)

vsd <- vst(dds_primary, blind = FALSE)
pca <- prcomp(t(assay(vsd)))
percent_variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_data <- cbind(
    as.data.frame(pca$x[, 1:10, drop = FALSE]),
    as.data.frame(colData(dds_primary))
)

write.csv(
    pca_data,
    file.path(project_folder, "PCA_sample_coordinates.csv"),
    row.names = FALSE
)

write.csv(
    data.frame(
        PC = paste0("PC", seq_along(percent_variance)),
        percent_variance = percent_variance
    ),
    file.path(project_folder, "PCA_percent_variance.csv"),
    row.names = FALSE
)

