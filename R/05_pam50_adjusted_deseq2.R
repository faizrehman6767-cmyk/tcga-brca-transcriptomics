# PAM50 common-support DESeq2 models ---------------------------------------

library(DESeq2)
library(dplyr)

if (!exists("project_folder")) {
    source("config/project_config.R")
}

dds_primary <- readRDS(
    file.path(project_folder, "DESeq2_primary_model_966_cases.rds")
)

analysis_metadata <- as.data.frame(colData(dds_primary))

# BRCA1/2 mutant cases were represented within Basal, LumA, and LumB.
# Restricting to these levels avoids extrapolation across unsupported subtypes.
common_subtypes <- c("LumA", "LumB", "Basal")

keep_common <- analysis_metadata$paper_BRCA_Subtype_PAM50 %in% common_subtypes
dds_common <- dds_primary[, keep_common]

colData(dds_common)$PAM50_model <- factor(
    colData(dds_common)$paper_BRCA_Subtype_PAM50,
    levels = c("LumA", "LumB", "Basal")
)

design(dds_common) <- ~ PAM50_model + Mutation_Group_model
dds_common_adjusted <- DESeq(dds_common)

dds_common_unadjusted <- dds_common
design(dds_common_unadjusted) <- ~ Mutation_Group_model
dds_common_unadjusted <- DESeq(dds_common_unadjusted)

saveRDS(
    dds_common_adjusted,
    file.path(project_folder, "DESeq2_PAM50_adjusted_common_support_845_cases.rds"),
    compress = FALSE
)

saveRDS(
    dds_common_unadjusted,
    file.path(project_folder, "DESeq2_unadjusted_common_support_845_cases.rds"),
    compress = FALSE
)

cat("Adjusted coefficients:\n")
print(resultsNames(dds_common_adjusted))

