# Construct strict BRCA1/2 mutation groups ---------------------------------

library(dplyr)
library(stringr)

if (!exists("project_folder")) {
    source("config/project_config.R")
}

# Expected input: GDCprepare(query_maf, directory = maf_download_folder)
# The large MAF object is intentionally not version-controlled.

stopifnot(exists("maf"))

brca_variants <- maf %>%
    filter(
        Hugo_Symbol %in% c("BRCA1", "BRCA2"),
        Mutation_Status == "Somatic",
        is.na(GDC_FILTER)
    ) %>%
    mutate(
        Case_ID = str_sub(Tumor_Sample_Barcode, 1, 12)
    )

strict_classes <- c(
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Nonsense_Mutation",
    "Splice_Site"
)

strict_brca_variants <- brca_variants %>%
    filter(Variant_Classification %in% strict_classes)

strict_mutant_records <- strict_brca_variants %>%
    distinct(Case_ID, Hugo_Symbol) %>%
    mutate(
        Mutation_Group = case_when(
            Hugo_Symbol == "BRCA1" ~ "BRCA1-mutant",
            Hugo_Symbol == "BRCA2" ~ "BRCA2-mutant"
        )
    ) %>%
    select(Case_ID, Mutation_Group)

double_group_cases <- strict_mutant_records %>%
    count(Case_ID, name = "groups_per_case") %>%
    filter(groups_per_case > 1L)

stopifnot(nrow(double_group_cases) == 0L)

all_maf_cases <- maf %>%
    transmute(Case_ID = str_sub(Tumor_Sample_Barcode, 1, 12)) %>%
    distinct()

strict_case_groups <- all_maf_cases %>%
    left_join(strict_mutant_records, by = "Case_ID") %>%
    mutate(
        Mutation_Group = coalesce(Mutation_Group, "BRCA1/2-WT")
    )

print(count(strict_case_groups, Mutation_Group))

saveRDS(
    strict_case_groups,
    file.path(project_folder, "strict_mutation_groups.rds")
)

write.csv(
    strict_brca_variants,
    file.path(project_folder, "strict_brca_variants.csv"),
    row.names = FALSE
)

