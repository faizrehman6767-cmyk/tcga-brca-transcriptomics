# Expanded pathogenic-variant sensitivity grouping ------------------------

library(dplyr)

if (!exists("project_folder")) {
    source("config/project_config.R")
}

annotation_final <- readRDS(
    file.path(project_folder, "BRCA1_BRCA2_variant_annotation_final.rds")
)

strict_case_groups <- readRDS(
    file.path(project_folder, "strict_mutation_groups.rds")
)

expanded_additions <- annotation_final %>%
    filter(Include_in_expanded_group == "Yes") %>%
    transmute(
        Case_ID,
        Expanded_Group = case_when(
            Hugo_Symbol == "BRCA1" ~ "BRCA1-mutant",
            Hugo_Symbol == "BRCA2" ~ "BRCA2-mutant",
            TRUE ~ NA_character_
        ),
        Added_Variant = paste(Hugo_Symbol, HGVSc, HGVSp_Short)
    )

expanded_case_groups <- strict_case_groups %>%
    rename(Strict_Mutation_Group = Mutation_Group) %>%
    left_join(expanded_additions, by = "Case_ID") %>%
    mutate(
        Mutation_Group = if_else(
            !is.na(Expanded_Group) &
                Strict_Mutation_Group == "BRCA1/2-WT",
            Expanded_Group,
            Strict_Mutation_Group
        ),
        Grouping_Definition = if_else(
            !is.na(Expanded_Group) &
                Strict_Mutation_Group == "BRCA1/2-WT",
            "Added by pathogenic missense annotation",
            "Retained from strict classification"
        )
    )

strict_counts <- count(strict_case_groups, Mutation_Group, name = "Strict_n")
expanded_counts <- count(expanded_case_groups, Mutation_Group, name = "Expanded_n")

group_count_comparison <- full_join(
    strict_counts,
    expanded_counts,
    by = "Mutation_Group"
) %>%
    mutate(Change = Expanded_n - Strict_n)

print(group_count_comparison)

stopifnot(
    nrow(expanded_case_groups) == nrow(strict_case_groups),
    n_distinct(expanded_case_groups$Case_ID) == nrow(expanded_case_groups)
)

saveRDS(
    expanded_case_groups,
    file.path(project_folder, "expanded_pathogenic_mutation_groups.rds")
)

