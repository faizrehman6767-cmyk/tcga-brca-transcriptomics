# Package installation and loading -----------------------------------------

cran_packages <- c(
    "dplyr", "tidyr", "stringr", "tibble", "ggplot2", "ggrepel",
    "msigdbr"
)

bioconductor_packages <- c(
    "TCGAbiolinks", "SummarizedExperiment", "DESeq2", "apeglm", "fgsea"
)

missing_cran <- cran_packages[
    !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0L) {
    install.packages(missing_cran)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

missing_bioc <- bioconductor_packages[
    !vapply(bioconductor_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_bioc) > 0L) {
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

invisible(lapply(
    c(cran_packages, bioconductor_packages),
    library,
    character.only = TRUE
))

cat("Required packages loaded.\n")

