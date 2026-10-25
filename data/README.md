# Data directory

Raw and processed TCGA data are not committed to this repository.

The workflow queries TCGA-BRCA through TCGAbiolinks and uses:

- STAR - Counts gene-expression quantification;
- Primary Tumor RNA samples;
- open-access Masked Somatic Mutation MAF files; and
- TCGA molecular annotations, including published PAM50 subtype labels.

Expected local-only directories and objects are excluded through `.gitignore`.
Users must comply with GDC data-use requirements and should not commit
controlled-access data.

