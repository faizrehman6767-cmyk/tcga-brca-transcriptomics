# TCGA-BRCA BRCA1/BRCA2 transcriptomics

An R-based, hypothesis-driven analysis of whether primary breast tumors with
somatic **BRCA1** and **BRCA2** loss-of-function variants show distinct
transcriptional profiles despite disruption of the same homologous
recombination repair pathway.

## Study questions

1. How does gene expression differ between BRCA1-mutant and BRCA1/2
   mutation-negative tumors?
2. How does gene expression differ between BRCA2-mutant and BRCA1/2
   mutation-negative tumors?
3. How do BRCA1-mutant and BRCA2-mutant tumors differ from each other?
4. Which findings persist after accounting for PAM50 intrinsic subtype?

The label `BRCA1/2-WT` is used as a short analysis label. More precisely, it
means **no qualifying BRCA1 or BRCA2 somatic coding mutation was detected under
the study definition**. It does not exclude germline variants, copy-number
loss, structural variants, promoter methylation, or other mechanisms of
homologous-recombination deficiency.

## Data source and cohort construction

- Project: TCGA-BRCA
- Reference genome: GRCh38/hg38
- RNA-seq workflow: STAR - Counts
- Mutation data: open-access Masked Somatic Mutation MAF
- RNA sample type: Primary Tumor

Metadata inspection identified:

| Stage | Number |
|---|---:|
| STAR-count RNA files | 1,231 |
| Primary-tumor RNA files | 1,111 |
| Unique primary-tumor sample IDs | 1,106 |
| Unique primary-tumor RNA cases | 1,095 |
| MAF files | 992 |
| Cases represented by MAF metadata | 969 |
| RNA and mutation intersection | 966 |

Technical RNA aliquots from the same biological sample were summed at the raw
count level. When a case contained more than one primary-tumor sample, sample
selection prioritized non-FFPE material and preserved one biological sample per
case. The final RNA matrix contained 60,660 genes across 966 cases before gene
filtering. Genes with at least 10 counts in at least 8 tumors were retained,
leaving 35,720 genes.

## Mutation definitions

### Primary, strict definition

The primary mutant groups include somatic exonic variants classified as:

- frameshift deletion;
- frameshift insertion;
- nonsense mutation; or
- qualifying splice-site mutation.

The RNA-linked primary cohort contained:

| Group | Cases |
|---|---:|
| BRCA1-mutant | 9 |
| BRCA2-mutant | 8 |
| BRCA1/2 mutation-negative | 949 |

### Expanded sensitivity definition

Thirty excluded BRCA1/2 missense or in-frame variants were reviewed using
ClinVar, BRCA Exchange/ENIGMA, OncoKB evidence when available, and the existing
VEP annotations. Computational predictions alone were not considered
sufficient for inclusion.

One additional variant met the high-confidence inclusion rule:

- BRCA1 c.5363G>T, p.Gly1788Val in TCGA-AN-A0XU

The expanded mutation grouping has been constructed and validated. The
RNA-linked DESeq2 and GSEA sensitivity rerun is the next planned analysis and is
not presented as completed in this repository version.

## Completed analyses

- GDC metadata querying and cohort audit with TCGAbiolinks
- RNA/MAF case linkage
- technical-aliquot count aggregation
- selection of one primary-tumor sample per case
- DESeq2 primary differential-expression models
- PAM50-adjusted DESeq2 models on common subtype support
- matched unadjusted-versus-adjusted comparisons
- variance-stabilizing transformation and PCA
- protein-coding Hallmark ranked GSEA with `fgsea`
- pathway stability and leading-edge analyses
- manual review of excluded missense and in-frame variants

PAM50 labels were obtained from the TCGA molecular annotations added by
TCGAbiolinks. PAM50 subtypes were not independently recalculated. A PAM50
risk-of-recurrence score was not calculated.

## Main completed results

The primary 966-case model identified 947, 440, and 238 genes at FDR < 0.05
for BRCA1 vs WT, BRCA2 vs WT, and BRCA1 vs BRCA2, respectively. In the
PAM50 common-support cohort, the adjusted models identified 131, 77, and 68
genes. Full interpretation must account for the small mutant groups and the
strong association between PAM50 subtype and the leading expression PCs.

Summary tables are available in [`results/`](results/README.md).

## Repository structure

```text
.
├── R/                  Reproducible analysis scripts
├── config/             Local path configuration template
├── data/               Documentation only; raw data are not tracked
├── docs/               Study design and interpretation notes
├── figures/            Figure documentation/placeholders
├── results/            Small summary tables suitable for version control
├── .gitignore
├── CITATION.cff
└── README.md
```

## Reproducing the workflow

1. Install R and the packages listed in `R/00_setup.R`.
2. Copy `config/project_config.example.R` to
   `config/project_config.R` and set a local project directory.
3. Run the scripts in numerical order.
4. Do not commit GDC downloads or large fitted model objects.

The scripts are organized as a documented reconstruction of the completed
workflow. Before treating the repository as a frozen research compendium, the
scripts should be compared with the original RStudio history and saved local
objects.

## Data availability

TCGA data are obtained from the NCI Genomic Data Commons. Raw RNA-seq counts,
MAF files, count matrices, and fitted DESeq2 objects are intentionally excluded
from this repository because they are large and can be regenerated from GDC.

## Author

Faizur Rehman

