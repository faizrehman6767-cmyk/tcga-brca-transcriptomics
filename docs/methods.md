# Methods summary

## Cohort

TCGA-BRCA RNA-seq and somatic mutation metadata were queried through
TCGAbiolinks. Primary-tumor RNA files were retained. RNA and mutation records
were linked at the TCGA case level.

Technical RNA aliquots from the same biological sample were combined by
summing raw STAR counts. Separate primary-tumor samples from the same case were
not combined. One sample per case was selected using specimen-preservation and
collection information, with non-FFPE material prioritized.

## Mutation grouping

The primary definition classified frameshift, nonsense, and qualifying
splice-site variants in BRCA1 or BRCA2 as loss-of-function candidates. Silent
variants were excluded. Missense and in-frame variants were not included on the
basis of SIFT or PolyPhen predictions alone.

For sensitivity analysis, excluded missense and in-frame variants were reviewed
using ClinVar, BRCA Exchange/ENIGMA, and somatic oncogenicity evidence when
available. Only variants with high-confidence pathogenic/likely pathogenic or
oncogenic loss-of-function evidence were eligible for addition.

## RNA analysis

Genes were retained when they had at least 10 raw counts in at least 8 tumors.
DESeq2 was used for the three prespecified contrasts. Variance-stabilized counts
were used for PCA, not for differential-expression testing.

## PAM50 adjustment

PAM50 labels were taken from the molecular annotations provided with the TCGA
sample metadata. The adjusted analysis was restricted to PAM50 levels shared by
the mutant groups and controls (Basal, LumA, and LumB), producing an 845-case
common-support cohort. Models with and without PAM50 were fitted to the same
cohort to separate sample-restriction effects from covariate-adjustment effects.

## Pathway analysis

Protein-coding genes were ranked using the DESeq2 Wald statistic. Hallmark gene
sets were analyzed with fgsea. Adjusted and matched unadjusted normalized
enrichment scores were compared, and leading-edge genes were examined for
stable pathways.

## Limitations

- Mutant groups are small, particularly after PAM50 common-support restriction.
- PAM50 is strongly associated with the leading expression principal
  components and is also biologically related to BRCA mutation status.
- The mutation-negative label does not establish normal BRCA function.
- The MAF-based definition does not capture all germline, copy-number,
  structural, methylation, or regulatory mechanisms of BRCA loss.

