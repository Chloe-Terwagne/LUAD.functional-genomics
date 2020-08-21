# LUAD.functional-genomics
## BINF-F401 : Computational Methods for Functional Genomics
### Chloé Terwagne
#### Registration number: 000409683

This repository contains the R code supporting the report for the BINF-F401 course.
The report present results of the study of the lung adenocarcinoma cancer associated with the estimated DNA methylation age using the publicly available data from The Cancer Genome Atlas (TCGA) project.

The files used for this project are download via Firehose (http://gdac.broadinstitute.org/) using TCGA data version 01.28.2016

- LUAD cancer: 'Clinical\_pick\_Tier1', Patients clinical annotation (including chronological age, demographic information, treatment information, survival data, etc).
- LUAD cancer: 'illuminahiseq\_rnaseqv2\_RSEM\_genes\_normalized', mRNA genes expression estimated from TCGA mRNA-seq data.
- LUAD cancer: 'Mutation\_Packager\_Calls', Mutation Annotation Format (MAF).
- LUAD cancer: 'genome\_wide\_snp\_6\_segmented\_scna\_minus\_germline\_cnv\_hg19', file containing the copy number variation (CNV).
- LUAD-7.Rda, contains the estimated DNAm age computed by Vincent Detours.

