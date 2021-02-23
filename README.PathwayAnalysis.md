# Pathway analysis based on GWAS data using Polygenic Risk Score



## Introduction
Despite the considerable progress made in unraveling the genetic causes of amyotrophic lateral sclerosis (ALS), we do not fully understand the molecular mechanisms underlying the disease. We analyzed genome-wide data involving nearly 80,000 individuals using a polygenic risk score approach to identify the biological pathways and cell types involved in ALS.
## Overview
We have performed a **three-stage study design** to identify pathways and cell types relevant to ALS risk

**Reference Dataset.** A published GWAS study involving 12,577 ALS cases and 23,475 controls (*Van Rheenen et al., 2016*)
https://www.ncbi.nlm.nih.gov/pubmed/27455348


**Training Dataset**. 	Plink Binary files containing individual-level data  genotype consisting in 5,605 ALS cases and 24,110 control subjects *(Nicolas et al., 2018)*.
https://www.ncbi.nlm.nih.gov/pubmed/29566793

**Replication Dataset.** Plink Binary files containing individual-level data consisting in 2,411 ALS cases and 10,322 controls *(Nicolas et al., 2018)*.
https://www.ncbi.nlm.nih.gov/pubmed/29566793


# Pathway analysis based on Polygenic Risk Scores using PRSice2

A detailed information about PRSice2 can be found in the following link:
https://choishingwan.github.io/PRSice/

###### Citation: Choi SW, and O’Reilly PF. "PRSice-2: Polygenic Risk Score Software for Biobank-Scale Data." GigaScience 8, no. 7 (July 1, 2019). [https://doi.org/10.1093/gigascience/giz082](https://doi.org/10.1093/gigascience/giz082).
## Requirements

 - [ ] PRSice2_Linux executable (v2.1.1)
 - [ ] PLINK v1.90b4.4 64-bit   (www.cog-genomics.org/plink/1.9/)
 - [ ] PLINK binary files with PCAs and a covariate file. Due to the nature of PRS analysis, it is recommendable to include 20 PCs plus sex plus age.
 - [ ] [R](https://www.r-project.org/) (**version 3.2.3+**)

#### PRSet Specific Input

- [ ] **MSigDB file**: File containing name of each gene sets and the ID of genes within the gene set on each individual line. If MSigDB is provided, GTF file is required.
- [ ] **GTF file**: A file contain the genome boundary of each individual gene


## Base/Reference file
Summary Stats without overlapping samples with training dataset. We have used a previous published GWAS summary Stats (Van Rheenen et al., 2016) https://www.ncbi.nlm.nih.gov/pubmed/27455348

## Training file
**Training dataset**: Raw genotype data of "target phenotype" in the form of PLINK binary
https://www.ncbi.nlm.nih.gov/pubmed/29566793

**Testing/Replication dataset**: Raw genotype data of "target phenotype" in the form of PLINK binary
https://www.ncbi.nlm.nih.gov/pubmed/29566793
