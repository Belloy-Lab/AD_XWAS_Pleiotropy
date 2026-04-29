## Coloc ABF base script (Applicable for both trait pleiotropy and xQTL COLOC analyses)

# Libraries
library(coloc)
library(tidyverse)
library(locuscomparer)


###########################################################################################
## Prepare data for coloc.abf

# D1 (AD XWAS data from merged input file)
D1 <- list(
  pvalues = merged$AD_P, # AD GWAS p-value
  type = "cc", # "cc" for case-control outcome variable
  snp = merged$SNP, # SNP
  position = merged$BP, # SNP basepair position
  N = merged$AD_N, # AD GWAS sample size
  MAF = merged$AD_EAF, # SNP effect allele frequency in AD GWAS
  s = merged_c$s # Number of cases / (number of cases + number of controls)
)

check_dataset(D1)

# D2 (QTL data from merged input file)
D2 <- list(
  pvalues = merged$QTL_P, # QTL p-value
  type = "quant", # "quant" for quantitative outcome variable
  snp = merged$SNP, # SNP
  position = merged$BP, # SNP basepair position
  N = merged$QTL_N, # QTL data sample size
  MAF = merged$QTL_EAF # SNP effect allele frequency in QTL data
)

check_dataset(D2)


###########################################################################################
## Run Colocalization

# Run coloc.abf (default parameters)
coloc_results_default = coloc.abf(dataset1 = D1, dataset2 = D2)


###########################################################################################
## Create Locus Compare plot

# Get variables for outcome variable 1 (Alzheimer's disease) - extract rsID and p-value
ad_comp <- merged %>%
  dplyr::select(SNP, AD_P) %>%
  dplyr::rename(rsid = SNP, pval = AD_P) %>% 
  arrange(pval)

# Get variables for outcome variable 2 (QTL) - extract rsID and p-value
qtl_comp <- merged %>%
  dplyr::select(SNP, QTL_P) %>%
  dplyr::rename(rsid = SNP, pval = QTL_P) %>% 
  arrange(pval)

# Define lead SNP (Top SNP in AD dataset)
rsID1 = ad_comp$rsid[1]

# Create locus compare plot
locuscompare(in_fn1 = ad_comp, in_fn2 = qtl_comp, snp = rsID1, title = "AD GWAS", title2 = "QTL", genome = "hg38", lz_ylab_linebreak = T)

