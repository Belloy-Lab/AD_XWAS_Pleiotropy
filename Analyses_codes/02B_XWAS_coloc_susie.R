## Coloc SuSiE base script (Applicable for both trait pleiotropy and xQTL COLOC analyses)

# Libraries
library(coloc)
library(tidyverse)
library(locuscomparer)


###########################################################################################
## Prepare data for coloc.susie

# AD GWAS coloc.susie preperation
ADl <- as.list(merged[,c("AD_BETA","varbeta_ad","SNP","BP","AD_EAF")])
names(ADl)[1:5] <- c("beta","varbeta","snp","position","MAF")
names(ADl$beta) <- ADl$snp
names(ADl$varbeta) <- ADl$snp
ADl$N <- median(merged$AD_N, na.rm = TRUE)

matching_snps_ADl <- ADl$snp %in% snps
ADl$beta <- ADl$beta[matching_snps_ADl]
ADl$varbeta <- ADl$varbeta[matching_snps_ADl]
ADl$snp <- ADl$snp[matching_snps_ADl]
ADl$position <- ADl$position[matching_snps_ADl]
ADl$MAF <- ADl$MAF[matching_snps_ADl]
ADl$LD <- LD
ADl$type <- "cc"
cc <- ADl

common_snps <- intersect(cc$snp, colnames(cc$LD))
cc_index <- match(common_snps, cc$snp)
ld_index <- match(common_snps, colnames(cc$LD))

cc$beta <- cc$beta[cc_index]
cc$varbeta <- cc$varbeta[cc_index]
cc$snp <- cc$snp[cc_index]
cc$position <- cc$position[cc_index]
cc$MAF <- cc$MAF[cc_index]
cc$LD <- cc$LD[ld_index, ld_index]


# QTL data coloc.susie preperation
traitl <- as.list(merged[,c("QTL_BETA","varbeta_t2","SNP","BP","QTL_EAF")])
names(traitl)[1:5] <- c("beta","varbeta","snp","position","MAF")
names(traitl$beta) <- traitl$snp
names(traitl$varbeta) <- traitl$snp
traitl$N <- round(median(merged$QTL_N, na.rm = TRUE))

matching_snps_traitl <- traitl$snp %in% snps
traitl$beta <- traitl$beta[matching_snps_traitl]
traitl$varbeta <- traitl$varbeta[matching_snps_traitl]
traitl$snp <- traitl$snp[matching_snps_traitl]
traitl$position <- traitl$position[matching_snps_traitl]
traitl$MAF <- traitl$MAF[matching_snps_traitl]
traitl$LD <- LD
traitl$type <- "quant"
eq <- traitl

common_snps <- intersect(eq$snp, colnames(eq$LD))
eq_index <- match(common_snps, eq$snp)
ld_index <- match(common_snps, colnames(eq$LD))

eq$beta <- eq$beta[eq_index]
eq$varbeta <- eq$varbeta[eq_index]
eq$snp <- eq$snp[eq_index]
eq$position <- eq$position[eq_index]
eq$MAF <- eq$MAF[eq_index]
eq$LD <- eq$LD[ld_index, ld_index]


###########################################################################################
## Run coloc.susie with sequenced coverages
coverage_values <- seq(0.95, 0.05, by = -0.1)
cs_length_S1 <- 0
cs_length_S2 <- 0
for (coverage in coverage_values) {
  tryCatch({
    S1 <- runsusie(cc, coverage = coverage, max_iter = 10000)
    cs_length_S1 <- length(S1$sets$cs)
    
    S2 <- runsusie(eq, coverage = coverage, max_iter = 10000)
    cs_length_S2 <- length(S2$sets$cs)
    
    if (cs_length_S1 > 0 & cs_length_S2 > 0) {
      break
    }
  }, error = function(e) {
    # Print the error message if you want to
    message("Error encountered: ", e)
    # Continue to the next iteration
  })
}


# Run coloc.susie and return the maximum posterior probability
out <- coloc.susie(S1, S2)


###########################################################################################
## Get results for maximum PP4
print(out$summary)
index_row <- which.max(out$summary$PP.H4.abf)
coloc_PP0 = out$summary$PP.H0.abf[index_row]
coloc_PP1 = out$summary$PP.H1.abf[index_row]
coloc_PP2 = out$summary$PP.H2.abf[index_row]
coloc_PP3 = out$summary$PP.H3.abf[index_row]
coloc_PP4 = out$summary$PP.H4.abf[index_row]

coloc_hit1 = out$summary$hit1[index_row]
coloc_hit2 = out$summary$hit2[index_row]



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