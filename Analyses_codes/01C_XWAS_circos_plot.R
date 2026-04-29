## Create circular Manhattan Plot for 4 AD XWAS samples

## libraries
library(data.table)
library(tidyverse)
library(stringr)
library(CMplot)

###################################
# Read in AD XWAS summary stats (hg38): Variables should be: SNP, Chromosome, Basepair Position, P-value

rxci = fread("/path/to/rxci/summary/stats") %>% 
  mutate(CHR = "X") %>% 
  dplyr::rename(Chromosome = CHR, Position = BP, trait1 = P)

exci = fread("/path/to/exci/summary/stats") %>% 
  mutate(CHR = "X") %>% 
  dplyr::rename(Chromosome = CHR, Position = BP, trait2 = P)

female = fread("/path/to/female/summary/stats") %>% 
  mutate(CHR = "X") %>% 
  dplyr::rename(Chromosome = CHR, Position = BP, trait3 = P)

male = fread("/path/to/male/summary/stats") %>% 
  mutate(CHR = "X") %>% 
  dplyr::rename(Chromosome = CHR, Position = BP, trait4 = P)


## Merge 4 summary stat files to contain SNP, Chromosome, BP, trait1, .., traiti
df1 = merge(rxci, aexci, by = "SNP") %>% 
  dplyr::select(SNP, Chromosome.x, Position.x, trait1, trait2) %>% 
  dplyr::rename(Chromosome = Chromosome.x, Position = Position.x)


df2 = merge(df1, female, by = "SNP") %>% 
  dplyr::select(SNP, Chromosome.x, Position.x, trait1, trait2, trait3) %>% 
  dplyr::rename(Chromosome = Chromosome.x, Position = Position.x)


df3 = merge(df2, male, by = "SNP") %>% 
  dplyr::select(SNP, Chromosome.x, Position.x, trait1, trait2, trait3, trait4) %>% 
  dplyr::rename(Chromosome = Chromosome.x, Position = Position.x) 


# Create Circos plot using CMplot function from CMplot package - file output is saved to working directory
ad = CMplot(df3, 
            type = "p", 
            plot.type = "c", 
            multraits = T,
            r = 0.4, 
            H = 2,
            col = matrix(c("gray20", NA,  NA, NA,
                           NA, "gray70", NA, NA,
                           NA, NA, "darkolivegreen4", NA, 
                           NA, NA, NA, "mediumpurple1"), nrow=4, byrow=TRUE),
            chr.labels = paste("Chr", "X", sep = ""), 
            threshold=c(1e-5),
            cir.chr.h = 1.3,
            ylim = c(-log10(1), -log10(1e-8)),
            threshold.col=c("red"),
            signal.col=c("red"),
            cir.band = 2,
            cir.axis.grid = F,
            amplify = T,
            signal.cex = 0.85,
            signal.line = 2,
            bin.size = 1e6,
            main = c("Alzheimer's Disease XWAS"),
            mar.between = 5,
            legend.ncol = 1, 
            legend.pos = "middle",
            outward = T, 
            file = "jpg", 
            file.name = "plot_name", 
            dpi = 320, 
            file.output = TRUE, 
            verbose = TRUE,
            width = 11, 
            height = 11)







