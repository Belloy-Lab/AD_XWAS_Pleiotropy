## Run sex- and APOE4-stratified AD XWAS Meta-analyses using GWAMA (rXCI, Female, Male)

## Libraries
library(data.table)
library(tidyverse)
library(devtools)
library(Rfast)

## Establish output directory as working directory
out_dir = "/path/to/output/directory/"
setwd(out_dir)


################################################################################
## Create custom map and key file using input summary stats from contributing cohorts
# possibel cohorts: adgc, adsp, ukb, fg, mvp1, mvp2
map = rbind(adgc, adsp, ukb, fg, mvp1, mvp2) %>% distinct(MARKER, .keep_all = T) %>% arrange(CHR, POS) 

key = map %>% dplyr::select(MARKER, SNP)
map = map %>% dplyr::select(CHR, MARKER, EMP, POS)


## Create map and key file names that will be referenced in GWAMA and write out these files
mapn <- paste0(out_dir, "gwama.map.txt")
fwrite(map, mapn, col.names = F, sep="\t")

keyn <- paste(out_dir, "gwama.key.txt", sep="")
fwrite(key, keyn, col.names = T, sep="\t")


################################################################################
## Run GWAMA
# list of files to meta-analyze - filepath is a list of filepaths to the corresponding summary files for GWAMA
gwamaf = paste(filepath, collapse = "\n")
gwaman = paste(out_dir, "gwama.merge_list.txt", sep="")
writeLines(gwamaf, gwaman)


# Execute this in interactive session command line - Use dmr07083/general-utility:1.0 docker
# Be sure to chanmge directories to match to output directory
/usr/bin/GWAMA \
--map "/path/to/output/directory/gwama.map.txt" \
--name_marker MARKER \
--indel_alleles \
--filelist "/path/to/output/directory/gwama.merge_list.txt" \
--output "/path/to/output/name_of_GWAMA_output"




