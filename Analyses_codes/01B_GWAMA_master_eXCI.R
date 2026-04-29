## Example script of runnning eXCI XWAS with Male and Female XWAS summary stats

## Libraries
library(data.table)
library(tidyverse)
library(devtools)
library(Rfast)

## Establish output directory as working directory
out_dir = "/path/to/output/directory/"
setwd(out_dir)

########################################################################################################################################################################################################################
### GWAMA for eXCI model

# Read in Female Summary stats
f = fread("/path/to/clean/female/summary/stats.gen090")

# Read in Male summary stats (If males coded as 0/2, multiply BETA and SE by 2. If males coded as 0/1, no need to adjust. Check your coding before proceeding!)
m = fread("//path/to/clean/male/summary/stats.gen090")
m = m %>% 
  mutate(BETA = 2*BETA,
         SE = 2*SE)

#########################
## save individual sumstats - write out
fwrite(f, paste0(out_dir, "XWAS_female.txt"), col.names = T, row.names = F, quote = F, sep = '\t')
fwrite(m, paste0(out_dir, "XWAS_male.txt"), col.names = T, row.names = F, quote = F, sep = '\t')


################################################################################
## Create custom map and key file to be referenced later in script
map = rbind(f, m) %>% distinct(MARKER, .keep_all = T) %>% arrange(CHR, POS) 

key = map %>% dplyr::select(MARKER, SNP)
map = map %>% dplyr::select(CHR, MARKER, POS)


## write out these files
mapn <- paste0(out_dir, "gwama.map.txt")
fwrite(map, mapn, col.names = F, sep="\t")

keyn <- paste(out_dir, "gwama.key.txt", sep="")
fwrite(key, keyn, col.names = T, sep="\t")



################################################################################
## create gwama.in file that contains list of all study files each on a separate row
gwamad = as.data.frame(c(paste0(out_dir, "XWAS_female.txt"), paste0(out_dir, "XWAS_male.txt")))
names(gwamad) = NULL

# Write out list file
gwamaf = paste0(out_dir, "gwama_list.txt")
fwrite(gwamad, gwamaf, col.names = F, row.names = F, quote = F, sep = '\t')


################################################################################################################################################################
# Execute this in interactive session command line - Use Danielle's dokcer with GWAMA pre-installed: dmr07083/general-utility:1.0 
/usr/bin/GWAMA \
--map "/path/to/output/gwama.map.txt" \
--name_marker MARKER \
--indel_alleles \
--random \
--filelist "/path/to/output/gwama_list.txt" \
--output "/path/to/output/name_of_GWAMA_output"



