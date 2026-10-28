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
# Execute this in interactive session command line - Use Docker with GWAMA pre-installed: dmr07083/general-utility:1.0 

#/usr/bin/GWAMA \
#--map "/path/to/output/gwama.map.txt" \
#--name_marker MARKER \
#--indel_alleles \
#--random \
#--filelist "/path/to/output/gwama_list.txt" \
#--output "/path/to/output/name_of_GWAMA_output"


################################################################################
# Update GWAMA output to correct format for downstream analyses
meta = fread("/path/to/output/name_of_GWAMA_output.out")

colnames(meta) <- c("CHR","BP","posID","ALLELE1","ALLELE0","A1FREQ","BETA","SE",
                    "BETA_95L","BETA_95U","Z","P","P_neg_log10","Q_V","Q_P","i2",
                    "N_studies","N_incl","META_DIR")

# Add back rsIDs
meta2 = merge(meta, key, by.x = "posID", by.y = "MARKER") 

# filter to 90% genotype cut-off
meta90 = meta2 %>% 
  mutate(DIR1 = substr(META_DIR, 1, 1)) %>% 
  filter(DIR1 != "?",
        N_incl > max(N_incl)*0.9) %>% 
  dplyr::select(-DIR1)

# Write out file
fwrite(meta90, paste0(out_dir, 'stats.gen090'))



