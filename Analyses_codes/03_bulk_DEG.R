## Differential expression analysis for top XWAS genes

## Libraries
library(data.table)
library(tidyverse)
library(patchwork)
library(scales)


#####################################################
# Read in gene reference
gene_ref = fread("/path/to/gene/reference/file/gene_reference.csv") %>% 
  dplyr::select(gene_name, ens_id) 


#####################################################
### Read in DEG results - downloaded from Synapse
deg_results = fread("/path/to/DEG/results/DEG_results.csv") %>% 
  filter(chromosome_name == "X") %>% 
  dplyr::rename(gene_name = hgnc_symbol,
                ens_id = ensembl_gene_id) %>% 
  # Filter to comaprison of interest (AD vs Control in tissue)
  filter(Comparison == "diagnosis_tissue*age_death.AD_CBE - diagnosis_tissue*age_death.CT_CBE") %>% 
  mutate(abs_t = abs(t))


# Merge the prioritized genes with deg_results
prioritized_genes_df = gene_ref %>%
  inner_join(deg_results, by = "ens_id") 



######################################################################
## Define variables for enrichment
sample_n = n_distinct(prioritized_genes_df$gene_name)
sample_mean = mean(prioritized_genes_df$abs_t)
sample_sd = sd(prioritized_genes_df$abs_t)
pop_sd = sd(deg_results$abs_t)
standard_e = pop_sd/sqrt(sample_n)

######################################################################
## Set up simulation to calculate population mean

# Get number of intersectiong genes in prioritized to reference for random sample size
sample_n = n_distinct(prioritized_genes_df$gene_name)
num_simulations = 1000
simulated_means_abs_t = numeric(num_simulations)

# Run simulations
set.seed(123)  # For reproducibility
for (i in 1:num_simulations) {
  # Randomly select 59 genes from deg_results
  random_genes <- deg_results[sample(nrow(deg_results), sample_n), ]
  # Compute the average absolute t-statistic of the randomly selected genes
  simulated_means_abs_t[i] <- mean(abs(random_genes$t))
}


### Calculate Z-score and p-value from z-test
z = (sample_mean-mean(simulated_means_abs_t))/standard_e
p = 2*pnorm(-abs(z))

print(z)
print(p)



######################################################################
### Create plots

## Create null distribution histogram
hist <- ggplot(simulated_means_abs_t_df, aes(x = mean_abs_t)) +
  geom_histogram(binwidth = 0.05,
                 fill = "grey65") +
  coord_flip() +
  scale_x_continuous(limits = c(0, 8)) +
  geom_vline(xintercept = pm, color = "red", linetype = "solid", linewidth = 3) +
  theme_classic() +
  labs(x = "Absolute Effect Size", y = "Null") +
  theme(
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12))



## Make violin plot
violin = ggplot(prioritized_genes_df, aes(x = "", y = abs_t)) +
  geom_violin(fill = "grey70", color = "black", alpha = 0.5, trim = FALSE) +
  geom_jitter(width = 0.08, alpha = 0.5, size = 2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.8, linetype = "solid", color = "red", linewidth = 1.3) +
  scale_y_continuous(limits = c(0,8)) +  
  theme_classic() +
  labs(x = "Prioritized",
       y = "|t-stat|") +
  theme(axis.title = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 14, face = "bold"),
        axis.text.y  = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(), 
        strip.text.y = element_blank())


### Use patch work to put the two plots together
comb = hist + violin + plot_layout(ncol = 2, widths = c(1, 0.6))
print(comb)
