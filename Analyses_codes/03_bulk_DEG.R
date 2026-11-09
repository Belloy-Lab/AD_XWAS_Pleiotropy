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
region = "ACC"   # brain region; change for each of the eight regions
deg_results = fread("/path/to/DEG/results/DEG_results.csv") %>%
  filter(chromosome_name == "X") %>%
  dplyr::rename(gene_name = hgnc_symbol,
                ens_id = ensembl_gene_id) %>%
  # Filter to comaprison of interest (AD vs Control in tissue)
  filter(Comparison == sprintf("diagnosis_tissue*age_death.AD_%s - diagnosis_tissue*age_death.CT_%s", region, region)) %>%
  filter(!is.na(t), !is.na(AveExpr)) %>%
  mutate(abs_t = abs(t))


# Merge the prioritized genes with deg_results
prioritized_genes_df = gene_ref %>%
  inner_join(dplyr::select(deg_results, -gene_name), by = "ens_id")

# Background genes (prioritized genes excluded)
background_df = deg_results %>%
  filter(!ens_id %in% gene_ref$ens_id)


######################################################################
## Define variables for enrichment
sample_n = n_distinct(prioritized_genes_df$gene_name)
sample_mean = mean(prioritized_genes_df$abs_t)
B = 100000   # permutations
K = 30       # closest background genes (by AveExpr) for the matched background

######################################################################
## Permutation 1 - non-matched background
set.seed(20261103)
null_unmatched = replicate(B, mean(sample(background_df$abs_t, sample_n)))
p_unmatched = (1 + sum(null_unmatched >= sample_mean)) / (B + 1)

######################################################################
## Permutation 2 - background matched on average expression
# K closest background genes to each prioritized gene
closest = t(sapply(prioritized_genes_df$AveExpr, function(e) order(abs(background_df$AveExpr - e))[1:K]))

set.seed(20261103)
draw_ids = function(nb) {
  idx = matrix(sample.int(K, sample_n * nb, replace = TRUE), nrow = sample_n)
  matrix(closest[cbind(rep(seq_len(sample_n), nb), as.vector(idx))], nrow = sample_n)
}
ids = draw_ids(B)
# Redraw permutations that sample the same gene twice
dup = function(m) apply(m, 2, anyDuplicated) > 0
bad = which(dup(ids))
while (length(bad) > 0) {
  ids[, bad] = draw_ids(length(bad))
  bad = bad[dup(ids[, bad, drop = FALSE])]
}
null_matched = colMeans(matrix(background_df$abs_t[ids], nrow = sample_n))
p_matched = (1 + sum(null_matched >= sample_mean)) / (B + 1)

print(p_unmatched)
print(p_matched)

# Create a dataframe for null plot
simulated_means_abs_t_df <- data.frame(mean_abs_t = null_unmatched)
pm = mean(simulated_means_abs_t_df$mean_abs_t)
sm = mean(prioritized_genes_df$abs_t)

######################################################################
### Create plots (non-matched background)

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
