# Ensure all necessary libraries are loaded
library(AnnotationDbi)
library(org.Bt.eg.db)
library(dplyr)
library(tidyr)
library(ggplot2)

ncdf <- normalized_counts
# Save original Ensembl IDs
original_ensembl_ids <- row.names(ncdf)

# Map
gene_symbols <- mapIds(org.Bt.eg.db,
                       keys = original_ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multivals = 'first')

# Replace NAs with the original Ensembl ID
gene_symbols[is.na(gene_symbols)] <- original_ensembl_ids[is.na(gene_symbols)]

# Update row names in ncdf
row.names(ncdf) <- gene_symbols
ncdf <- as.data.frame(ncdf)
ncdf$gene <- rownames(ncdf)


# Pivot the data to a long format
ld_goi <- pivot_longer(ncdf, cols = -gene, names_to = "group", values_to = "value")

# Clean group names by removing numeric suffixes
ld_goi$group <- gsub("(_\\d+)", "", ld_goi$group)

# Specify gene of interest
goi <- 'CSN1S1'

# Filter for gene of interest
ld_gene_of_interest <- ld_goi %>% 
  filter(gene == goi) %>%
  filter(complete.cases(.))

# Order the 'group' factor 
ld_gene_of_interest$group <- factor(ld_gene_of_interest$group, levels = c("v", "mp", "lp", "el", "l"))

# Plotting
ggplot(ld_gene_of_interest, aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  theme_minimal() +
  geom_jitter(color="black", size=2, alpha=0.9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = c("V", "MP", "LP", "EL", "L")) +
  labs(title = paste("Expression of", unique(ld_gene_of_interest$gene)),
       x = "Group",
       y = "Expression")

contrast_results$`v to mp`[rownames(contrast_results$`v to mp`) == "ENSBTAG00000019262", ]
averages_df[rownames(averages_df) == "ENSBTAG00000019262", ]
normalized_counts_df[rownames(normalized_counts_df) == "ENSBTAG00000004943", ]
df[rownames(df) == "ENSBTAG00000004943", ]

