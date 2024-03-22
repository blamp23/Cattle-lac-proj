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


myColors <- c(
  "v" = "#F8766D",
  "mp" = "#00BA38",
  "lp" = "#619CFF",
  "el" = "#F564E3", 
  "l" = "#FFC000" 
)


# Pivot the data to a long format
ld_goi <- pivot_longer(ncdf, cols = -gene, names_to = "group", values_to = "value")

# Clean group names by removing numeric suffixes
ld_goi$group <- gsub("(_\\d+)", "", ld_goi$group)

# Specify gene of interest######################################################################
goi <- 'MRPL16'
goi_model_vector <- ""

# Iterate through mapped_motif_index to find the model vector for goi
for(model in names(mapped_motif_index)) {
  if(goi %in% mapped_motif_index[[model]]) {
    goi_model_vector <- model
    break
  }
}

# Filter for gene of interest
ld_gene_of_interest <- ld_goi %>% 
  filter(gene == goi) %>%
  filter(complete.cases(.))

# Order the 'group' factor 
ld_gene_of_interest$group <- factor(ld_gene_of_interest$group, levels = c("v", "mp", "lp", "el", "l"))

#Log
ld_gene_of_interest$value_log10 <- log10(ld_gene_of_interest$value)
ggplot(ld_gene_of_interest, aes(x = group, y = value_log10, fill = group)) +
  geom_boxplot() +
  theme_minimal() +
  geom_point(color="black", size=2, position = position_dodge(width = 0.75)) +
  geom_smooth(aes(group = 1, fill = NULL), method = "loess", color = "navy", se = FALSE, show.legend = FALSE, size = 2.25) +
  stat_summary(fun = mean, geom = "point", color = "red", size = 3, shape = 18, show.legend = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 20),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15)) +
  scale_x_discrete(labels = c("V", "MP", "LP", "EL", "L")) +
  labs(title = paste("Expression of", goi,"|",goi_model_vector),
       x = "Time Point",
       y = bquote(~Log[10]~ 'Normalized Counts')) +
  scale_fill_manual(values = myColors)


#nolog
ggplot(ld_gene_of_interest, aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  theme_minimal() +
  geom_point(color="black", size=2, position = position_dodge(width = 0.75)) +
  geom_smooth(aes(group = 1, fill = NULL), method = "loess", color = "navy", se = FALSE, show.legend = FALSE, size = 2.25) +
  stat_summary(fun = mean, geom = "point", color = "red", size = 3, shape = 18, show.legend = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 20),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15)) +
  scale_x_discrete(labels = c("V", "MP", "LP", "EL", "L")) +
  labs(title = paste("Expression of", goi,"|",goi_model_vector),
       x = "Time Point",
       y = bquote('Normalized Counts')) +
  scale_fill_manual(values = myColors)

#ISDS ISDD


















contrast_results$`v to mp`[rownames(contrast_results$`v to mp`) == "ENSBTAG00000013239", ]
averages_df[rownames(averages_df) == "ENSBTAG00000013239", ]
normalized_counts_df[rownames(normalized_counts_df) == "ENSBTAG00000013239", ]

df[rownames(df) == "ENSBTAG00000013239", ]

