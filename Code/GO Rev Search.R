
# MYC MAX complex
file_path <- "C:/Users/12142/Downloads/QuickGO-annotations-1711462600951-20240326.tsv"
chrom_genes <- read.delim(file_path, sep = "\t")
head(chrom_genes)

chrom_genes$ensembl <- mapIds(org.Bt.eg.db,
                       keys = chrom_genes$SYMBOL,
                       column = "ENSEMBL",
                       keytype = "SYMBOL",
                       multivals = 'first')
gene_of_interest <- c(chrom_genes$ensembl)
filtered_motif_index <- lapply(motif_index, function(genes) {
  genes[genes %in% genes_of_interest]
})

smi <- summary(motif_index)
sfmi <- summary(filtered_motif_index)

# Convert to data frames and convert 'Length' to numeric
smi_df <- data.frame(Motif = rownames(smi), Length = as.numeric(smi[, "Length"]))
sfmi_df <- data.frame(Motif = rownames(sfmi), Length = as.numeric(sfmi[, "Length"]))

# Ensure the motifs are in the same order for both data frames
smi_df <- smi_df[order(smi_df$Motif),]
sfmi_df <- sfmi_df[order(sfmi_df$Motif),]

#Calculate proportions
proportions <-sfmi_df$Length / smi_df$Length

# Combine results
final_df <- data.frame(Motif = smi_df$Motif, Proportion = proportions)
final_df$Proportion <- final_df$Proportion * 100
final_df <- final_df[order(-final_df$Proportion),]
print(final_df)

filtered_motif_index[['SSIS']]






