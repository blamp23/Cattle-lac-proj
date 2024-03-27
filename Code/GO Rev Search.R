library(jsonlite)
# Import From JSON #############################################################
# MYC-19030024-MESC-mouse
url <- "https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/MYC-19030024-MESC-mouse/CHEA+Transcription+Factor+Binding+Site+Profiles"
data <- fromJSON(url)
str(data) # Structure may change based on import site, adjust GOI path as needed
GOI <- data$associations$gene$symbol
GOI <- mapIds(org.Bt.eg.db,
              keys = GOI,
              column = "ENSEMBL",
              keytype = "SYMBOL",
              multivals = 'first')

# Import From TSV ##############################################################
# Taxon ID. 9913, If you want to filter based on BosTau
# rRNA metabolic process
file_path <- "C:/Users/12142/Downloads/QuickGO-annotations-1711470871009-20240326.tsv"
chrom_genes <- read.delim(file_path, sep = "\t")
head(chrom_genes)
chrom_genes$ensembl <- mapIds(org.Bt.eg.db,
                       keys = chrom_genes$SYMBOL, # Adjust keys for Symbol
                       column = "ENSEMBL",
                       keytype = "SYMBOL",
                       multivals = 'first')
GOI <- c(chrom_genes$ensembl)

##########################################################
# Myc Targets V1 hallmark MSigDB
url <- "C:/Users/12142/Downloads/HALLMARK_MYC_TARGETS_V1.v2023.2.Hs.json"
data <- fromJSON(url)
str(data) # Structure may change based on import site, adjust GOI path as needed
GOI <- data$HALLMARK_MYC_TARGETS_V1$geneSymbols
GOI <- mapIds(org.Bt.eg.db,
              keys = GOI,
              column = "ENSEMBL",
              keytype = "SYMBOL",
              multivals = 'first')


# Filter for GOI ###############################################################
filtered_motif_index <- lapply(motif_index, function(genes) {
  genes[genes %in% GOI]
})

smi <- summary(motif_index)
sfmi <- summary(filtered_motif_index)

smi
sfmi


# Proportionality Table ########################################################
# Convert to data frames and convert 'Length' to numeric
# Length rep # of genes
smi_df <- data.frame(Motif = rownames(smi), Length = as.numeric(smi[, "Length"]))
sfmi_df <- data.frame(Motif = rownames(sfmi), Length = as.numeric(sfmi[, "Length"]))

# Ensure the motifs are in the same order for both data frames
smi_df <- smi_df[order(smi_df$Motif),]
sfmi_df <- sfmi_df[order(sfmi_df$Motif),]

# Calculate proportions
proportions <-sfmi_df$Length / smi_df$Length

# Combine results
final_df <- data.frame(Motif = smi_df$Motif, Proportion = proportions)
final_df$Proportion <- final_df$Proportion * 100 # Proportion is a % of 100
final_df <- final_df[order(-final_df$Proportion),]
print(final_df)
smi
filtered_motif_index[['ISDS']]







