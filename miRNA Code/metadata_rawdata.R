library(readr)
library(DESeq2)

# Read the CSV file
Filtered_Raw_Counts_fix <- read_csv("Filtered_Raw_Counts_fix.csv")
Filtered_Raw_Counts_fix <- as.data.frame(Filtered_Raw_Counts_fix)
row.names(Filtered_Raw_Counts_fix) <- Filtered_Raw_Counts_fix$...1
Filtered_Raw_Counts_fix <- Filtered_Raw_Counts_fix[, -1]
df <- as.data.frame(t(Filtered_Raw_Counts_fix))


samples <- data.frame(
  row.names = c(colnames(df)),
  condition = c((substr(colnames(df), 5, nchar(colnames(df)))))
)



df <- as.matrix(df)
dds <- DESeqDataSetFromMatrix(
  countData = df,
  colData = samples,
  design = ~ condition
)

dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized = TRUE)
