library(readr)
library(DESeq2)

# Read the CSV file
file_path <- 'https://raw.githubusercontent.com/blamp23/Cattle-lac-proj/main/Raw_data/Filtered_Raw_Counts_fix.csv'
df <- read.csv(file_path)
row.names(df) <- df$X
df <- df[-1]
df <- t(df)


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

