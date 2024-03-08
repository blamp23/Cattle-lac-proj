# Libraries #####################################################################
library(DESeq2)
library(sva)
library(dplyr)
library(readxl)
library(tidyr)
# Raw Count Data From Galaxy #####################################################
file_path <- 'https://raw.githubusercontent.com/blamp23/Cattle-lac-proj/main/Raw_data/Raw_Cattle_Counts.tabular'
df <- read.delim(file_path, header = TRUE, row.names = 1)
df <- df[-1]

samples
samples <- data.frame(
  row.names = c("v_468", "v_502", "v_504", "v_507", "v_509", "v_598", "v_610", 
                "mp_468", "mp_502", "mp_504", "mp_507", "mp_509", "mp_598", "mp_610",
                "lp_468", "lp_502", "lp_504", "lp_509", "lp_507", "lp_598", "lp_610",
                "el_468", "el_502", "el_504", "el_507", "el_509", "el_598", "el_610",
                "l_468", "l_502", "l_504", "l_507", "l_509", "l_598", "l_610"),
  condition = c(rep("v", 7), rep("mp", 7), rep("lp", 7), rep("el", 7), rep("l", 7)),
  animal = c(468, 502, 504, 507, 509, 598, 610,
             468, 502, 504, 507, 509, 598, 610,
             468, 502, 504, 509, 507, 598, 610,
             468, 502, 504, 507, 509, 598, 610,
             468, 502, 504, 507, 509, 598, 610),
  batch = c(1, 1, 1, 1, 1, 1, 1, 
            0, 1, 0, 1, 0, 0, 0,
            0, 1, 0, 0, 1, 0, 0,
            0, 1, 0, 1, 0, 0, 0, 
            1, 1, 1, 1, 1, 1, 1 ),
  BCM = c(1, 1, 1, 1, 1, 1, 1,
          0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0,
          1, 1, 1, 1, 1, 1, 1)
)

new_colnames <- c(row.names(samples))
colnames(df) <- new_colnames

colnames(df)[28] <- "lp_610"  # lp_610 and el_610 are swapped samples
colnames(df)[21] <- "el_610"
samples <- samples[colnames(df), ]
# Filter Method A ##############################################################
# originaly 27607 genes
# filtered to 17344 genes
# Step 1: Splitting Data by Condition
conditions <- unique(samples$condition)
split_data <- lapply(conditions, function(cond) {
  df[, samples$condition == cond]
})

# Step 2: Define Filtering Function
count_high_expressions <- function(data) {
  apply(data, 1, function(gene) {
    sum(gene >= 3) >= 4  # Returns TRUE if at least 4 out of 7 samples have counts ≥ 3
  })
}
high_exp_counts <- lapply(split_data, count_high_expressions)

# Step 3: Aggregate Results for inclusion
gene_inclusion_filter <- rowSums(do.call(cbind, high_exp_counts)) >= 1 # If 1/5 timepoints meet crit, gene is included

# Step 4: Filter Genes
df <- df[gene_inclusion_filter, ]
# Deseq2 - ComBat #########################################################################
count_matrix <- as.matrix(df)
cdf <- ComBat_seq(count_matrix, batch = samples$batch, group = samples$condition)

dds <- DESeqDataSetFromMatrix(
  countData = cdf,
  colData = samples,
  design = ~ condition
)


dds <- DESeq(dds)
res <- results(dds)


dds_count <- estimateSizeFactors(dds)
normalized_counts <- counts(dds_count, normalized=TRUE)
ordered_column_indices <- order(colnames(normalized_counts))
normalized_counts <- normalized_counts[, ordered_column_indices]


