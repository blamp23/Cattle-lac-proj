# Libraries #####################################################################
library(DESeq2)
library(apeglm)
library(ggplot2)
library(sva)
library(ggrepel)
library(VennDiagram)
library(tidyverse)
library(org.Bt.eg.db)
library(dplyr)
library(reshape2)
library(clusterProfiler)
library(enrichplot)
library(grid)
library(readxl)
library(tidyr)
library(plotly)
library(plot3D)
library(KEGGREST)
library(AnnotationDbi)
# Raw Count Data From Galaxy #####################################################
file_path <- 'C:/Users/12142/Downloads/Galaxy37-[Multi-Join_on_data_35,_data_34,_and_others].tabular'
df <- read.delim(file_path, header = TRUE, row.names = 1)
df <- df[-1]

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

# Time-point Averages of Normalized Counts ######################################
# Define/assign column names
col_names <- c(rep("el", 7), rep("l", 7), rep("lp", 7), rep("mp", 7), rep("v", 7))
names(normalized_counts_df) <- col_names

# Prepare an empty data frame for storing averages, with correct dimensions and column names
averages_df <- data.frame(matrix(ncol = length(unique(col_names)), nrow = nrow(normalized_counts_df)))
names(averages_df) <- unique(col_names)
rownames(averages_df) <- rownames(normalized_counts_df)

# Loop through each time-point to calculate and store row means
for (timepoint in unique(col_names)) {
  # Identify columns belonging to the current time-point
  cols <- which(col_names == timepoint)
  # Calculate and store
  averages_df[[timepoint]] <- rowMeans(normalized_counts_df[, cols], na.rm = TRUE)
}
head(averages_df)


# pairwise comparisons, hypothesis testing, and assignment of expression pattern ####
contrasts <- list(
  c("condition", "mp", "v"),
  c("condition", "lp", "mp"),
  c("condition", "el", "lp"),
  c("condition", "l", "el")
)

contrast_results <- list()

# DE for pairs in contrasts
for (i in seq_along(contrasts)) {
  contrast = contrasts[[i]]
  res <- results(dds, contrast=contrast)
  
  # Determine if the mean expression increases or decreases
  condition1_mean <- averages_df[[contrast[3]]]  # Mean for the first condition
  condition2_mean <- averages_df[[contrast[2]]]  # Mean for the second condition
  expressionDirection <- ifelse(condition1_mean < condition2_mean, "I", "D")
  
  # Hypothesis test, FTR null assign 'S', Rej null assign I/D
  res$expressionChange <- ifelse(res$padj > 0.05, "S", expressionDirection)
  
  contrast_results[[paste(contrast[3], "to", contrast[2])]] <- res
}
# Generation of expression motif ###############################################

# Initialize a dataframe with gene names
combined_results <- data.frame(gene=rownames(contrast_results[[1]]))

# Loop through each contrast result to combine them
for (i in seq_along(contrast_results)) {
  contrast_name <- names(contrast_results)[i]  # Get the name of the current contrast
  # Create a temporary dataframe with genes and their expression changes for the current contrast
  temp_df <- data.frame(gene=rownames(contrast_results[[i]]), 
                        expressionChange=contrast_results[[i]]$expressionChange)
  colnames(temp_df)[2] <- contrast_name  # Rename the second column to the current contrast name
  # Merge the temporary dataframe with the combined_results dataframe
  combined_results <- merge(combined_results, temp_df, by="gene", all=TRUE)
}

# Concatenate expression patterns to make model vector
combined_results$modelVector <- apply(combined_results[, -1], 1, 
                                      function(x) paste(x, collapse = ""))
head(combined_results)

# Genes indexed by model vector ################################################
motif_index <- list()
# Loop through each unique modelVector
for(model in unique(combined_results$modelVector)) {
  # Subset the genes that match the current modelVector
  genes <- combined_results$gene[combined_results$modelVector == model]
  # Store
  motif_index[[model]] <- genes
}
summary(motif_index)






