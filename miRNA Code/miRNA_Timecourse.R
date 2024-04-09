# Libraries ##########################################################################
library(DESeq2)
library(dplyr)

# Normalized Count table and dds object from DESeq are needed to run code
col_names <- c(rep("V", 11),
               rep("MP", 8),
               rep("LP", 8),
               rep("EL", 8),
               rep("PL", 13))
names(normalized_counts) <- col_names


# Time-point Averages of Normalized Counts ######################################
normalized_counts_df <- as.data.frame(normalized_counts)


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


# pairwise comparisons, hypothesis testing, and assignment of expression pattern ###################################
contrasts <- list(
  c("condition", "MP", "V"),
  c("condition", "LP", "MP"),
  c("condition", "EL", "LP"),
  c("condition", "PL", "EL")
)

contrast_results <- list()

# DE for pairs in contrasts
for (i in seq_along(contrasts)) {
  contrast = contrasts[[i]]
  res <- results(dds,
                 independentFiltering = FALSE,
                 contrast = contrast)
  
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

# Search commands ##############################################################
summary(motif_index)
contrast_results$`LP to EL`[rownames(contrast_results$`LP to EL`) == "bta.miR.2285d", ]
motif_index$ISNAS

