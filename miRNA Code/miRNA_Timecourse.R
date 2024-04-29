

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
  c("condition", "PL", "V"),
  #
  c("condition", "EL", "V"),
  c("condition", "PL", "MP"),
  #
  c("condition", "LP", "V"),
  c("condition", "EL", "MP"),
  c("condition", "PL", "LP"),
  #
  c("condition", "MP", "V"),
  c("condition", "LP", "MP"),
  c("condition", "EL", "LP"),
  c("condition", "PL", "EL")
)

contrast_results <- list()

for (i in seq_along(contrasts)) {
  contrast = contrasts[[i]]
  res <- results(dds,
                 independentFiltering = FALSE,
                 contrast = contrast)
  
  # Get the means for the conditions, ensure they are single values
  condition1_mean <- mean(averages_df[[contrast[3]]])  # Calculate mean if not already a single value
  condition2_mean <- mean(averages_df[[contrast[2]]])  # Calculate mean if not already a single value
  
  # Check if both means are less than 100
  if (condition1_mean < 100 & condition2_mean < 100) {
    expressionDirection <- "S"
  } else {
    # Determine if the mean expression increases or decreases
    expressionDirection <- ifelse(condition1_mean < condition2_mean, "I", "D")
  }
  
  # Hypothesis test: if FTR null assign 'S', else assign I/D based on expressionDirection
  res$expressionChange <- ifelse(res$padj > 0.05, "S", expressionDirection)
  
  # Store results in a list
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
# for Quaternary tree structure ###############################################################
combined_results$primary <- combined_results[,2]
combined_results$secondary <- paste0(combined_results[,3], combined_results[,4])
combined_results$tertiary <- paste0(combined_results[,5], combined_results[,6], combined_results[,7])
combined_results$quaternary <- paste0(combined_results[,8], combined_results[,9], combined_results[,10], combined_results[,11])

head(combined_results)

# Genes indexed by model vector ################################################
motif_index <- list()
x <- combined_results
vector <- paste0(x$`V to MP`,
                 x$`V to LP`,
                 x$`MP to EL`,
                 x$`LP to EL`,
                 x$`MP to PL`,
                 x$`LP to PL`,
                 x$`EL to PL`)

# Loop through each unique modelVector
for(model in unique(vector)) {
  # Subset the genes that match the current modelVector
  genes <- combined_results$gene[vector == model]
  # Store
  motif_index[[model]] <- genes
}
summary(motif_index)
mi_miRNA <- motif_index

# Search commands ##############################################################
summary(motif_index)
contrast_results$`LP to EL`[rownames(contrast_results$`LP to EL`) == "bta.miR.2285d", ]
motif_index$ISNAS

# Add vector Column to DF

# Convert row names to a column for joining
normalized_counts_df$gene <- rownames(normalized_counts_df)

# Assuming combined_results also has gene names as the first column
# Create a new data frame with only gene and modelVector columns
model_vector_df <- data.frame(gene = combined_results$gene, modelVector = combined_results$modelVector)

# Merge the data frames
normalized_counts_df1 <- merge(normalized_counts_df, model_vector_df, by = "gene", all.x = TRUE)

# Now normalized_counts_df should have the modelVector column added
############################
# Create a new workbook
wb <- createWorkbook()

# Loop through each item in the list of lists and create a sheet
for (name in names(motif_index)) {
  # Create a new sheet in the workbook with the name of the list key
  addWorksheet(wb, name)
  
  # Write the vector of genes into the sheet
  writeData(wb, name, x = motif_index[[name]], startCol = 1, startRow = 1)
}

# Define the file path and name of your Excel file
file_path <- "C:/Users/12142/Downloads/miRNA_motif_index.xlsx"

# Save the workbook
saveWorkbook(wb, file_path, overwrite = TRUE)




