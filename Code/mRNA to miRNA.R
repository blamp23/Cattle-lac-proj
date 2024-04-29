# mi_mrna = motif index from mRNA dataset
# mi_miRNA = motif index from miRNA dataset

# Assuming that an increase in expression of miRNA would suppress expression of mRNA,
# we can 'flip' the expression patterns to get a better understanding of 
# what mirna patterns would correspond to mrna

# In order for this to work, motif structure of mi_mrna and mi_mirna need to be the same. 

# Lib
library(openxlsx)
######
str(mi_mrna)
str(mi_miRNA)
fmi_miRNA <- mi_miRNA

flip_characters <- function(key) {
  # Replace all 'D' with a temporary character, here we use 'X'
  temp_key <- gsub("D", "X", key)
  # Replace all 'I' with 'D'
  temp_key <- gsub("I", "D", temp_key)
  # Replace all 'X' with 'I'
  flipped_key <- gsub("X", "I", temp_key)
  return(flipped_key)
}

# Apply the function to the keys of your miRNA list
names(fmi_miRNA) <- sapply(names(fmi_miRNA), flip_characters)


# Find common keys
common_keys <- intersect(names(mi_mrna), names(fmi_miRNA))
tempck <- intersect(names(fmi_miRNA), names(mi_mrna))

# Create a new list to store the data structure
combined_data <- list()

# Populate the data structure with keys, miRNAs, and mRNAs
for (key in common_keys) {
  combined_data[[key]] <- list(
    miRNAs = fmi_miRNA[[key]],
    mRNAs = mi_mrna[[key]]
  )
}
############################################################################
pattern_summary <- data.frame()
pattern_summary <- data.frame(
  Pattern = character(),
  miRNA = integer(),
  mRNA = integer(),
  stringsAsFactors = FALSE 
)

# Loop through the combined_data to populate the DataFrame
for (key in names(combined_data)) {
  miRNA <- length(combined_data[[key]]$miRNAs)
  mRNA <- length(combined_data[[key]]$mRNAs)
  
  # Append to the data frame
  pattern_summary <- pattern_summary %>%
    add_row(Pattern = key, miRNA = miRNA, mRNA = mRNA)
}


print(pattern_summary)
View(pattern_summary)

####################################################################
nested_list <- combined_data

# Loop through each item and update mRNAs with gene symbols
for (item in names(nested_list)) {
  # Extract Ensembl IDs
  ensembl_ids <- nested_list[[item]]$mRNAs
  
  # Map Ensembl IDs to gene symbols
  gene_symbols <- mapIds(org.Bt.eg.db,
                         keys = ensembl_ids,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multivals = 'first')
  
  # Update the mRNAs list with gene symbols
  nested_list[[item]]$mRNAs <- gene_symbols
}

str(nested_list)


# Create a new workbook
wb <- createWorkbook()

# Add data from nested_list to the workbook
for (key in names(nested_list)) {
  # Create a sheet for each key
  addWorksheet(wb, key)
  
  # Get miRNA and mRNA data
  miRNAs <- nested_list[[key]]$miRNAs
  mRNAs <- nested_list[[key]]$mRNAs
  mRNANames <- names(mRNAs) # Assuming you want to include mRNA Ensembl IDs as well
  
  # Create a data frame for miRNAs and mRNAs
  # If miRNA and mRNA lists are not the same length, adjust by filling shorter list with NA
  max_length <- max(length(miRNAs), length(mRNAs))
  miRNAs <- c(miRNAs, rep(NA, max_length - length(miRNAs)))
  mRNAs <- c(mRNAs, rep(NA, max_length - length(mRNAs)))
  mRNANames <- c(mRNANames, rep(NA, max_length - length(mRNANames)))
  
  data <- data.frame(miRNA = miRNAs, mRNA = mRNAs, mRNA_ID = mRNANames)
  
  # Write the data frame to the currently active sheet
  writeData(wb, key, data)
}

# Save the workbook
# Save the workbook, specifying the full path where it should be saved
saveWorkbook(wb, "C:\\Users\\12142\\Downloads\\tmp_mRNA to miRNA.xlsx", overwrite = TRUE)



