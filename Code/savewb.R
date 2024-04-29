# Create an empty list to hold the data frames.
summary(motif_index)
list_of_dfs <- list()

# Loop through each motif in 'motif_index'
for (motif in names(motif_index)) {
  # Extract the list of Ensembl gene IDs for the current motif
  ensembl_ids <- motif_index[[motif]]
  
  # Use mapIds to get gene symbols for these Ensembl IDs
  gene_symbols <- tryCatch({
    mapIds(org.Bt.eg.db,
           keys = ensembl_ids,
           column = "SYMBOL",
           keytype = "ENSEMBL",
           multivals = 'first')
  }, error = function(e) rep(NA, length(ensembl_ids)))  # Handling missing symbols
  
  # Use mapIds to get gene names for these Ensembl IDs
  gene_names <- tryCatch({
    mapIds(org.Bt.eg.db,
           keys = ensembl_ids,
           column = "GENENAME",
           keytype = "ENSEMBL",
           multivals = 'first')
  }, error = function(e) rep(NA, length(ensembl_ids)))  # Handling missing gene names
  
  # Create a data frame for the current motif
  df <- data.frame(EnsemblID = ensembl_ids, GeneSymbol = gene_symbols, GeneName = gene_names, stringsAsFactors = FALSE)
  
  # Assign the data frame to the list using the motif name
  list_of_dfs[[motif]] <- df
}

# list_of_dfs now contains a data frame for each motif with Ensembl IDs, gene symbols, and gene names
library(openxlsx)
# Create a new workbook
wb <- createWorkbook()

# Add each data frame as a separate sheet
for (motif in names(list_of_dfs)) {
  addWorksheet(wb, motif)  # Create a worksheet named after the motif
  writeData(wb, motif, list_of_dfs[[motif]])  # Write data frame to the worksheet
}

# Save the workbook to a file
saveWorkbook(wb, "Motif_Gene_Data_add_el.xlsx", overwrite = TRUE)
