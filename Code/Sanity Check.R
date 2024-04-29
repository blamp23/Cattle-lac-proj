# Sanity check

totalgene <- row.names(combined_results)
test_stable <- motif_index$SSSSSSS
true_stable <- motif_index$SSSSSSSSSS
unique_in_true <- setdiff(test_stable, true_stable)
mapped_unique <- mapIds(org.Bt.eg.db,
                        keys = unique_in_true,
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multivals = 'first')
write.table(mapped_unique, file = "mapped_motif_index_SII.txt", sep = "\t", quote = FALSE, row.names = FALSE)
unique_in_true


nons_motif_index <- lapply(mapped_motif_index, function(genes) {
  # Filter genes where the names (identifiers) are in unique_in_true
  genes[names(genes) %in% unique_in_true]
})

# Optionally, remove empty elements if there are any
nons_motif_index <- nons_motif_index[sapply(nons_motif_index, length) > 0]

summary(nons_motif_index)
