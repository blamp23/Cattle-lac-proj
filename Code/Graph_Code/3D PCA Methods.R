# 3D PCA Plots / Loadings
# C:/Users/12142/Bioinformatics/Spring_lactation_proj/3D PCA Methods.R
# Edited 2/16/2024
# Libraries ####################################################################
library(DESeq2)
library(plotly)
library(dplyr)
library(org.Bt.eg.db)
library(AnnotationDbi)
# 3D PCA #######################################################################
vsd <- vst(dds) # variance stabilizing transformation of Deseq object
pcaResult <- prcomp(t(assay(vsd))) # transpose vsd and perform PCA
pcaData3D <- as.data.frame(pcaResult$x[, 1:3]) # Extract the first three principal components
pcaData3D$condition <- vsd$condition # New condition column is added

# Calculate the percent variance for the first three components
# squares the standard deviations of the PCA results, divides by the total variance, and multiplies by 100
percentVar <- round(100 * pcaResult$sdev^2 / sum(pcaResult$sdev^2)) 

# Create a 3D plot with Plotly
fig <- plot_ly(data = pcaData3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~condition, 
               type = 'scatter3d', mode = 'markers', text = rownames(pcaData3D)) %>%
  layout(title = "3D PCA Plot",
         scene = list(
           xaxis = list(title = paste0("PC1: ", percentVar[1], "% variance")),
           yaxis = list(title = paste0("PC2: ", percentVar[2], "% variance")),
           zaxis = list(title = paste0("PC3: ", percentVar[3], "% variance"))
         ))
fig



# 3d PCA DESeq Method ##########################################################
vsd <- vst(dds)
top_var_genes <- order(rowVars(assay(vsd)), decreasing = TRUE)[1:500] # Calculate most var genes and add them to a new obj
vsd_top_var <- vsd[top_var_genes, ] # pull vsd from only top var genes

pcaResult <- prcomp(t(assay(vsd_top_var)))
pcaData3D <- as.data.frame(pcaResult$x[, 1:3])
pcaData3D$condition <- vsd_top_var$condition  
percentVar <- round(100 * pcaResult$sdev^2 / sum(pcaResult$sdev^2))

fig <- plot_ly(data = pcaData3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~condition, 
               type = 'scatter3d', mode = 'markers', text = rownames(pcaData3D)) %>%
  layout(title = "3D PCA Plot",
         scene = list(
           xaxis = list(title = paste0("PC1: ", percentVar[1], "% variance")),
           yaxis = list(title = paste0("PC2: ", percentVar[2], "% variance")),
           zaxis = list(title = paste0("PC3: ", percentVar[3], "% variance"))
         ))
fig

# lil 2D PCA moment ############################################################
vsd <- vst(dds)
# Get the PCA data
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)

# Extract the percent variance captured by each principal component
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Plot with ggplot2
ggplot(pcaData, aes(x=PC1, y=PC2, color=condition)) +
  geom_point(size=3) +
  geom_text_repel(aes(label=row.names(pcaData))) +  # using ggrepel to prevent overlap
  stat_ellipse(aes(group=condition), level=0.85) +  # adding ellipses around each condition group
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA Plot")

# Loadings #####################################################################
# Extracting loadings from PCA result
loadings <- pcaResult$rotation
loadings_PC123 <- loadings[, 1:3]
# Convert to data frame
loadings_PC123 <- as.data.frame(loadings_PC123)

# Using mapIds to get gene symbols
symbols <- mapIds(org.Bt.eg.db,
                  keys=row.names(loadings_PC123),
                  column="SYMBOL",
                  keytype="ENSEMBL",
                  multivals='first')
symbols_vector <- as.character(symbols)
loadings_PC123$Symbol <- symbols_vector

# Double check and save
head(loadings_PC123)
write.csv(loadings_PC123, "PCA_Loadings_PC123.csv")
