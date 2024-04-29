ndf <- as.data.frame(normalized_counts_ordered)
meta <- data.frame(sample = samples$sample,
                   condition = samples$condition,
                   time = samples$time)

library(dplyr)
library(tidyr)
library(ggplot2)

# Convert row names to a column
ndf$gene_id <- rownames(ndf)

# Reshape the data to a long format
data_long <- pivot_longer(ndf, cols = -gene_id, names_to = "sample", values_to = "expression")

# Merge with metadata to align with timepoints
data_long <- left_join(data_long, meta, by = c("sample" = "sample"))


# After merging, check the data
print(head(data_long))
print(summary(data_long))
# Calculate the variance or standard deviation for each gene across timepoints
variance_data <- data_long %>%
  group_by(gene_id) %>%
  summarise(variance = var(expression))

# Filter out genes with variance below a certain threshold
threshold = quantile(variance_data$variance, 0.25) # for example, lower quartile
dynamic_genes <- variance_data %>%
  filter(variance > threshold) %>%
  pull(gene_id)

# Use only dynamic genes for further analysis
data_dynamic <- data_long %>%
  filter(gene_id %in% dynamic_genes)

# Fit the polynomial regression and handle potential errors or insufficient data warnings
poly_fits <- data_dynamic %>%
  group_by(gene_id) %>%
  do({
    # Ensure enough points are available and data is clean
    if(nrow(.) >= 5 && sum(is.na(.$expression)) == 0 && sum(is.na(.$time)) == 0) {
      tidy(lm(expression ~ poly(time, 4), data = .))
    } else {
      data.frame(term=character(), estimate=numeric(), std.error=numeric(), statistic=numeric(), p.value=numeric())
    }
  })


# Extracting coefficients for clustering, handling missing data
coefficients <- poly_fits %>%
  filter(grepl("poly", term)) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  group_by(gene_id) %>%
  summarise(across(starts_with("poly"), mean, na.rm = TRUE)) # summarizing to ensure one row per gene

# Clustering based on the polynomial coefficients
set.seed(123) # for reproducibility
clusters <- kmeans(dplyr::select(coefficients, starts_with("poly")), centers = 5) # adjust the number of centers as needed

# Add cluster labels to the coefficients
coefficients$cluster <- clusters$cluster

# Review the results
print(head(coefficients))
View(coefficients
     )
summary(coefficients)


