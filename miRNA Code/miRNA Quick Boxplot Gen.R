# Function to check if a package is installed, and install it if it's not
check_and_install <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

# List of required packages
required_packages <- c("dplyr", "tidyr", "ggplot2", "plotly", "ggrepel")

# Check and install each required package
for (package in required_packages) {
  check_and_install(package)
}

# Required Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(ggrepel)



#####################################################################################
# miRNA must be input using the below structure
gene_input <- "bta.let.7a.5p"
#####################################################################################



# Download and load data
url <- "https://raw.githubusercontent.com/blamp23/Cattle-lac-proj/main/Raw_data/normalized_mir_counts_df.rds"
destfile <- "normalized_mir_counts_df.rds"
download.file(url, destfile, method = "auto")
ncdf1 <- readRDS(destfile)

url2 <- "https://raw.githubusercontent.com/blamp23/Cattle-lac-proj/main/Raw_data/mi_miRNA_git.rds"
destfile2 <- "mi_miRNA_git.rds"
download.file(url2, destfile_mir, method = "auto")
mi_mRNA_git <- readRDS(destfile_mir)

ncdf <- ncdf1
mapped_motif_index <- mi_mRNA_git

col_names <- c(rep("V", 11),
               rep("MP", 8),
               rep("LP", 8),
               rep("EL", 8),
               rep("PL", 13))

colnames(ncdf) <- col_names
head(ncdf)
head(mi_mRNA_git)



ncdf$gene <- rownames(ncdf)
ld_goi <- pivot_longer(ncdf, cols = -gene, names_to = "group", values_to = "value")
ld_goi$group <- gsub("(_\\d+)", "", ld_goi$group)

myColors <- c(
  "V" = "#F8766D",
  "MP" = "#00BA38",
  "LP" = "#619CFF",
  "EL" = "#F564E3", 
  "PL" = "#FFC000" 
)
# Plot the data
goi <- gene_input
goi_model_vector <- ""

if (any(sapply(mi_mRNA_git, function(x) goi %in% x))) {
  for (model in names(mi_mRNA_git)) {
    if (goi %in% mi_mRNA_git[[model]]) {
      goi_model_vector <- model
      break
    }
  }
} else {
  print("Gene not found. Please check the gene name and try again.")
}

ld_gene_of_interest <- ld_goi %>% 
  filter(gene == goi) %>%
  filter(complete.cases(.)) %>%
  mutate(group = factor(group, levels = c("V", "MP", "LP", "EL", "PL")))

ggplot(ld_gene_of_interest, aes(x = group, y = log10(value), fill = group)) +
  geom_boxplot() +
  theme_minimal() +
  geom_point(color="black", size=2, position = position_dodge(width = 0.75)) +
  geom_smooth(aes(group = 1), method = "loess", color = "navy", se = FALSE, size = 2.25) +
  stat_summary(fun = mean, geom = "point", color = "red", size = 3, shape = 18) +
  theme(plot.title = element_text(size = 18, face = 'bold'), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Adjust x-axis text size
        axis.text.y = element_text(size = 12),  # Adjust y-axis text size
        axis.title = element_text(size = 14),   # Adjust axis title size
        legend.title = element_text(size = 14), # Adjust legend title size
        legend.text = element_text(size = 12),  # Adjust legend text size
        legend.key.size = unit(3, "lines")) +   # Adjust legend key size
  scale_x_discrete(labels = c("V", "MP", "LP", "EL", "L")) +
  labs(title = paste("Expression of", goi, "|", goi_model_vector),
       x = "Time Point",
       y = bquote(~Log[10]~ 'Normalized Counts'),
       fill = "Time Point") +   # Rename the legend
  scale_fill_manual(values = myColors,
                    labels = c("Virgin", "Mid-Pregnant", "Late-Pregnant", "Early Lactation", "Peak Lactation"))
