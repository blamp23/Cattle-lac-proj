# motif model 
library(ggplot2)

# Example motifs input
input_motifs <- c("SISD", "SIID", "ISID", "SSSD", "DIIS")


# Initialize an empty data frame to store motif trends
motif_data <- data.frame(Timepoint = integer(), Value = numeric(), Motif = character())

# Function to process each motif and convert to numeric values
process_motif <- function(motif) {
  # Starting value
  value <- 0
  values <- value
  
  # Process each character in the motif
  for (char in strsplit(motif, "")[[1]]) {
    if (char == "I") {
      value <- value + 1
    } else if (char == "D") {
      value <- value - 1
    }
    # "S" does nothing, so it's omitted
    values <- c(values, value)
  }
  
  return(values)
}

# Process each motif and add to the data frame
for (motif in input_motifs) {
  values <- process_motif(motif)
  timepoints <- 0:(length(values) - 1) # Adjust timepoint labels as needed
  motif_data <- rbind(motif_data, data.frame(Timepoint = as.factor(timepoints), Value = values, Motif = motif))
}

# Plot the motifs
ggplot(motif_data, aes(x = Timepoint, y = Value, group = Motif, color = Motif)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Motif Trends",
       x = "Timepoint",
       y = "Value",
       color = "Motif") +
  scale_color_brewer(palette = "Set1")
