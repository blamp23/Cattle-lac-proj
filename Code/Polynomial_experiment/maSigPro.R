# Load the maSigPro library
library(maSigPro)

# Retrieve column names from the normalized count matrix
column_names <- colnames(normalized_counts)

# Extract prefix and suffix from column names using regular expression
parts <- strcapture("([a-z]+)(_[0-9]+)", column_names, 
                    proto = list(prefix = character(), suffix = character()))

# Define the order of experimental timepoints
prefix_order <- c("v", "mp", "lp", "el", "l")

# Compute indices for reordering based on timepoints and replicate numbers
ordered_indices <- order(match(parts$prefix, prefix_order), 
                         as.numeric(sub("_", "", parts$suffix)))

# Reorder columns in the normalized count matrix according to the computed indices
normalized_counts_ordered <- as.matrix(normalized_counts[, ordered_indices])

# Print the reordered data frame to check the output
print(normalized_counts_ordered)

# Construct the experimental design data frame
design <- data.frame(
  Sample = colnames(normalized_counts_ordered),
  Time = factor(rep(c("V", "MP", "LP", "EL", "L"), each = 7))
)

# Map timepoints to numerical values for analysis
time_mapping <- c(V = 1, MP = 2, LP = 3, EL = 4, L = 5)
design$Time <- as.factor(design$Time)
design$nTime <- as.numeric(time_mapping[as.character(design$Time)])

# Print the design data frame to verify its structure
print(design)

# Generate the design matrix using 'nTime' column
edesign <- make.design.matrix(design, time.col = "nTime", degree = 4)
colnames(edesign$dis) <- c("nTime", "nTime2")
# Assigning correct row names to edesign$dis based on the sample identifiers
rownames(edesign$dis) <- colnames(normalized_counts_ordered)


# Fit the model using the maSigPro function 'p.vector'
fit <- p.vector(normalized_counts_ordered,
                edesign, Q = 0.05,
                MT.adjust = "BH",
                min.obs = 20)

tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)

print(colnames(edesign$dis))
colnames(edesign$dis) <- c("nTime", "nTime2", "nTime3", "nTime4")
fit <- p.vector(normalized_counts_ordered,
                edesign, Q = 0.05,
                MT.adjust = "BH",
                min.obs = 20)
print(fit$design)
str(tstep)

get<- get.siggenes(tstep, vars = 'all')
see.genes(gene_matrix, design_data, time.var = "Time")

