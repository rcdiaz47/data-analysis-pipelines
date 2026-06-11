# ===================================
# Random Forest Classification for omics data 
# Demonstrates an end-to-end ML workflow for 
# classifying samples by experimental group 

# ----- Libraries -----
library(tidyverse) # data manipulation
library(randomForest) # Random Forest algorithm 
library(caret) # train/test splitting and model evaluation 
library(pheatmap) # visualizing important features 

# =====================================
# USER CONFIGURATION - Edit this section only
# =====================================

# ----- Data Ingestion (dataset specific)-----
# Raw metabolomics workbench named metabolite data file 
raw_input_file <- "MSdata_ST003506_1.txt"

# Which sample sourse to keep for a clean comparison
# This dataset mixes Blood serum and interstitial fluid
# set to the source you wnt, or null to keep all samples

sample_source_filter <- "Blood serum"

# ----- ML workflow (general/reproducible) ----- 

# Random seed for reproducibility
# Setting a seed means the random train/test split and model gives identical
# results everytime the script is run 
set.seed(42)

# Proportion of data that will be used for training (the rest for testing)
train_proportion <- 0.7

# Number of trees in the Random Forest (more trees = more stable but slower)
n_trees <- 500 

# Number of top important features to display
top_features <- 20

# ===================================
# Section 1: DATA INGESTION (data-specific)
# Parses Metabolomics Workbench named metabolite
# data into a clean matrix + metadata table.
# This section is specific to this file format. 
# ====================================

# ----- Read the raw file ----- 
# sep = \t because its tab-delimited 
# header = false because we will handle manually
#(file has metabolite header row and factors row with group assignments)
raw <- read.delim(raw_input_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)


# ----- Extract the sample IDs -------
# row 1, column 3 onwards
# columns 1-2 are metabolite name and refmet_name

sample_ids <- as.character(raw[1, 3:ncol(raw)])

# ----- Extract the factors (row 2) and parse out the groups ----- 
# Each entry looks like "sample source: blood serum | group:control" 
# we need to pull out the sample source and the group seperately 
factors_row <- as.character(raw[2, 3:ncol(raw)])

# Extract the group, everything after Group:
#sub() finds the pattern and replaces, capturing the part we need 
group_labels <- sub(".*Group:", "", factors_row)

# Extract the sample source 
# text between sample source: and |
source_labels <- sub(".*Sample source:", "", factors_row) # remove everything before and including sample source:
source_labels <-sub("\\|.*", "", source_labels) # remove everything after and including | operator
source_labels <- trimws(source_labels)
#----- Build the metadata table ----- 
meta <- data.frame(
  Sample = sample_ids,
  Group = group_labels,
  Source = source_labels,
  stringsAsFactors = FALSE
)

# -----Extract the metabolite data matrix -----
# Metabolite names are in column 1, row 3+
metabolite_names <- as.character(raw[3:nrow(raw), 1])

# The numeric values are rows 3+ , columns 3+
data_matrix <- raw[3:nrow(raw), 3:ncol(raw)]
data_matrix <- as.data.frame(lapply(data_matrix, as.numeric))
data_matrix <- as.matrix(data_matrix)

# Set rownames (metabolites) and column names (samples)
rownames(data_matrix) <- metabolite_names
colnames(data_matrix) <- sample_ids

# ----- Filter by the sample source if specified in config -----
if(!is.null(sample_source_filter)){
  keep_samples <- meta$Sample[meta$Source == sample_source_filter]
  data_matrix <- data_matrix[, keep_samples]
  meta <- meta[meta$Source== sample_source_filter, ]
  cat(" Filtered to" , sample_source_filter, "samples\n")
}

cat("Samples:", ncol(data_matrix), "\n")
cat("Metabolites", nrow(data_matrix), "\n")
cat("Groups:\n")
print(table(meta$Group))









