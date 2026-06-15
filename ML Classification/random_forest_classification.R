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

# Which sample source to keep for a clean comparison
# This dataset mixes Blood serum and interstitial fluid
# set to the source you want, or null to keep all samples

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

# ======================================
# SECTION 2: ML WORKFLOW (general/reproducible)
# Prepares the data, trains a Random Forest classifier
# Evaluates performance, and extracts feature importance
# Works on any standard matrix + metadata 
# ========================================

# ----- Prepare the data for Random Forest ----- 
# randomForest expects samples as rows and features as columns
# right now we have metabolites as rows, so we transpose it 
x_ml <- t(data_matrix)

# ----- Handle the missing values ----- 
# randomForest cannot handle NAs 
# check how many missing values we have

cat("Missing values:", sum(is.na(x_ml)), "\n")

# Remove the metabolites that are entirely NA
x_ml <- x_ml[ , colSums(is.na(x_ml)) < nrow(x_ml)]


# Simple approach: impute NAs with the column (metabolite) median
# For each metabolite, replace any NA with that metabolites median value 

for(i in 1:ncol(x_ml)){
  na_positions <- is.na(x_ml[, i])
  if(any(na_positions)){
    x_ml[na_positions, i] <- median(x_ml[, i], na.rm = TRUE)
  }
}

cat("Missing values after imputation", sum(is.na(x_ml)), "\n")

# ----- Attach the group labels as the classification target ----- 
# The Group column from meta becomes what we are predicting 
# Must be a factor for randomForest to do classification (not regression)
# Make sure the sample order in meta matches the rows of x_ml

y <- factor(meta$Group[match(rownames(x_ml), meta$Sample)])

cat("Class distribution:\n")
print(table(y))


# ----- Train/Test split (stratified) ----- 
# Split data into training (build the model) and test (evaluate it) set
# Stratified = preserves the class proportions in both sets
# Important for imbalanced data so we dont end up with too few of one class
# train_proportion comes from config 

# createDataPartition from caret does stratified splitting automatically
# It returns the row indices to use for training 
train_index <- createDataPartition(y, p = train_proportion, list = FALSE)

# Split the features 
x_train <- x_ml[train_index,]
x_test <- x_ml[-train_index,]

# Split the labels
y_train <- y[train_index]
y_test <- y[-train_index]

# ----- Confirm the split ----- 
cat("Training samples:", nrow(x_train), "\n")
cat("Test samples:", nrow(x_test), "\n")

cat("\nTraining Class Distribution:\n")
print(table(y_train))

cat("\nTest Class Distribution:\n")
print(table(y_test))

#----- Train the random forest model -----
# Builds the classifier using ONLY the training data
# The model learns patterns in the metabolites values that distinguish the groups
# n_trees comes from the user config section

rf_model <- randomForest(
  x = x_train, # training features (metabolite values)
  y = y_train, # training groups (control/lymphedema)
  ntree = n_trees, # number of trees in the forest 
  importance = TRUE # calculate feature importance so we can extract it later
)

# ----- View the model summary ----- 
# Prints the model details including out-of-bag error (OOB) estimate
print(rf_model)

# ----- Extract feature importance ----- 
# Random Forest tracks  how much each metabolite contributed to classifications
# MeanDecreasedGini: how much each feature improves node purity across all trees
# Higher values = more important for distinguishing the groups 
# top_features comes from user config section

# Pull importance scores into dataframe
importance_df <- as.data.frame(importance(rf_model))
importance_df$Metabolite <- rownames(importance_df)


# Sort by MeanDecreasedGini and take the top N
top_importance <- importance_df %>%
  arrange(desc(MeanDecreaseGini)) %>%
  slice_head(n=top_features)

# ----- Plot the feature importance ----- 
pdf("rf_feature_importance.pdf", width = 8, height = 6)
print(
  ggplot(top_importance, aes(x = reorder(Metabolite, MeanDecreasedGini), y = MeanDecreasedGini)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(
      title = paste0("Top", top_features, "Metabolites by Importance"),
      subtitle = "Rabdom Forest - Control vs Lymphedema",
      x = "",
      y = "Mean Decrease in Gini"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 8)
    )
)
 dev.off()

# Print the top features to console
 cat("\nTop", top_features, "most important metabolites:\n")
 print(top_importance[, c("Metabolite", "MeanDecreaseGini")])









