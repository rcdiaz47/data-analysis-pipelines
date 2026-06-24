# ===================================================
# Metabolomics Factorial Analysis Pipeline
# Two-way ANOVA with interaction effects for experimental design with 2 factors 
# ===================================================

#----- Libraries -----
library(tidyverse)
library(readxl)
library(pheatmap)
library(ggplot2)

# ===================================================
# USER CONFIGURATION - Edit this section only 
# ===================================================

# Input file (Compound Discoverer output)
input_file <- "your_data.xlsx"

# ----- Factor definitions -----
# A factorial design has 2 or more factors
# Example: Factor1= Genotype (WT/KO), Factor2 = Treatment (Control/Treatment)

# Names of the two factors
factor1_name <- "Genotype"
factor2_name <- "Treatment"

# The levels within each factor
factor1_levels <- c("WT", "KO")
factor2_levels <- c("Control", "Treated")

# ----- Sample Metadata -----
# For each sample, specify its values for both factors
# These must be in the same order as the sample 

sample_names <- c("WT_Control_1", "WT_Control_2", "WT_Control_3", "WT_Control_4", 
                  "KO_Control_1", "KO_Control_2", "KO_Control_3", "KO_Control_4", 
                  "WT_Treated_1", "WT_Treated_2", "WT_Treated_3", "WT_Treated_4",
                  "KO_Treated_1", "KO_Treated_2", "KO_Treated_3", "KO_Treated_4")

sample_factor1 <- c("WT", "WT", "WT", "WT", "KO", "KO", "KO", "KO")
sample_factor2 <- c("Control", "Control", "Treated", "Treated",
                    "Control", "Control", "Treated", "Treated")

# ----- Output and Thresholds ----- 
output_dir <- "./"
fdr_threshold <- 0.05
logfc_threshold <- 1

# =====================================
# END USER CONFIGURATION
# =====================================

# =======================================
# Simulated data for testing purposes only, generates a 2x2 factorial dataset with known 
# effects baked in, so we can verify the 2 way ANOVA correctly detects them
# Replace this section with real data ingestion later. 
# ========================================

set.seed(42) # For reproducibility 

# ----- Simulation settings ------
n_metabolites <- 100
n_per_group <- 5 # samples per factor combination (5x4 groups = 20 samples )

# Build the sample design. 2x2 genotype x treatment.
sim_meta <- expand.grid(
  Replicate = 1:n_per_group,
  Genotype = c("WT", "KO"),
  Treatment = c("Control", "Treated")
  
)

sim_meta$Sample <- paste(sim_meta$Genotype, sim_meta$Treatment, sim_meta$Replicate, sep = "_")


n_samples <- nrow(sim_meta)

# ----- Create the data matrix -----
# Start with baseline random noise for all metabolites/samples 

sim_data <- matrix(
  rnorm(n_metabolites * n_samples, mean = 10, sd = 1),
  nrow = n_metabolites,
  ncol = n_samples
)

rownames(sim_data) <- paste0("Metabolite_", 1:n_metabolites)
colnames(sim_data) <- sim_meta$Sample

# ----- Plant Known effects to verify detection from model ----- 
# Metabolites 1-20: Genotype effect ( KO Higher than WT)
ko_samples <- which(sim_meta$Genotype == "KO")
sim_data[1:20, ko_samples] <- sim_data[1:20, ko_samples] + 3

# Metabolites 21-40: Treatment effect (Treatment higher than Control)
treated_samples <- which(sim_meta$Treatment == "Treated")
sim_data[21:40, treated_samples] <- sim_data[21:40, treated_samples] + 3

# Metabolites 41-60 interaction effect
# effect of treatment depends on genoytype - only KO+treated is elevated 
ko_treated <- which(sim_meta$Genotype == "KO" & sim_meta$Treatment == "Treated")
sim_data[41:60, ko_treated] <- sim_data[41:60, ko_treated] + 4

# Metabolites 60-100 will have no effect. Pure noise, shouldnt be significant

# ----- Output: clean matrix + metadata ----- 
data_matrix <- sim_data
meta <- data.frame(
  Sample = sim_meta$Sample,
  Factor1 = sim_meta$Genotype,
  Factor2 = sim_meta$Treatment,
  stringsAsFactors = FALSE
)

# =============================================
# TWO-WAY ANOVA
# For each metabolite we will test 3 things :
# 1. Main effect of Factor1, Genotype
# 2. Main effect of Factor2, Treatment,
# 3. Interaction between Genotype and Treatment 

# ----- Set up the factors ----- 
# Convert factor columns to R factors so aov() treats them as categorical 
factor1 <- factor(meta$Factor1, levels = factor1_levels)
factor2 <- factor(meta$Factor2, levels = factor2_levels)

# ----- Run Two-Way ANOVA for every metabolite ----- 
# apply() runs the function on every metabolite (row) of the data matrix
# For each metabolite we fit: value ~ Factor1 * Factor2
# The * means "Both main effects and their interaction"
anova_results <- apply(data_matrix, 1 , function(z){
  
  # Build a small dataframe for this one metabolite
  df <- data.frame(
    value = z,
    factor1 = factor1,
    factor2= factor2
  )
  
  # Fit the two-way ANOVA model
  # value ~ factor1 * factor2 expands to :
  # value ~ factor1 + facotr 2 + factor1:factor2
  fit <- aov(value ~ factor1 * factor2, data = df)
  
  # Extract the p-values from the ANOVA table
  # summary(fit)[[1]] is the ANOVA table; "Pr(>F)" is the p value column 
  p_values <- summary(fit)[[1]][["Pr(>F)"]]
  
  
  # p_values will come back in this order
  # [1] factor1 main effect
  # [2] factor2 main effect
  # [3] interaction
  # [4] residuals (NA - ignore)
  data.frame(
    p_factor1 = p_values[1],
    p_factor2 = p_values[2],
    p_interaction = p_values[3]
  )
})

# ----- Combine all the results into one dataframe ----- 
# apply returned a list of small dataframes; bind them into one
anova_df <- do.call(rbind, anova_results)
anova_df$Metabolite <- rownames(data_matrix)
rownames(anova_df) <- NULL


# ----- FDR correction for each effect separately ----- 
# Each effect (factor1, factor2, interaction) gets its own multiple testing correction
anova_df$padj_factor1 <- p.adjust(anova_df$p_factor1, method = "BH")
anova_df$padj_factor2 <- p.adjust(anova_df$p_factor2, method = "BH")
anova_df$padj_interaction <- p.adjust(anova_df$p_interaction, method = "BH")

# ----- Check how many are significant for each effect ----- 
cat("Significant metabolites (FDR <", fdr_threshold, "):\n")
cat(" ", factor1_name, "main effect:", sum(anova_df$padj_factor1 < fdr_threshold, na.rm = TRUE), "\n")
cat(" ", factor2_name, "main effect:", sum(anova_df$padj_factor2 < fdr_threshold, na.rm = TRUE), "\n")
cat(" Interaction:", sum(anova_df$padj_interaction < fdr_threshold, na.rm = TRUE), "\n")






































