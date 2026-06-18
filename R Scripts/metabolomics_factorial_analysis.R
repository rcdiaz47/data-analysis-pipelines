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
factor2_levels <- c("Control", "Treatment")

# ----- Sample Metadata -----
# For each sample, specify its values for both factors
# These must be in the same order as the sample 

sample_names <- c("WT_Ctrl_1", "WT_Ctrl_2", "WT_Trt_1", "WT_Trt_2", 
                  "KO_Ctrl_1", "KO_Ctrl_2", "KO_Trt_1", "KO_Trt_2")

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


















