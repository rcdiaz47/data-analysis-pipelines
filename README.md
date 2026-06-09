# Data Analysis Pipelines
A collection of reproducible R pipelines for downstream analysis of omics data. 

## Metabolomics Downstream Analysis Pipeline
A flexible, end-to-end R pipeline for downstream analysis of LC-MS/MS metabolomics data (Compound Discoverer output). This pipeline dynamically adapts to the number of experiment groups, automatically running the appropriate statistical tests and generating publication-quality visualizations and result tables.

## Features
- **Dynamic group handling:** automatically runs a t-test for 2 groups, or ANOVA with Tukey HSD post-hoc for 3 or more groups.
- **Automatic pairwise comparisons:** generates every possible group comparison without manual hardcoding.
- **Reproducibile preprocessing:** log2 transoformation, missing-value filtering, and median-centering normalization.
- **Publication-quality figures:** PCA plots, volcano plots, fold-change bar plots, and clustered heatmaps
- **Organized output:** exports all the relevant results as clearly names CSV files for downstream use (including MetaboAnalyst pathway input)
- **Configurable thresholds:** FDR and log2 fold-change cutoffs set in one place

## How to use
All settings are controlled from a single configuration section at the top of the script. Edit only this section:

'input_file <- "your_data.xlsx" # Compound discoverer result file
n_groups <- 3 #2 = t test, 3+ = Anova + Tukey
group_names <- c("WT", "WRN", "shPGC") # must match sample column names
output_dir <- "./" # Location to output result files 
fdr_threshold <- 0.05 # Significane cutoff
logfc_threshold <- 1 # log2 fold-change cutoff for volcano plots'

Then run the full script. Group names are matched against the sample column headers automatically, so they must appear in the column names exactly as written.

## Workflow 
1. Preprocessing - reads the input file, cleans metabolite names, builds the numeric matrix, log2-transforms, filters metabolites missing in more than 30% of samples, and normalizes by median-centering.
2. Exploratory analysis - density and boxplots before/after normalization, plus a PCA plot to visualize group seperation.
3. Statistical testing - t-test (2 groups) of ANOVA + Tukey HSD (3+ groups), with Benjamini Hochberg FDR correction.
4. Pairwise comparisons - automatic t-tests for every group pair, with significane labeling,
5. Visualization - volcano plots, fold-change bar plots, and heatmaps generated automatically for each comparison.
6. Export - all results written to CSV, plus a printed run summary.

## Outputs 
Normalized and raw data matrices (CSV)
Sample metadata (CSV)
ANOVA and Tukey results (CSV, 3+ groups)
Full and cleaned pairwise results per comparison (CSV)
Significant metabolite lists formatted for MetaboAnalyst pathway analysis input (CSV)
Summary statistics (CSV)
PCA, volcano, bar plot, and heatmap figures (PDF)

## Requirements 
R with the following packages:
tidyverse, readxl, pheatmap, ggplot2, dplyr, tidyr, patchwork, ggrepel, vegan

## Proteomics Analysis Pipeline 
Note: This proteomics pipeline is actively under development, Core functionality is in progress - see the metabolomics pipeline above for completed work. 

