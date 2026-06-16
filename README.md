# Data Analysis Pipelines
A collection of reproducible R pipelines for downstream analysis of omics data. 

## Metabolomics Downstream Analysis Pipeline
A flexible, end-to-end R pipeline for downstream analysis of LC-MS/MS metabolomics data (Compound Discoverer output). This pipeline dynamically adapts to the number of experiment groups, automatically running the appropriate statistical tests and generating publication-quality visualizations and result tables.

### Features
- **Dynamic group handling:** automatically runs a t-test for 2 groups, or ANOVA with Tukey HSD post-hoc for 3 or more groups.
- **Automatic pairwise comparisons:** generates every possible group comparison without manual hardcoding.
- **Reproducibile preprocessing:** log2 transoformation, missing-value filtering, and median-centering normalization.
- **Publication-quality figures:** PCA plots, volcano plots, fold-change bar plots, and clustered heatmaps
- **Organized output:** exports all the relevant results as clearly names CSV files for downstream use (including MetaboAnalyst pathway input)
- **Configurable thresholds:** FDR and log2 fold-change cutoffs set in one place

### How to use
All settings are controlled from a single configuration section at the top of the script. Edit only this section:
```r
input_file <- "your_data.xlsx" # Compound discoverer result file
n_groups <- 3 #2 = t test, 3+ = Anova + Tukey
group_names <- c("WT", "WRN", "shPGC") # must match sample column names
output_dir <- "./" # Location to output result files 
fdr_threshold <- 0.05 # Significane cutoff
logfc_threshold <- 1 # log2 fold-change cutoff for volcano plots
```
Then run the full script. Group names are matched against the sample column headers automatically, so they must appear in the column names exactly as written.

### Workflow 
1. Preprocessing - reads the input file, cleans metabolite names, builds the numeric matrix, log2-transforms, filters metabolites missing in more than 30% of samples, and normalizes by median-centering.
2. Exploratory analysis - density and boxplots before/after normalization, plus a PCA plot to visualize group seperation.
3. Statistical testing - t-test (2 groups) of ANOVA + Tukey HSD (3+ groups), with Benjamini Hochberg FDR correction.
4. Pairwise comparisons - automatic t-tests for every group pair, with significane labeling,
5. Visualization - volcano plots, fold-change bar plots, and heatmaps generated automatically for each comparison.
6. Export - all results written to CSV, plus a printed run summary.

### Outputs 
Normalized and raw data matrices (CSV)
Sample metadata (CSV)
ANOVA and Tukey results (CSV, 3+ groups)
Full and cleaned pairwise results per comparison (CSV)
Significant metabolite lists formatted for MetaboAnalyst pathway analysis input (CSV)
Summary statistics (CSV)
PCA, volcano, bar plot, and heatmap figures (PDF)

### Requirements 
R with the following packages:
tidyverse, readxl, pheatmap, ggplot2, dplyr, tidyr, patchwork, ggrepel, vegan

## Proteomics Analysis Pipeline 
A reproducibile R pipeline for downstream analysis of label-free quanitifcation (LFQ) prteomics data. Like the metabolomics pipeline, it dynamically adapts to the number of experimental groups and runs the appropriate statistical test, extending the analysis through to pathway enrichment for biological interpretation.

### Features
- **Dynamic group handling:** t-test for 2 groups, one-way ANOVA with Tukey HSD post-hoc for 3+ groups
- **Automatic pairwise comparison:** generates every group comparison without hardcoding
- **Reproducibile preprocessing:** log2 transformation, median centering normalization, and detection-rate filtering
- **Pathway enrichment:** GO (Biological process) and KEGG pathway enrichment via clusterProfiler, with automatic UniProt-to-Entrez ID mapping
- **Organism support:** configurable for mouse or human
- **Publication-quality figures:** PCA plots, volcano plots, clustered heatmaps, and enrichment dotplots 
- **Configurable thresholds:** FDR, fold-change, and enrichment cutoffs set in one place 

### How to Use
All settings are controlled from the configuration section at the top of the script:
```r
input file <- "proteinGroups.txt" # Maxquant output 
intensity_pattern <- "LFQ" # Intensity column pattern 
n_groups <- 2 # 2 = t-test, 3+ = ANOVA + Tukey
group_names <- c("KO", "WT")
sample_names <- c("KO_5","KO_6","WT_1", "WT_2") # renamed sample columns
sample_groups <- c("KO". "KO", "WT", "WT")
organism <- "mouse" # Mouse or human
fdr_threshold <- 0.1
logfc_threshold <- 0.58
```
### Outputs 
- Normalized data matrix and sample metadata (CSV)
- Statistical results per comparison (CSV)
- ANOVA and Tukey results (3+groups)
- GO and KEGG enrichment results per comparison (CSV)
- PCA, volcano, heatmap, and enrichment dotplot (PDF)

### Requirements
R with the following packages:

```r
dplyr, ggplot2, ggrepel, pheatmap, clusterProfiler, org.Mm.eg.db, org.Hs.eg.db, enrichPlot, AnnotationDbi
```

## Machine Learning: Metabolomics Classification

A Random Forest classification workflow demonstrating supervised machine learning applied to metabolomics data. The pipeline predicts disease status (Control vs Lymphedema) from serum metabolite profiles and identifies the most discriminating metabolites via feature importance.

### Dataset

Public NMR metabolomics data from Metabolomics Workbench study ST003506 (breast cancer related lymphedema vs healthy controls). Data is not inlcuded in this repository and can be downloaded directly from the source (DOI: 10.21228/M8FR6S).

### Workflow 

- **Dataset-specific ingestion** - parses through the Metabolomics Workbench format, extracting group labels, filters to blood serum samples, and outputs a clean matrix and metadata
- **Data preparation** - transposes to samples as rows, removes and empty features, and imputes missing values with the metabolite median
- **Stratified train/test split** - preserves class proportions across training and test sets, important for imbalance data
- **Random Forest training** - trains a classifier with user configurable tree count, and OOB estimation
- **Feature importance** - ranks metabolites by their contribution to classification (Mean Decrease in Gini)
- **Evaluation** - confusion matrix with accuracy, sensitivity, and specificity on the test set that was held out

### Key Result 

The model's feature importance independently identified 5 of the 7 metabolomic biomarkers reported in the original study: 3-methyl-2-oxovalerate, pyruvate, 2-ketoisovalerate, ketoleucine, and tryptophan. This suggests that the model captured biologically meaningful signal. Several of these are branched chain amino acid degradation products, consistent with BCAA altered metabolism. 

### Limitations

With 49 features and 34 samples, the model achieves perfect separation, which reflects the high feature-to-sample ratio rather than a genuinely robust classifier. This project demonstrates an end-to-end ML workflow. Reliable performance claims would need feature selection and validation in an independent cohort. 

### Requirements 

R: with tidyverse, randomForest, caret, pheatmap 











