library(dplyr)
# Visualize Ms.Ms identified and peptide sequences identified
summary <- read.delim("summary.txt", stringsAsFactors = FALSE)
print(summary[,c("Raw.file", "MS.MS.identified....", "Peptide.sequences.identified")])

#Read in maxQuant proteinGroups output for analysis
proteins <- read.delim("proteinGroups.txt", stringsAsFactors = FALSE)
head(proteins)

# Extract the sample columns for testing (Label Free quantification)
sample_cols <- grep("LFQ", colnames(proteins), value = TRUE)

cat("Found", length(sample_cols), "samples:\n")
print(sample_cols)


# Extract as sample columns as a matrix
x <- as.matrix(proteins[, sample_cols])
rownames(x) <- proteins$Majority.protein.IDs
sum(duplicated(rownames(x)))

# Clean the column names for better readability
colnames(x) <- c("KO_5", "KO_6", "KO_7", "KO_8", "WT_1", "WT_2", "WT_3", "WT_4")

# Replace zeroes with NA to not skew PCA analysis 
x[x==0] <- NA

#Log transform the data with log2
x_log <- log2(x)
head(x_log)

# Median Centering normalization

# Calculate the column median for each sample
column_median <- apply(x_log, 2, median, na.rm = TRUE)
cat("Sample medians BEFORE centering:\n")
print(column_median)

# apply median centering to each column
x_norm <- sweep(x_log, 2, column_median, FUN = "-")

# Verify all column medians are now 0
cat("Sample medians after centering:\n")
print(apply(x_norm, 2, median, na.rm = TRUE))


# Boxplot to visually confirm normalization
boxplot(x_norm,
        las = 2,
        main = "Median Centered LFQ Intensities",
        ylab = "Median Centered log2 (LFQ intensities)")


# filter proteins detected in at least 70% of samples (>= 6 out of 8)
detection_rate <- rowMeans(!is.na(x_norm))
keep <- detection_rate >= 0.7
x_filt <- x_norm[keep, ]

cat("Proteins Before Filtering:", nrow(x_norm), "\n")
cat("Proteins After Filtering:", nrow(x_filt), "\n")
cat("Proteins removed:", sum(!keep), "\n")


# Define sample metadata 
meta <- data.frame(
  Sample = colnames(x_filt),
  Group = c("KO", "KO", "KO", "KO", "WT", "WT", "WT", "WT")
)

# Run PCA (Transpose so samples are rows, proteins are columns)

#Remove zero variance proteins before PCA
x_complete <- na.omit(x_filt)

# Check the variance of each protein
protein_var <- apply(x_complete, 1, var)

# Keep only the proteins with variance > 0 for PCA
x_complete <- x_complete[protein_var > 0, ]

cat("Proteins used for PCA:", nrow(x_complete))


# Transpose and run PCA 
x_pca <- t(x_complete)
pca_res <- prcomp(x_pca, scale. = TRUE)

#Calculate the variance explained 
variance_explained <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100

# Build a PCA dataframe for plotting
pca_data <- data.frame(
  PC1 = pca_res$x[ , 1],
  PC2 = pca_res$x[ , 2],
  Sample = rownames(pca_res$x),
  Group = meta$Group
)


library(ggplot2)

ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point( size = 4) +
  geom_text( vjust = -1, show.legend = FALSE) +
  labs(
    title = "PCA Proteomics LFQ WT vs KO",
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  theme_minimal() +
  theme(legend.position = "right")


# Define our groups 
ko_cols <- c("KO_5", "KO_6", "KO_7", "KO_8")
wt_cols <- c("WT_1", "WT_2", "WT_3", "WT_4")

results <- data.frame(
  Protein = rownames(x_filt),
  log2FC = NA,
  pvalue = NA
)


for( i in 1:nrow(x_filt)){
  
  ko_vals <- as.numeric(x_filt[i, ko_cols])
  wt_vals <- as.numeric(x_filt[i, wt_cols])
  
  # Need at least 2 non NA values per group to run t test
  if (sum(!is.na(ko_vals)) >= 2 & sum(!is.na(wt_vals)) >= 2 &
    var(ko_vals, na.rm = TRUE) > 0 & var(wt_vals, na.rm = TRUE) > 0)
  {
    t_res <- t.test(ko_vals, wt_vals)
    results$pvalue[i] <- t_res$p.value
    results$log2FC[i] <- mean(ko_vals, na.rm = TRUE) - mean(wt_vals, na.rm = TRUE)
  }
}

# Multiple testing correction
results$padj <- p.adjust(results$pvalue, method = "BH")

cat("Proteins tested:", sum(!is.na(results$pvalue)), "\n")
cat("Significant (FDR < 0.05):", sum(results$padj < 0.05, na.rm = TRUE), "\n")
cat("Significant (FDR < 0.1):", sum(results$padj < 0.1, na.rm = TRUE), "\n")



## Only 13 signficant proteins, WT1 and KO7 may be contributing to this

# lets check the raw p value distribution
hist(
  results$pvalue,
  breaks = 50,
  main ="Raw p values",
  xlab = "p val",
  col = "steelblue" 
)

# Check significant proteins before correction
cat("Raw p < 0.05", sum(results$pvalue < 0.05, na.rm = TRUE))
cat("Raw p < 0.01", sum(results$pvalue < 0.01, na.rm = TRUE))

# Volcano plot to view statistically significant proteins 

# Add significane labels 
results$significance <- "Not Significant"
results$significance[results$padj < 0.1 & results$log2FC >= 0.58] <- "Up in KO"
results$significance[results$padj < 0.1 & results$log2FC <= - 0.58] <- "Down in KO"

table(results$significance)

# Volcano Plot
ggplot(results, aes(x = log2FC, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Not Significant" = "grey60", "Up in KO" = "red", "Down in KO" = "blue")) +
  geom_vline(xintercept = c(-0.58,0.58), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "black") +
  
  geom_text_repel(data = top_labeled,
                  aes(label = Protein),
                  size = 3,
                  max.overlaps = 20,
                  box.padding = 0.3,
                  show.legend = FALSE) +
  
  labs(
    title = "Volcano Plot KO vs WT",
    x = "log2 Fold change (KO/WT)",
    y = "-log10 (p-value)",
    color = "Significance"
  ) + theme_minimal() + theme(legend.position = "right")

# Get the top 20 up and top 20 down by adj p value

library(ggrepel)

top_up <- results[results$significance == "Up in KO", ]
top_up <- top_up[order(top_up$padj), ][1:min(20, nrow(top_up)), ]

top_down<- results[results$significance == "Down in KO", ]
top_down <- top_down[order(top_down$padj), ][1:min(20, nrow(top_down)), ]

top_labeled <- rbind(top_up, top_down)


library(pheatmap)

# Get the significant proteins
sig_proteins <- results$Protein[results$significance%in% c("Up in KO", "Down in KO")]

# Subset the normalized data to significant proteins
heatmap_data <- x_filt[rownames(x_filt) %in% sig_proteins, ]

# Must remove NAs since we are using "Correlation" distance as our argument
heatmap_data_complete <- na.omit(heatmap_data)

cat("Proteins before NA removal:", nrow(heatmap_data), "\n")
cat("Proteins after NA removal", nrow(heatmap_data_complete), "\n")




# Create annotation for columns (samples)
annotation_col <- data.frame(
  Group = c("KO", "KO","KO", "KO", "WT", "WT", "WT", "WT"),
  row.names = colnames(heatmap_data)
)

# Plot the heatmap
pheatmap(
  heatmap_data_complete,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "average",
  fontsize_row = 7,
  fontsize_col = 7,
  treeheight_row = 100,
  treeheight_col = 100,
  main = "Significant Proteins KO vs WT",
  filename = "heatmap_proteomics_LFQ.pdf",
  width = 14,
  height = 12
)


# Pathway enrichment analysis with cluster profiler
# We must first convert our Uniprot Ids to Entrez Gene IDs
# Sig proteins we created earlier 


library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(AnnotationDbi)


# Clean up our protein ids
# Some may have multiple ids seperated by semicolons e.g. "P12334;P123342"
# Take the first one

sig_proteins_clean <- gsub(";,", "", sig_proteins)

cat("Significant proteins before cleaning:", length(sig_proteins), "\n")
cat("Significant proteins after cleaning:", length(sig_proteins_clean), "\n")

# Convert the Uniprot Ids to the Entrez Gene Ids 
entrez_ids <- mapIds(
  org.Mm.eg.db,
  keys = sig_proteins_clean,
  column = "ENTREZID",
  keytype = "UNIPROT",
  multiVals = "first"
)

# Remove NAS (proteins that couldnt be mapped)
entrez_ids <- entrez_ids[!is.na(entrez_ids)]


#HELLO




