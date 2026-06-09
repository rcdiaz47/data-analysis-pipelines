library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(enrichplot)
library(AnnotationDbi)

#===========================================
# User Configuration- Edit this section only
#===========================================

# Input file name
input_file <- "proteinGroups.txt"

# Number of groups in experiment 
# 2 = T-test only
# 3+ = ANOVA + Tukey + pairwise comparisons 
n_groups <- 2

# Group names must match how they appear in your samples
group_names <- c("KO", "WT")

# Sample names in order - these replace the raw LFQ column names
# Must be in the same order as the LFQ columns appear in the data
sample_names <- c("KO_5", "KO_6", "KO_7", "KO_8", "WT_1", "WT_2", "WT_3", "WT_4")

# Group assignments for each sample (same order as sample_names)
sample_groups <- c("KO", "KO", "KO", "KO", "WT", "WT", "WT", "WT")

# Intensity column pattern to search for
# Maxquant LFQ = "LFQ", or use "Intensity" for raw intensities 
# Change this to match your datas intensity column naming
intensity_pattern <- "LFQ"

# Output directory
output_dir <- "./"

# Thresholds 
fdr_threshold <- 0.1 # FDR significance cutoff 
logfc_threshold <- 0.58 # log2 fold-change cutoff (0.58 = 1.5x fold change)

# Detection filter - minimum fraction of samples a protein must appear in 
detection_threshold <- 0.7

# Organism database for pathway analysis
# org.Mm.eg.db = mouse, org.Hs.eg.db = human
# (loaded above- change the library() call too if not mouse)

#===========================================================
# END OF USER CONFIGURATION
#===========================================================


# ----- Read in proteomics data -----\
proteins <- read.delim(input_file, stringsAsFactors = FALSE)
head(proteins)

# ----- Extract the intensity columns ----- 
sample_cols <- grep(intensity_pattern, colnames(proteins), value = TRUE)

cat("Found", length(sample_cols), "samples:\n")
print(sample_cols)

#----- Build the numeric matrix -----
x <- as.matrix(proteins[, sample_cols])
rownames(x) <- proteins$Majority.protein.IDs

# Check for duplicate protein IDs
cat("Duplicate protein IDs:", sum(duplicated(rownames(x))), "\n")

# ----- Rename columns using sample_names from config ----- # 
# Safety check - number of sample_names must match the number of detected columns 
if(length(sample_names) != length(sample_cols)){
  stop("Number of sample_names in config (", length(sample_names),") does not match number of detected intensity columns (", length(sample_cols), ")")
}

colnames(x) <- sample_names


#----- Replace zeroes with NA to not skew PCA analysis -----
x[x==0] <- NA

#-----Log transform -----
x_log <- log2(x)
head(x_log)

# ----- Median Centering normalization ----- 
column_median <- apply(x_log, 2 , median, na.rm = TRUE)
x_norm <- sweep(x_log, 2, column_median, FUN = "-")

# Verify medians are centered at 0
cat("Sample medians after centering:\n")
print(round(apply(x_norm, 2, median, na.rm = TRUE)))



# ----- Boxplot to confirm normalization ----- 
pdf("normalization_boxplot.pdf", width = 10, height = 6)
boxplot(x_norm,
        las = 2,
        main = "Median Centered Intensities",
        ylab = "Median Centered log2 Intensities")
dev.off()

# -----  Filter proteins by detection rate ----- 
# detection threshold from config (e.g. 0.7 = present in 70% of samples)
detection_rate <- rowMeans(!is.na(x_norm))
keep <- detection_rate >= detection_threshold
x_filt <- x_norm[keep, ]

cat("Proteins Before Filtering:", nrow(x_norm), "\n")
cat("Proteins After Filtering:", nrow(x_filt), "\n")
cat("Proteins removed:", sum(!keep), "\n")


# ----- Build metadata dynamically -----  
meta <- data.frame(
  Sample = sample_names,
  Group = sample_groups,
  stringsAsFactors = FALSE
)

meta$Group <- factor(meta$Group, levels = group_names)

# ----- PCA -----

# Remove proteins with any missing values (PCA cant handle NAs)
x_complete <- na.omit(x_filt)

# Remove zero-variance proteins (they break scaling)
protein_var <- apply(x_complete, 1, var)
x_complete <- x_complete[protein_var > 0, ]

cat("Proteins used for PCA:", nrow(x_complete), "\n")

# Transpose so samples are rows, then run PCA 
x_pca <- t(x_complete)
pca_res <- prcomp(x_pca, scale. = TRUE)

#Calculate the variance explained by each PC 
variance_explained <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100

# Build a PCA data frame for plotting
pca_data <- data.frame(
  PC1 = pca_res$x[ , 1],
  PC2 = pca_res$x[ , 2],
  Sample = rownames(pca_res$x),
  Group = meta$Group
)

pdf("pca_plot.pdf", width = 8, height = 6)
print(
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point( size = 4) +
  geom_text( vjust = -1, show.legend = FALSE) +
  labs(
    title = "PCA Proteomics",
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  theme_minimal() +
  theme(legend.position = "right")
)
dev.off()

# ----- Reusable function for pairwise comparison between two groups -----
# Takes two group names, the data matrix, and metadata
# Returns a dataframe with log2FC, pvalue, and padj for every protein

run_pairwise_proteomics <- function(group1, group2, x_mat, meta){
  
  # Get the sample names belonging to each group
  group1_samples <- meta$Sample[meta$Group == group1]
  group2_samples <- meta$Sample[meta$Group == group2]
  
  # Build a results dataframe, one row per protein
  res <- data.frame(
    feature_id = rownames(x_mat),
    log2FC = NA,
    pvalue = NA
  )
  
  # Loop through every protein
  for(i in 1:nrow(x_mat)){
    
    # Extract values for each group for this protein
    g1_vals <- as.numeric(x_mat[i, group1_samples])
    g2_vals <- as.numeric(x_mat[i, group2_samples])
    
    # Safety check, need at least 2 non-NA values per group
    # and non-zero variance in both groups to run a t-test
    if(sum(!is.na(g1_vals)) >= 2 & sum(!is.na(g2_vals)) >= 2 &
       var(g1_vals, na.rm =TRUE) > 0 & var(g2_vals, na.rm = TRUE) > 0){
      
      t_res <- t.test(g1_vals, g2_vals)
      res$pvalue[i] <- t_res$p.value
      
      # log2FC = mean of group1 minus mean of group 2
      res$log2FC[i] <- mean(g1_vals, na.rm = TRUE) - mean(g2_vals, na.rm = TRUE)
      
    }
  }
  
  # FDR correction across all proteins
  res$padj <- p.adjust(res$pvalue, method = "BH")
  
  res
}

# ----- Statistical Testing -----
# t-test for 2 groups, ANOVA + Tukey for 3+ groups
# n_groups and group_names come from config section

group <- factor(meta$Group)

if(n_groups == 2){
  
  # For 2 groups, run a single pairwise t-test comparison
  res_pairwise_list <- list(
    run_pairwise_proteomics(group_names[1], group_names[2], x_filt, meta)
  )
  
  # Name the comparison for downstream reference
  names(res_pairwise_list) <- paste(group_names[1], "vs", group_names[2])
  
} else{
  
  # ----- One-way ANOVA per protein ----- 
  # apply() runs the function on every row (protein) of x_filt
  anova_res <- apply(x_filt, 1, function(z){
    
    df <- data.frame(value = z, group = group)
    df <- df[!is.na(df$value), ]
    
    # Safety check - need enough data and at least 2 groups present
    if(nrow(df) < 3 || length(unique(df$group)) < 2){
      
      return(data.frame(pvalue = NA))
      
    }
    
    # Fit ANOVA and extract p-value
    fit <- aov(value ~ group, data = df)
    p <- summary(fit)[[1]][["Pr(>F)"]][1]
    
    data.frame(pvalue = p)
    
  })
  
  # Combine results and add protein IDs
  anova_res <- do.call(rbind, anova_res)
  anova_res$feature_id <- rownames(x_filt)
  rownames(anova_res) <- NULL
  
  # FDR correction 
  anova_res$padj <- p.adjust(anova_res$pvalue, method = "BH")
  
  # Significant proteins only 
  anova_sig_features <- anova_res$feature_id[anova_res$padj < fdr_threshold]
  anova_significant_features_clean <- anova_sig_features[!is.na(anova_sig_features)]
  
  #----- Tukey HSD Post Hoc ----- 
  # ANOVA says something differs; Tukey says which group differs 
  tukey_list <- lapply(anova_significant_features_clean, function(fid){
    z <- x_filt[fid, ]
    df <- data.frame(value = z, group = group)
    df <- df[!is.na(df$value), ]
    s
    fit <- aov(value ~ group, data = df)
    tk <- TukeyHSD(fit)
    
    out <- as.data.frame(tk$group)
    out$comparison <- rownames(out)
    out$feature_id <- fid
    rownames(out) <- NULL
    out
  })
  
  tukey_res <- do.call(rbind, tukey_list)
  tukey_sig <- tukey_res %>% filter(`p adj` < fdr_threshold)
  
  # Generate all pairwise combination automatically ----- 
  pairs <- combn(group_names, 2, simplify = FALSE)
  
  res_pairwise_list <- lappy(pairs, function(pair){
    
    run_pairwise_proteomics(pair[1], pair[2], x_filt, meta)
    
  })
  
  names(res_pairwise_list) <- sapply(pairs, function(pair){
    
    paste(pair[1], "vs", pair[2])
  })
  
}


# ----- Volcano Plot -----

# ----- Add significance labels to each pairwise comparison ----- 
# Loops through every comparison and labels each protein
# Up = higher in the first group, Down = lower in the first group 
# Thresholds come from congif (fdr_threshold, logfc_threshold)

res_pairwise_list <- lapply(res_pairwise_list, function(res){
  res$significance <- "Not Significant"
  res$significance[res$padj < fdr_threshold & res$log2FC >= logfc_threshold] <- "Up"
  res$significance[res$padj < fdr_threshold & res$log2FC <= -logfc_threshold] <- "Down"
  res
})


#----- Reusable volcano plot function ----- 
# Builds a volcano plot for one comparison
# Lables the top proteins by adjust p-value
# comparison_name used for title and axis labels
plot_volcano_proteomics <- function(res, comparison_name, top_n = 30){
  
  # Split the comparison name into the two group names for labeling
  groups <- strsplit(comparison_name, " vs ")[[1]]
  group1 <- groups[1]
  group2 <- groups[2]
  
  # Get top up and top down protein by adjusted p value for labeling
  top_up <- res[res$significance == "Up", ]
  top_up <- top_up[order(top_up$padj), ][1:min(top_n, nrow(top_up)), ]
  
  top_down <- res[res$significance == "Down",]
  top_down <- top_down[order(top_down$padj), ][1:min(top_n, nrow(top_down)), ]
  
  top_labeled <- rbind(top_up, top_down)
  
  # Build the plot 
  ggplot(res, aes(x = log2FC, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Not Significant" = "grey60", 
                                 "Up" = "red4",
                                 "Down" = "blue4"),
                       labels = c(
                         "Up" = paste("Up in", group1),
                         "Down" = paste("Down in", group1),
                         "Not Significant" = "Not Signficant"
                       )) +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed", color = "black")+
    geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "grey70") +
    geom_text_repel(
      data = top_labeled,
      aes(label = feature_id),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.3, 
      show.legend = FALSE
    ) +
    labs(
      title = paste("Volcano Plot:", comparison_name),
      x = paste0("log2 Fold Change (", group1, "/", group2, ")"),
      y = "-log10 (FDR)",
      color = "Significance"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  
}

# ----- Generate volcano plots for each comparison ----- 
for(comparison_name in names(res_pairwise_list)){
  res <- res_pairwise_list[[comparison_name]]
  filename <- paste0("volcano_", gsub(" ", "_", comparison_name), ".pdf")
  
  pdf(filename, width = 8, height = 6)
  print(plot_volcano_proteomics(res, comparison_name))
  dev.off()
}


#----- Build annotation for heatmap columns (samples) -----
# Maps each sample to its group for the colored annotation bar
# Built from metadata so it works for any number of groups 

annotation_col <- data.frame(
  Group = meta$Group,
  row.names = meta$Sample 
)

# ----- Generate heatmaps for each pairwise comparison -----
# Loops through res_pairwise_list and creates one heatmap per comparison 
# Shows only significant proteins for that comparison 
# File name automatically generated based on comparison_name 

for(comparison_name in names(res_pairwise_list)){
  
  res <- res_pairwise_list[[comparison_name]]
  
  # Get significant proteins for this comparison (Up and Down)
  sig_proteins <- res$feature_id[res$significance %in% c("Up", "Down")]
  
  # Subset normalized data to just these proteins 
  heatmap_data <- x_filt[rownames(x_filt) %in% sig_proteins, ]
  
  # Remove proteins with NAs (correlation discant cannot handle them)
  heatmap_data_complete <- na.omit(heatmap_data)
  
  # Skip this comparison if too few proteins are present
  if(nrow(heatmap_data_complete) < 2){
    cat("Not enough significant proteins present for:", comparison_name, "-skipping heatmap\n")
    next
  }
  
  cat(comparison_name, "-", nrow(heatmap_data_complete), "proteins in heatmap\n")
  
  # Create the filename automatically 
  filename <- paste0("heatmap_", gsub(" ", "_", comparison_name), ".pdf")
  
  # Plot and save directly via pheatmaps filename argument 
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
    main = paste("Significant Proteins:", comparison_name),
    filename = filename,
    width = 14,
    height = 12
  )
  
}


# Organism for pathway analysis
# "mouse" or "human"
organism <- "mouse"

# Enrichment p-value cutoff
enrichment_pvalue <- 0.05

# ----- Load the correct organism database based on conifg -----
if(organism == "mouse"){
  org_db <- org.Mm.eg.db
  kegg_organism <- "mmu"
} else if(organism == "human"){
  org_db <- org.Hs.eg.db
  kegg_organism <= "hsa"
} else{
  stop("organism must be either mouse or human")
}


#----- Pathway Enrichment Analysis (GO + KEGG) ----- 
# Loops through every comparison and runs enrichment on significant proteins
# Converts Uniprot IDs to Entrez IDs first (required by clusterProfiler)
# Organism and database come from config (org_db, kegg_organism)

for(comparison_name in names(res_pairwise_list)){
  res <- res_pairwise_list[[comparison_name]]
  clean_comparison <- gsub(" ", "_", comparison_name)
  
  cat("Pathway analysis for:", comparison_name, "===\n")
  
  # Get significant proteins for this comparison
  sig_proteins <- res$feature_id[res$significance %in% c("Up","Down")]
  
  # Skip if too few significant proteins to bother with enrichment
  if(length(sig_proteins) < 5){
    cat("Fewer than 5 significant proteins - skipping enrichment \n")
    next
  }
  
  # ----- Clean protein IDs ----- 
  # Maxquant IDs may have multiple IDs seperated by semi colons (e.g. "P12345;P67890)
  # gsub(";.*", "", ...) keeps everything before the first semicolon
  sig_proteins_clean <- gsub(";.*", "", sig_proteins)
  
  #----- Convert UniProt Ids to Entrez Gene Ids ----- 
  # clusterProfiler requires Entrez IDs for enrichment 
  entrez_ids <- mapIds(
    org_db,
    keys = sig_proteins_clean,
    column = "ENTREZID",
    keytype = "UNIPROT",
    multiVals = "first"
  )
  
  # Remove proteins that couldnt be mapped (NA)
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  
  cat("Mapped", length(entrez_ids), "of", length(sig_proteins_clean), "proteins to Entrez IDs\n")
  
  # Skip if mapping failed for too many 
  if(length(entrez_ids) < 5) {
    cat("Too few mapped IDs - skipping enrichment")
    next
  }
    
  
  #----- GO Enrichment (Biological Process) ----- 
  # enrichGo tests if significant proteins are enriched for biological processes
  go_results <- enrichGO(
    gene = entrez_ids,
    OrgDb = org_db,
    ont = "BP", # BP = Biological Process (also "MF", "CC", or "ALL"
    pAdjustMethod = "BH",
    pvalueCutoff = enrichment_pvalue,
    readable = TRUE # Converts Entrez Ids back into gene symbols in output
    )
  # Save GO results if any were found 
  if(!is.null(go_results) && nrow(go_results) > 0){
    write.csv(as.data.frame(go_results), 
              paste0("GO_enrichment_", clean_comparison, ".csv"), row.names = FALSE)
    
    pdf(paste0("GO_dotplot_", clean_comparison, ".pdf"), width = 10, height = 8)
    print(dotplot(go_results, showCategory = 20, title = paste("GO BP:", comparison_name)))
    dev.off()
    
    cat("GO Enrichment:", nrow(go_results), "significant terms\n")
  } else{
    cat("No significant GO terms found\n")
  }
  
  #----- KEGG pathway enrichment ----- 
  # enrichKEGG tests if significant proteins are enriched for KEGG pathways 
  kegg_results <- enrichKEGG(
    gene = entrez_ids,
    organism = kegg_organism,
    pAdjustMethod = "BH",
    pvalueCutoff = enrichment_pvalue
  )
  
  # Save KEGG results if any were found 
  if(!is.null(kegg_results) && nrow(kegg_results) > 0){
    write.csv(as.data.frame(kegg_results), 
              paste0("KEGG_enrichment_", clean_comparison, ".csv"),
              row.names = FALSE)
    
    pdf(paste0("KEGG_dotplot_", clean_comparison, ".pdf"), width = 10, height = 9)
    print(dotplot(kegg_results, showCategory = 20, title = paste("KEGG:", comparison_name)))
    dev.off()
    
    cat("KEGG Enrichment:", nrow(kegg_results), "significant pathways\n")
    
  }else{
    cat("No significant KEGG pathways found\n")
  }
  

}




















