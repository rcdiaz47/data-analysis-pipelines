# ----- Statistical Analysis ----- #
library(tidyverse)
library(readxl)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggrepel)
library(vegan)



#============================================
# USER CONFIGURATION - Edit this section only
#=============================================

# Input file path

input_file <- "2Group_test_metabolomics.xlsx"

# Number of groups to be tested in experiment
# 2 = t-test only
# 3+ = ANOVAA + Tukey + Pairwise analysis

n_groups <- 2

# Group names must match sample column names
group_names <- c("WT", "WRN")

# Output directory for plots and csv files ("./) for current directory
output_dir <- "./results"

# FDR significant threshold
fdr_threshold <- 0.05

# Log2fold change threshold for plots
logfc_threshold <- 1

# Column containing the metabolite names
name_column <- "Name"


#=========================================
# END USER CONFIGURATION
#=========================================

# Create output directory if it does not exist already
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}

# Helper function for creating output file paths 
output_path <- function(filename){
  
  file.path(output_dir, filename)
  
}

# Validate the input file 

if(!file.exists(input_file)){
  stop(
    paste0(
      "Input file not found: ",
      input_file,
      "\n Please check the input_file path in the user configuration section."
    )
  )
}


# Read input data

cd <- read_xlsx(input_file)

head(cd)

# Validate the metabolite name column 
if(!name_column %in% colnames(cd)){
  stop(
    paste0(
      "The configured metabolite name column '",
      name_column,
      "'was not found in the input file.\n",
      "Avaliable columns are: ",
      paste(colnames(cd), collapse = " ,")
    )
  )
}


# Validate that the name column exists 

# ------ Clean and shorten the metabolite names for readability ------ #
process_metabolite_names <- function(x, max_length = 40){
  
  x %>%
    
    gsub("\\.", " ", .) %>% #Replace dots with spaces
    
    gsub("\\.\\d+$", "", .) %>% #Remove trailing .1, .2, etc. 
    
    gsub("\\s+", " ", .) %>% # Replace multiple spaces with single space
    
    trimws() %>% #Remove leading/trailing whitespace
  
  ifelse(nchar(.) > max_length,
         
         paste0(substr(., 1, max_length - 3), "..."), .)
  
}


cd[[name_column]] <- process_metabolite_names(cd[[name_column]], max_length = 40)


# ----- Build the numeric matrix with the columns we want to use for analyzing ----- #
area_cols <- grep("^Group Area:", colnames(cd), value = TRUE)

# Exclude blanks and QCs from the analysis
area_cols <- area_cols[!grepl("Blank|blank|BLANK|QC|Qc", area_cols)]

# Area columns check
if (length(area_cols) == 0){
  stop(
    paste0(
      "No sample area columns were found.\n",
      "Expected column names beginning with 'Group Area:'."
    )
  )
}


x <- cd %>%
  
  dplyr::select(all_of(area_cols)) %>% as.matrix()

rownames(x) <- cd[[name_column]]

# ----- Replace missing values with Nas ----- #
sum(is.na(x))

x[x== ""] <- NA

# ----- Log transform the data with log2 ------ #
x_log <- log2(x)

# ----- Drop metabolites that are missing in too many samples ----- #
keep <- rowMeans(!is.na(x_log)) >= 0.7

x_filt <- x_log[keep,]


# ----- Data normalization using Median Centering ----- #
column_median <- apply(x_filt,2,median, na.rm = TRUE)

# Subtract each columns median from the value 
x_norm <- sweep(
  
  x_filt, 2, column_median, FUN = "-"
  
)

summary(x_norm[,1])

any(is.na(x_norm))

nrow(x_norm)

# --- Prepare data for plotting; Box and whisker plots and Density plot --- #

# Before normalization 
before_df <- x %>%
  
  as.data.frame() %>% 
  
  tibble::rownames_to_column("metabolite") %>%
  
  pivot_longer(-metabolite, names_to = "sample", values_to = "value")

# After normalization 
after_df <- x_norm %>%
  
  as.data.frame() %>%
  
  tibble::rownames_to_column("metabolite") %>%
  
  pivot_longer(-metabolite, names_to = "sample", values_to = "values")

# --- Density Plots --- #
p1 <- ggplot(before_df, aes(x=value)) + 
  
  geom_density(fill = "blue", alpha = 0.4) +
  
  labs(title ="Before Normalization" , x = "Raw Intensity", y = "Density") +
  
  theme_minimal()

p2 <- ggplot(after_df, aes(x = values)) +
  
  geom_density(fill = "blue", alpha = 0.4) +
  
  labs(title = "After Normalization", x = "Normalized Intensity", y = "Density") +
  
  theme_minimal()

# --- Boxplots (Subset so they dont look all compressed) ---- #

before_df <- x[1:50, ] %>%
  
  as.data.frame() %>%
  
  tibble::rownames_to_column("metabolite") %>%
  
  pivot_longer(-metabolite, names_to = "sample", values_to = "value")

after_df <- x_norm[1:50, ] %>%
  
  as.data.frame() %>%
  
  tibble::rownames_to_column("metabolite") %>%
  
  pivot_longer(-metabolite, names_to = "sample", values_to = "value")


p3 <- ggplot(before_df, aes(x = metabolite, y = value)) +
  
  geom_boxplot(fill = "lightgreen") +
  
  coord_flip() + 
  
  labs(x = "", y = "Before Normalization") + 
  
  theme_minimal() +
  
  theme(axis.text.y = element_text(size = 6))

p4 <- ggplot(after_df, aes(x = metabolite, y = value)) +
  
  geom_boxplot(fill = "lightgreen") +
  
  coord_flip() + 
  
  labs(x = "", y = "After Normalization") + 
  
  theme_minimal() + 
  
  theme(axis.text.y = element_text(size = 6))

# ---- Combine all the plots into one figure ----#
data_boxplots <- (p3|p4)

data_density_plots <- (p1|p2)

pdf("data_boxplots_before_norm_2.pdf", width = 10, height = 6)

plot(data_boxplots)

dev.off()

pdf("data_density_plots_before_norm.pdf", width = 8, height = 6)

plot(data_density_plots)

dev.off()


##------ Exploratory Data Analysis ---------#

# ----- Build metadata  ----- #

if(n_groups == 2){
  meta <- data.frame(
    
    Sample = area_cols,
    
    Group = ifelse(grepl(group_names[1], area_cols), group_names[1], group_names[2]),
    
    stringsAsFactors = FALSE
  )
} else{
  meta <- data.frame(
    
    Sample = area_cols,
    
    Group = sapply(area_cols, function(s){
      
      matched <- group_names[sapply(group_names, function(g) grepl(g,s))]
      
      if (length(matched) == 1) matched else NA
      
    }),
    
    stringsAsFactors = FALSE
  )
}



# ---- PCA plot for all of the metabolites to visualize the group separation ----

# Transpose the data because PCA needs the samples as the rows instead of columns 

x_pca <- t(na.omit(x_norm))

# Run PCA
pca_res <- prcomp(x_pca, scale. = TRUE)

# Extract the PCA scores and add group information to prepare for plotting 
pca_scores <- as.data.frame(pca_res$x)

pca_scores$Group <- meta$Group

# Plot 
pdf(output_path("pca_plot.pdf"), width = 8, height = 6)

ggplot(pca_scores, aes(x = PC1, y=PC2, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill=Group),
               geom = "polygon",
               alpha = 0.2,
               color = NA) +
  labs(
    title = "PCA Scores Plot", 
    x = paste0("PC1 (" , round(summary(pca_res)$importance[2,1] * 100, 1), "%)"),
    y = paste0("PC2 (" , round(summary(pca_res)$importance[2,2] * 100, 1), "%)")
  ) +
  theme_minimal()

dev.off()

## NMDS Exploratory analysis. Skip if you dont have at least 20-30 samples 

## We need samples as rows and features as columns
#x_nmds <- t(na.omit(x_norm))

## Calculate the distance matrix using Bray-Curtis 
#dist_matrix <- vegdist(x_nmds, method = "bray")

## Run NMDS with 2 dimensions
#set.seed(123) # Reproducibility
#nmds_res <- metaMDS(dist_matrix, k = 2, trymax = 100)

## Check the stress value of the results (<0.2 is acceptable, < 0.1 is good)
#cat("NMDS Stress:", nmds_res$stress, "\n")

## Extract the scors and add group information
#nmds_scores <- as.data.frame(scores(nmds_res))
#nmds_scores$Group <- meta$Group

## Plot the NMDS
#ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Group)) +
 # geom_point(size = 3) + 
 # stat_ellipse(aes(fill= Group),
  #             geom = "polygon", 
  #             alpha = .2,
   #            level = 0.95,
   #            color = NA) +
  #labs(
   # title = "NMDS Ordination",
   # subtitle = paste0("Stress = ", round(nmds_res$stress, 3)) ,
   # x = "NMDS1",
   # y = "NMDS2"
 # ) + theme_minimal()



# ----- Statistical Testing ----- #

# ---- Reusable function for pairwise analysis of groups ---- #

run_pairwise <- function(group1, group2, x_mat, meta){
  
  samples_keep <- meta$Sample[meta$Group %in% c(group1, group2)]
  
  groups_keep <- meta$Group[meta$Group %in% c(group1, group2)]
  
  x_sub <- x_mat[,samples_keep]
  
  res <- apply(x_sub, 1, function(z){
    
    df <- data.frame(
      value = z,
      group = groups_keep
    )
    
    df <- df[!is.na(df$value), ]
    
    if(length(unique(df$group)) < 2){
      return(c(logFC = NA, pvalue = NA))
    }
    
    # LogFC (difference in means)
    m1 <- mean(df$value[df$group == group1], na.rm = TRUE)
    
    m2 <- mean(df$value[df$group == group2], na.rm = TRUE)
    
    logFC = m1 - m2
    
    # t test
    
    p <- t.test(value ~ group, data = df)$p.value
    
    c(logFC = logFC, pvalue = p)
    
    
  })
  
  res <- as.data.frame((t(res)))
  
  res$feature_id <- rownames(x_mat)
  
  res$padj <- p.adjust(res$pvalue, method = "BH")
  
  res
  
}


# Runs T test for 2 groups or ANOVA + Tukey for 3+ groups
# group_names and n_groups come from the user config section 
group <- factor(meta$Group)

if(n_groups == 2){
  # For 2 groups we skip ANOVA and Tukey entirely
  # run_pairwise already runs a t-test internally
  res_pairwise_list <- list(
    
    run_pairwise(group_names[1], group_names[2], x_norm, meta)
  )
  
  # Name the comparison so we can reference it later 
  names(res_pairwise_list) <- paste(group_names[1], "vs", group_names[2])
  
} else{
  
  # ----- One Way Anova Per Metabolite ----- #
  # apply() runs the function on every row (metabolite) of x_norm 
  # MARGIN = 1 means rows, 2 would mean columns 
  anova_res <- apply(x_norm, 1, function(z){
    # Build a small dataframe for each metabolite
    df <- data.frame(
      value = z,
      group = group
    )
    
    # Remove any missing values 
    df <- df[!is.na(df$value), ]
    
    # Safety check - need at least 3 rows and 2 groups to run ANOVA
    if(nrow(df) < 3 || length(unique(df$group)) < 2 || sum(!is.na(df$value)) < 3) {
      return(data.frame(pvalue = NA))
    }
    
    # Fit the ANOVA model
    # aov() is R's built in ANOVA function
    # value ~ group means "model value as a function of group"
    fit <- aov(value ~ group, data = df)
    
    # Extract the p value from the ANOVA summary 
    p <- summary(fit)[[1]][["Pr(>F)"]][[1]]
    
    data.frame(pvalue = p)
    
    
  })
  
  # Combine all the results into one dataframe
  # do.call(rbind) stacks all the individual results row by row
  anova_res <- do.call(rbind, anova_res)
  
  anova_res$feature_id <- rownames(x_norm)
  
  rownames(anova_res) <- NULL
  
  
  # ----- FDR Correction (Benjamini Hochberg) ----- #
  # p.adjust corrects for multiple testing
  # BH method controls the false discovery rate
  # Without this, running thousands of tests inflates the false positive rate 
  anova_res$padj <- p.adjust(anova_res$pvalue, method = "BH")
  
  # Filter for significant features only, remove NAs
  anova_sig_features <- anova_res$feature_id[anova_res$padj < fdr_threshold]
  
  anova_significant_features_clean <- anova_sig_features[!is.na(anova_sig_features)]
  
  #----- Tukey HSD Post Hoc ----- #
  # ANOVA tells us SOMETHING is different between groups
  # Tukey tells us specificanlly WHICH groups differ from each other 
  tukey_list <- lapply(anova_significant_features_clean, function(fid){
    z <- x_norm[fid,]
    
    df <- data.frame(value = z,
                     group = group)
    df <- df[!is.na(df$value), ]
    
    fit <- aov(value ~ group, data = df)
    
    # TukeyHSD performs all pairwise comparisons with correction
    tk <- TukeyHSD(fit)
    
    out <- as.data.frame(tk$group)
    
    out$comparison <- rownames(out)
    
    out$feature_id <- fid
    
    rownames(out) <- NULL
    out
  })
  
  tukey_res <- do.call(rbind, tukey_list)
  
  tukey_sig <- tukey_res %>% filter(`p adj` < fdr_threshold)
  
  # ----- Generate all the pairwise combinations automatically ----- 
  # combn(group_names, 2) generates every possible pair of groups
  # For 3 groups WT, WRN, shPGC it produces:
  # WT-WRN, WT-shPGC, WRN-shPGC
  # Do not need to hardcode each comparison manually
  # Simplify = FALSE returns as a list instead of a matrix
  
  pairs <- combn(group_names, 2, simplify = FALSE)
  
  # lapply loops through each pair and runs run_pairwise
  # The result is a named list of dataframes, one per comparison 
  res_pairwise_list <- lapply(pairs, function(pair){
    run_pairwise(pair[1], pair[2], x_norm, meta)
  })
  
  # Name each result by its comparison for easy reference later
  names(res_pairwise_list) <- sapply(pairs, function(pair){
    paste(pair[1], "vs", pair[2])
  })
  
  
}

# ----- Print summary of statistical testing completed ----- #
cat("\n======================================\n")
cat("Statitsical Testing Complete\n")
cat("========================================\n")
cat("Groups tested:", paste(group_names, collapse = ", "), "\n")
cat("Number of groups:", n_groups, "\n")
cat("Statistical test used:", ifelse(n_groups == 2,"T-test (pairwise)", "ANOVA + Tukey + Pairwise t tests"), "\n")


if(n_groups > 2){
  cat("\nANOVA Results:\n")
  cat(" Total metabolites tested:", nrow(anova_res), "\n")
  cat(" Significant (FDR <", fdr_threshold, "):", length(anova_significant_features_clean), "\n")
  cat(" Tukey Significant comparisons:", nrow(tukey_sig), "\n")
}

cat("\nPairwise Comparisons Run:\n")
for(comparison_name in names(res_pairwise_list)){
  res <- res_pairwise_list[[comparison_name]]
  cat(" -", comparison_name, ":", sum(res$padj < fdr_threshold, na.rm = TRUE), "significant metabolites\n")
}
cat("========================================\n")



cat("=========================================\n")
cat("Output Datasets Created\n")
cat("=========================================\n")

if(n_groups >2){
  cat("\nANOVA & Tuket Datasets:\n")
  cat(" - anova_res:", nrow(anova_res), "metabolites\n")
  cat(" - anova_significant_features_clean:", length(anova_significant_features_clean), "metabolites\n")
  cat(" - tukey_res:", nrow(tukey_res), "rows\n")
  cat(" -tukey_sig:", nrow(tukey_sig), "significant comparisons\n")
  
} else{
  cat("\nT-test Dataset:\n")
  cat(" - res_pairwise_list", nrow(res_pairwise_list[[1]]), "metabolites tested\n")
  cat(" -Significant (FDR <", fdr_threshold, "):", sum(res_pairwise_list[[1]]$padj < fdr_threshold, na.rm = TRUE), "\n")
  
}

cat("\nPairwise Datasets:\n")
for(comparison_name in names(res_pairwise_list)){
  res <- res_pairwise_list[[comparison_name]]
  cat(" - res_pairwise_list[['", comparison_name, "']]: ", nrow(res), " metabolites\n", sep="")
}

cat("\nCore Datasets\n")
cat(" - x_norm:", nrow(x_norm), "metabolites x", ncol(x_norm), "samples\n")
cat(" - meta:", nrow(meta), "samples\n")
cat("====================================\n")

# ------ Visualizations for ANOVA Results ------ #


# Annotation 
annotation_col <- meta %>% select(Group)

rownames(annotation_col) <- meta$Sample

if(n_groups >2){

cat("\n === Anova Summary === \n")
cat("Total significant metabolites tested: ", nrow(anova_res), "\n")
cat("Significant Metabolites (FDR < 0.05): ", length(anova_significant_features_clean), "\n")
cat("Percentage Significant: ", 
    round(length(anova_significant_features_clean)/ nrow(anova_res) * 100, 1), "%\n")

# Heatmap
top_anova_metabolites <-  anova_res %>%
  
  filter(!is.na(feature_id) & padj < fdr_threshold) %>%
  
  arrange(padj) %>%
  
  slice_head(n = 30) %>% # Top 30 for heatmap 
  
  pull(feature_id)
  
pdf(output_path("heatmap_anova_top30.pdf"), width = 14, height =10)

print(pheatmap(
  x_norm[top_anova_metabolites,],
  scale = "row",
  annotation_col = annotation_col,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "average",
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 8,
  
  main = "Top 30 Anova Significant metabolites"

  ))

dev.off()


# Boxplots
top6_anova <- anova_res %>%
  filter(!is.na(feature_id) & padj < fdr_threshold) %>%
  arrange(padj) %>%
  slice_head(n=6) %>%
  pull(feature_id)

plot_data_anova <- x_norm[top6_anova, ] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("metabolite") %>%
  pivot_longer(-metabolite, names_to = "sample", values_to = "abundance") %>%
  left_join(meta, by = c("sample" = "Sample"))


# Create boxplots
pdf(output_path("anova_boxplots_top6.pdf"), width = 8, height = 6)

print(ggplot(plot_data_anova, aes(x = Group, y = abundance, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.5, size = 2) +
  facet_wrap(~ metabolite, scales = "free_y", ncol = 3)+ 
  labs(
    title = "Top 6 Metabolites by ANOVA",
    y = "Log2 Normalized Abundance",
    x = ""
  ) + theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5)
  ))
dev.off()

}

# ----- Add significance labels to each pairwise comparison ----- #
# Loops through every comparison in the res_pairwise_list  
# and adds Up/Down/Not significant based on the fdr_threshold and logfc threshold 
# Both thresholds come from user config section

res_pairwise_list <- lapply(res_pairwise_list, function(res){
  res$significance <- "Not Significant"
  res$significance[res$padj < fdr_threshold & res$logFC > logfc_threshold] <- "Up"
  res$significance[res$padj < fdr_threshold & res$logFC < -logfc_threshold] <- "Down"
  res
})



# Volcano Plot function 
plot_volcano <- function(res, title, top_n = 35) {
  
  top_up <- res %>%
    filter(significance == "Up") %>%
    group_by(significance) %>%
    arrange(padj) %>%
    slice_head(n=top_n)
  
  top_down <- res %>%
    filter(significance == "Down") %>%
    arrange(padj) %>%
    slice_head(n=top_n)
  
  top_labels <- bind_rows(top_up, top_down)
  
  # Count signficant metabolites for subtitles
  n_up <- sum(res$significance == "Up")
  n_down <- sum(res$significance == "Down")
  n_total <- n_up + n_down

    
  # Create the plot 
  ggplot(res, aes(x = logFC, y = -log10(padj))) +
    
    # Background points (Not significant) Plotted first so theyre behind 
    geom_point(data = res %>% filter(significance == "Not Significant"),
               color = "grey70",
               alpha = 0.3,
               size = 1.8) +
    geom_point( data = res %>% filter(significance != "Not Significant"),
                aes(color = significance),
                alpha = 0.8,
                size = 2.5) +
    
    geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed",
               color = "grey50", linewidth = 0.4) +
    
    geom_vline(xintercept =c(-1,1) , linetype = "dashed",
               color = "grey50",
               linewidth = 0.4) +
    
    # Labels with better spacing
    geom_text_repel(
      data = top_labels,
      aes(label = feature_id, color = significance),
      size = 2.5,
      max.overlaps = 40, box.padding = 0.4, point.padding = 0.2, segment.color = "grey40",
      segment.size = 0.3, min.segment.length = 0, force = 3, show.legend = FALSE
    ) +
    
    
    scale_color_manual(values = c(
      "Up" = "#D62728",
      "Down" = "#1F77B4",
      "Not Significant" = "grey70"
    ),
    
    labels = c(
      "Up" = "Upregulated",
      "Down" = "Downregulated",
      "Not Significant" = "Not Significant"
    ), 
    name = ""
    
    ) +
    
    labs(
      title = title , 
      subtitle = paste0(n_total, " Significant metabolites (", n_up, " up, ", n_down, " down)"),
      x = "Log2 Fold Change",
      y = "-log10(FDR)"
    ) +
    
    theme_minimal() +
    theme(
      # Title styling
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey30"),
      
      # Axis Styling
      axis.title = element_text(size = 13, face = "bold"),
      axis.text = element_text(size = 11, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      
      # Legend styling
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      legend.title = element_blank(),
      
      # Grid and panel
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color="grey90", linewidth = 0.3),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      
      # Background
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill="white", color = NA)
    )
    
}

#----- Generate volcano plots for each pairwise comparison ----- #
# Loops through res_pairwise_list and creates on volcano plot per comparison
# File is named automatically based on comparison name

for(comparison_name in names(res_pairwise_list)){
  
  # Get the result dataframe for this comparison
  res <- res_pairwise_list[[comparison_name]]
  
  # Create a clean filename by replacing spaces with underscores
  # Example: "NL vs Y" becomes "NL_vs_Y"
  filename <- paste0("volcano_plot_", gsub(" ", "_", comparison_name), ".pdf")
  
  pdf(output_path(filename), width = 8, height = 6)
  
  print(plot_volcano(res, comparison_name))
  
  dev.off()
}


#------ Fold change Bar Plots (Pairwise) ----- #
# Reusable function to create FC Bar Plots
plot_fc_barplot <- function(res_df, comparison_name, n_metabolites = 10){
  
  # Get the top metabolites by absolute fold change 
  top_fc <- res_df %>%
    filter(padj < fdr_threshold) %>%
    arrange(desc(abs(logFC))) %>%
    slice_head(n = n_metabolites)
  
  # Determine groups for labeling
  groups <- strsplit(comparison_name, " vs ")[[1]]
  group1 <- groups[1]
  group2 <- groups[2]
  
  # Create the plot
  ggplot(top_fc, aes(x = reorder(clean_name, logFC), y = logFC, fill = logFC > 0)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    coord_flip() +
    scale_fill_manual(
      values = c("TRUE" = "#E64B35", "FALSE" = "#4DbbD5"),
      labels = c("TRUE" = paste("Higher in", group1), "FALSE" = paste("Higher in", group2)),
      name = ""
    ) +
    labs(
      title = paste("Top", n_metabolites , "Metabolites by Fold Change"),
      subtitle = comparison_name,
      x = "",
      y = "Log2 Fold Change"
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text (hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(size = 6),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}


#----- Clean each pairwise result ----- #
# Remove the NA entries, unnamed features (X1, X2, etc.)
# and duplicate feature IDs from every comparison
# Results stored in a new named list called res_pairwise_clean
res_pairwise_clean <- lapply(res_pairwise_list, function(res){
  res %>%
    filter(!is.na(feature_id)) %>%
    filter(!grepl("^NA[0-9\\.]", feature_id)) %>%
    filter(!grepl("^X[0-9]^X\\.", feature_id)) %>%
    distinct(feature_id, .keep_all = TRUE) %>%
    rename(clean_name = feature_id)
})

# ----- Generate barplots for each pairwise comparison ----- #
# Loops through res_pairwise_clean and creates one bar plot per comparison
# Uses the cleaned results so unnamed/NA metabolites are excluded
# File name created automatically based on comparison name

for(comparison_name in names(res_pairwise_clean)){
  # Get the cleaned result dataframe for this comparison
  res_clean <- res_pairwise_clean[[comparison_name]]
  
  # Create the filename automatically
  # Example: "NL vs Y" becomes "NL_vs_Y.pdf"
  filename <- paste0("barplot_", gsub(" ", "_", comparison_name), ".pdf")
  
  pdf(output_path(filename), width = 8, height = 6)
  print(plot_fc_barplot(res_clean, comparison_name, n_metabolites = 15))
  dev.off()
  }


# ----- Generate heat maps for each pairwise comparison ----- #
# Loops through res_pairwise_clean and creates one heatmap per comparison
# Top 50 significant metabolites per comparison
# File name automatically created based on comparison name 
 
for(comparison_name in names(res_pairwise_clean)){
  
  # Get cleaned results for this comparison
  res_clean <- res_pairwise_clean[[comparison_name]]
  
  # Get top 50 significant metabolites by adjusted p-value
  top_metabs <- res_clean %>%
    filter(padj < fdr_threshold) %>%
    arrange(padj) %>%
    slice_head(n = 50) %>%
    pull(clean_name)
  
  # Skip this comparison if there are no significant metabolites found
  if(length(top_metabs) == 0){
    cat("No significant metabolites for", comparison_name, "-skipping heatmpa\n")
    next
  }
  
  # Extract the normalized data for the metabolites
  heat_data <- x_norm[top_metabs, ]
  
  # Create the file name automatically
  filename <- paste0("heatmap_", gsub(" ", "_", comparison_name), ".pdf")
  
  pdf(output_path(filename), width = 12, height = 12)
  
  print(pheatmap(
    heat_data,
    scale = "row",
    annotation_col = annotation_col,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "average",
    fontsize_row = 7,
    fontsize_col = 7,
    treeheight_row = 120,
    treeheight_col = 120,
    main = paste("Top", nrow(heat_data), "Significant Metabolites (", comparison_name, ")")
  ))
  
  dev.off()
  
}


#----- Export pairwise results to CSV ----- #
# Loops through both lists and exports full and cleaned results
# Files named automatically based on comparison name 

for(comparison_name in names(res_pairwise_list)){
  # Create clean filename
  clean_comparison <- gsub(" ", "_", comparison_name)
  
  # Export full results
  write.csv(
    res_pairwise_list[[comparison_name]],
    output_path(paste0("pairwise_", clean_comparison, "_full.csv")),
    row.names = FALSE
  )
  
  # Export cleaned results
  write.csv(
    res_pairwise_clean[[comparison_name]],
    output_path(paste0("pairwise_", clean_comparison, "_clean.csv")),
    row.names = FALSE
  )
  
  # Export significant metabolites for Metaboanalyst Pathway Analysis
  sig_metabs <- res_pairwise_clean[[comparison_name]] %>%
    filter(padj < fdr_threshold) %>%
    pull(clean_name)
  
  write.csv(
    data.frame(metabolite = sig_metabs),
    output_path(paste0("pathway_input_", clean_comparison, ".csv")),
    row.names = FALSE
  )
  
}

# Export ANOVA and Tukey Results only if 3+ groups
if(n_groups > 2){
  write.csv(anova_res, output_path("anova_all_results.csv"), row.names = FALSE)
  write.csv(anova_res %>% filter(padj < fdr_threshold), output_path("anova_significant_results.csv"), row.names = FALSE)
  write.csv(tukey_res, output_path("tukey_all_results.csv"), row.names = FALSE)
  write.csv(tukey_sig, output_path("tukey_significant_results.csv"), row.names = FALSE)
}

# Export normalized data matrix and metadata
x_norm_export <- x_norm %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Metabolite")
write.csv(x_norm_export, output_path("normalized_data_matrix.csv"), row.names = FALSE)
write.csv(meta, output_path("sample_metadata.csv"), row.names = FALSE)


#----- Dynamic summary statistics ----- #
# Builds a summary for every pairwise comparison automatically
# Also includes ANOVA summary if n_groups > 2

# Pairwise summary for every comparison
pairwise_summary <- do.call(rbind, lapply(names(res_pairwise_list), function(comparison_name){
  res <- res_pairwise_list[[comparison_name]]
  res_clean <- res_pairwise_clean[[comparison_name]]
  
  data.frame(
    Comparison = comparison_name,
    Total_Tested = nrow(res),
    Significant_FDR = sum(res$padj < fdr_threshold, na.rm = TRUE),
    Upregulated = sum(res$significance == "Up", na.rm = TRUE),
    Downregulated = sum(res$significance == "Down", na.rm = TRUE),
    High_FC_Significant = sum(res$significance != "Not Significant", na.rm = TRUE),
    Identifiable_Metabolites = nrow(res_clean)
  )
  
}))

write.csv(pairwise_summary, output_path("summary_statistics_pairwise.csv"), row.names = FALSE)

# ANOVA summary only for 3+ groups 
if(n_groups > 2){
  
  anova_summary <- data.frame(
    Analysis = c(
      "Total Metabolites Tested",
      "Metabolites after filtering (>70% present)",
      "ANOVA Significant (FDR < 0.05)",
      "Tukey Significant Comparisons"
    ),
    
    Count = c(
      nrow(x),
      nrow(x_norm),
      sum(anova_res$padj < fdr_threshold, na.rm = TRUE),
      nrow(tukey_sig)
    )
  )
  
  write.csv(anova_summary, output_path("summary_statistics_anova.csv"), row.names = FALSE)
  
}
























