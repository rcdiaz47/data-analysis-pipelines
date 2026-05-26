# ----- Statistical Analysis ----- 
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

input_file <- "LH_9033517_20260212_RP_Pos_CD_Results.xlsx"

# Number of groups to be tested in experiment
# 2 = t-test only
# 3+ = ANOVAA + Tukey + Pairwise analysis

n_groups <- 2

# Group names must match sample column names
group_names <- c("NL", "Y")

# Output directory for plots and csv files ("./) for current directory
output_dir <- "./"

# FDR significant threshold
fdr_threshold <- 0.05

# Log2fold change threshold for plots
logfc_threshold <- 1


#=========================================
# END USER CONFIGURATION
#=========================================


cd <- read_xlsx(input_file)
head(cd)

## ------ Clean and shorten the metabolite names for readability ------ 
process_metabolite_names <- function(x, max_length = 40){
  x %>%
    gsub("\\.", " ", .) %>% #Replace dots with spaces
    gsub("\\.\\d+$", "", .) %>% #Remove trailing .1, .2, etc. 
    gsub("\\s+", " ", .) %>% # Replace multiple spaces with single space
    trimws() %>% #Remove leading/trailing whitespace
  
  ifelse(nchar(.) > max_length,
         paste0(substr(., 1, max_length - 3), "..."), .)
}

cd$Name <- process_metabolite_names(cd$Name, max_length = 40)


# ----- Build the numeric matrix with the columns we want to use for analyzing ----- 
area_cols <- grep("^Area:", colnames(cd), value = TRUE)

## Exclude blanks and QCs from the analysis
area_cols <- area_cols[!grepl("Blank|blank|BLANK|QC|Qc", area_cols)]


x <- cd %>%
  dplyr::select(all_of(area_cols)) %>% as.matrix()

rownames(x) <- cd$Name

# ----- Replace missing values with Nas ----- 
sum(is.na(x))
x[x== ""] <- NA

# ----- Log transform the data with log2 ------ 
x_log <- log2(x)

# ----- Drop metabolites that are missing in too many samples ----- 
keep <- rowMeans(!is.na(x_log)) >= 0.7
x_filt <- x_log[keep,]


# ----- Data normalization using Median Centering ----- 
column_median <- apply(x_filt,2,median, na.rm = TRUE)

# Subtract each columns median from the value 
x_norm <- sweep(
  x_filt, 2, column_median, FUN = "-"
)

summary(x_norm[,1])
any(is.na(x_norm))
nrow(x_norm)

# --- Prepare data for plotting; Box and whisker plots and Density plot --- 

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

# --- Density Plots --- 
p1 <- ggplot(before_df, aes(x=value)) + 
  geom_density(fill = "blue", alpha = 0.4) +
  labs(title ="Before Normalization" , x = "Raw Intensity", y = "Density") +
  theme_minimal()

p2 <- ggplot(after_df, aes(x = values)) +
  geom_density(fill = "blue", alpha = 0.4) +
  labs(title = "After Normalization", x = "Normalized Intensity", y = "Density") +
  theme_minimal()

# --- Boxplots (Subset so they dont look all compressed) ---- 

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

## ---- Combine all the plots into one figure ----
data_boxplots <- (p3|p4)
data_density_plots <- (p1|p2)

pdf("data_boxplots_before_norm_2.pdf", width = 10, height = 6)
plot(data_boxplots)
dev.off()

pdf("data_density_plots_before_norm.pdf", width = 8, height = 6)
plot(data_density_plots)
dev.off()


##------ Exploratory Data Analysis ---------

# ----- build metadata to Begin statistical testing (Anova, t test, PCA) ----- 
meta <- data.frame(
  Sample = area_cols, 
  Group = dplyr::case_when(
    grepl("NL", area_cols) ~ "NL",
    grepl("Y", area_cols) ~ "Y"
  ),
  stringsAsFactors = FALSE
)

meta$Group <- factor(meta$Group, levels = c("NL", "Y"))
meta



## ---- PCA plot for all of the metabolites to visualize the group separation ----

## Transpose the data because PCA needs the samples as the rows instead of columns 

x_pca <- t(na.omit(x_norm))

## Run PCA
pca_res <- prcomp(x_pca, scale. = TRUE)

## Extract the PCA scores and add group information to prepare for plotting 
pca_scores <- as.data.frame(pca_res$x)
pca_scores$Group <- meta$Group

## Plot 
pdf("pca_plot.pdf", width = 8, height = 6)

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



# ----- Run statistical test one way anova per metabolite (Are there significant differences between groups?) ----- 
## Only for more than 2 groups 
group <- factor(meta$Group)

anova_res <- apply(x_norm, 1, function(z){
  
  df <- data.frame(
    value = z,
    group = group
  )
  
  df <- df[!is.na(df$value), ]
  
  if(nrow(df) < 3 || length(unique(df$group)) < 2 || sum(!is.na(df$value)) < 3){
    return(data.frame(pvalue = NA))
  }
  
  fit <- aov(value ~ group, data = df)
  
  p <- summary(fit)[[1]][["Pr(>F)"]][1]
  
  data.frame(pvalue = p)
})

anova_res <- do.call(rbind, anova_res)
anova_res$feature_id <- rownames(x_norm)
rownames(anova_res) <- NULL



# ----- FDR correction (Benjamini Hochberg) ----- 

anova_res$padj <- p.adjust(anova_res$pvalue, method = "BH")
head(anova_res)
nrow(anova_res)


# --- Filter the anova results for significance and also remove NAs for tukey analyis 
anova_sig_features <- anova_res$feature_id[anova_res$padj < 0.05]
anova_significant_features_clean <- anova_sig_features[!is.na(anova_sig_features)]
length(anova_significant_features_clean)



# ----- Tukey HSD Post Hoc (Where are the significant differences?) ----- 
tukey_list <- lapply(anova_significant_features_clean, function(fid){
  z <- x_norm[fid,]
  
  df <- data.frame(
    value = z,
    group = group
  )
  df <- df[!is.na(df$value), ]
  
  fit <- aov(value ~ group, data = df)
  
  tk <- TukeyHSD(fit)
  
  out <- as.data.frame(tk$group)
  out$comparison <- rownames(out)
  out$feature_id <- fid
  
  rownames(out) <- NULL
  out
  
})

tukey_res <- do.call(rbind, tukey_list)
tukey_sig <- tukey_res %>% filter(`p adj` < 0.05)
nrow(tukey_sig)
head(tukey_sig)

## ------ Visualizations for ANOVA Results ------ ##
## ----------------------------------------------##
cat("\n === Anova Summary === \n")
cat("Total significant metabolites tested: ", nrow(anova_res), "\n")
cat("Significant Metabolites (FDR < 0.05): ", length(anova_significant_features_clean), "\n")
cat("Percentage Significant: ", 
    round(length(anova_significant_features_clean)/ nrow(anova_res) * 100, 1), "%\n")

# Get the top 30 metabolites by ANOVA for plotting 
top_anova_metabolites <-  anova_res %>%
  filter(!is.na(feature_id) & padj < 0.05) %>%
  arrange(padj) %>%
  slice_head(n = 30) %>% # Top 30 for heatmap 
  pull(feature_id)
  

#Annotation 
annotation_col <- meta %>% column_to_rownames("Sample") %>%
  select(Group)

## ----- Heatmap: Anova Significant metabolites (All 3 Groups) ----- 

pdf("heatmap_anova_top30.pdf", width = 14, height =10)

pheatmap(
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

  )

dev.off()

### ----- PCA: ANOVA Significant Metabolites only ----- 
## Probably a little redundant
# Use only ANOVA significant metabolites for PCA
x_pca_anova <- t(na.omit(x_norm[anova_significant_features_clean,]))

# Run PCA 
pca_anova <- prcomp(x_pca_anova, scale. = TRUE)


# Extract the scores
pca_anova_scores <- as.data.frame(pca_anova$x)
pca_anova_scores$Group <- meta$Group

# Plot

pdf("pca_anova_significant.pdf", width = 8, height = 6)

ggplot(pca_anova_scores, aes(x =PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(
    aes(fill = Group),
    geom = "polygon",
    alpha = 0.2,
    level = 0.95,
    color = NA
  ) +
  labs(
    title = " PCA: ANOVA- Significant Metabolites Only",
    subtitle = paste0(length(anova_significant_features_clean), " metabolites"),
    x = paste0("PC1 (", round(summary(pca_anova)$importance[2,1] * 100,1), "%)"),
    y = paste0("PC2 (", round(summary(pca_anova)$importance[2,2] * 100,1), "%)")
  ) + theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

dev.off()

## ------ Box plots: Top 6 Anova Metabolites ----- 
top6_anova <- anova_res %>%
  filter(!is.na(feature_id) & padj < 0.05) %>%
  arrange(padj) %>%
  slice_head(n=6) %>%
  pull(feature_id)

plot_data_anova <- x_norm[top6_anova, ] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("metabolite") %>%
  pivot_longer(-metabolite, names_to = "sample", values_to = "abundance") %>%
  left_join(meta, by = c("sample" = "Sample"))


# Create boxplots
pdf("anova_boxplots_top6.pdf", width = 8, height = 6)

ggplot(plot_data_anova, aes(x = Group, y = abundance, fill = Group)) +
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
  )
dev.off()

## ---- Reusable function for pairwise analysis of groups ---- 

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
    
    #LogFC (difference in means)
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

## Results for all the pairwise comparisons and added significance columns for the volcano plot

## WT vs WRN
res_NL_Y <- run_pairwise("NL", "Y", x_norm, meta)
res_NL_Y$significance <- "Not Significant"
res_NL_Y$significance[res_NL_Y$padj < 0.05 & res_NL_Y$logFC > 1] <- "Up"
res_NL_Y$significance[res_NL_Y$padj < 0.05 & res_NL_Y$logFC < -1] <- "Down"

cat("\n=== NL vs Y Summary ===\n")
cat(" Significant Metabolites:", sum(res_NL_Y$padj < 0.05, na.rm = TRUE), "\n")
cat("Upregulated in WT:", sum(res_NL_Y$significance == "Up"), "\n")
cat("Downregulated in WT:", sum(res_NL_Y$significance == "Down"), "\n")


## WT vs WRN
res_WT_WRN <- run_pairwise("WT", "WRN", x_norm, meta)
res_WT_WRN$significance <- "Not Significant"
res_WT_WRN$significance[res_WT_WRN$padj < 0.05 & res_WT_WRN$logFC > 1] <- "Up"
res_WT_WRN$significance[res_WT_WRN$padj < 0.05 & res_WT_WRN$logFC < -1] <- "Down"

cat("\n=== WT vs WRN Summary ===\n")
cat(" Significant Metabolites:", sum(res_WT_WRN$padj < 0.05, na.rm = TRUE), "\n")
cat("Upregulated in WT:", sum(res_WT_WRN$significance == "Up"), "\n")
cat("Downregulated in WT:", sum(res_WT_WRN$significance == "Down"), "\n")

# WT vs shPGC
res_WT_shPGC <- run_pairwise("WT", "shPGC", x_norm, meta)
res_WT_shPGC$significance <- "Not Significant"
res_WT_shPGC$significance[res_WT_shPGC$padj < 0.05 & res_WT_shPGC$logFC > 1] <- "Up"
res_WT_shPGC$significance[res_WT_shPGC$padj < 0.05 & res_WT_shPGC$logFC < -1] <- "Down"

cat("\n=== WT vs shPGC Summary ===\n")
cat(" Significant Metabolites:", sum(res_WT_shPGC$padj < 0.05, na.rm = TRUE), "\n")
cat("Upregulated in WT:", sum(res_WT_shPGC$significance == "Up"), "\n")
cat("Downregulated in WT:", sum(res_WT_shPGC$significance == "Down"), "\n")

# WRN vs shPGC
res_WRN_shPGC <- run_pairwise("WRN", "shPGC", x_norm, meta)
res_WRN_shPGC$significance <- "Not Significant"
res_WRN_shPGC$significance[res_WRN_shPGC$padj < 0.05 & res_WRN_shPGC$logFC > 1] <- "Up"
res_WRN_shPGC$significance[res_WRN_shPGC$padj < 0.05 & res_WRN_shPGC$logFC < -1] <- "Down"

cat("\n=== WRN vs shPGC Summary ===\n")
cat(" Significant Metabolites:", sum(res_WRN_shPGC$padj < 0.05, na.rm = TRUE), "\n")
cat("Upregulated in WRN:", sum(res_WRN_shPGC$significance == "Up"), "\n")
cat("Downregulated in WRN:", sum(res_WRN_shPGC$significance == "Down"), "\n")

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
    
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
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
      
      #Legend styling
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      legend.title = element_blank(),
      
      #Grid and panel
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color="grey90", linewidth = 0.3),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      
      #Background
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill="white", color = NA)
    )
    
}

## Create pdfs of volcano plots for each comparison

pdf("volcano_plot_NL_VS_Y.pdf", width = 8, height = 6)
plot_volcano(res_NL_Y, "NL vs Y", top_n = 35)
dev.off()



pdf("volcano_plot_WT_VS_WRN.pdf", width = 8, height = 6)
plot_volcano(res_WT_WRN, "WT vs WRN", top_n = 10)
dev.off()

pdf("volcano_plot_WT_vs_shPGC.pdf", width = 8, height = 6)
plot_volcano(res_WT_shPGC, "WT vs shPGC")
dev.off()

pdf("volcano_plot_WRN_shPGC.pdf", width = 8, height = 6)
plot_volcano(res_WRN_shPGC, "WRN vs shPGC")
dev.off()

##------ Fold change Bar Plots (Pairwise) ----- ##
## Reusable function to create FC Bar Plots
plot_fc_barplot <- function(res_df, comparison_name, n_metabolites = 10){
  
  # Get the top metabolites by absolute fold change 
  top_fc <- res_df %>%
    filter(padj < 0.05) %>%
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


## Now we will remove unusable entries like NAs and also X1...X2...for each comparison
## This is removing entries, making the data easier to work with
res_NL_Y_clean <- res_NL_Y %>%
  filter(!is.na(feature_id)) %>%
  filter(!grepl("^NA[0-9\\.]", feature_id)) %>%
  filter(!grepl("^X[0-9]|^X\\.", feature_id)) %>%
  distinct(feature_id, .keep_all = TRUE) %>%
  rename(clean_name = feature_id)



res_WT_WRN_clean <- res_WT_WRN %>%
  filter(!is.na(feature_id)) %>%
  filter(!grepl("^NA", feature_id)) %>%
  filter(!grepl("^X[0-9]|^X\\.", feature_id)) %>%
  distinct(feature_id, .keep_all = TRUE) %>%
  rename(clean_name = feature_id)

res_WT_shPGC_clean <- res_WT_shPGC %>%
  filter(!is.na(feature_id)) %>%
  filter(!grepl("^NA", feature_id)) %>%
  filter(!grepl("^X[0-9]|^X\\.", feature_id)) %>%
  distinct(feature_id, .keep_all = TRUE) %>%
  rename(clean_name = feature_id)

res_WRN_shPGC_clean <- res_WRN_shPGC %>%
  filter(!is.na(feature_id)) %>%
  filter(!grepl("^NA", feature_id)) %>%
  filter(!grepl("^X[0-9]|^X\\.", feature_id)) %>%
  distinct(feature_id, .keep_all = TRUE) %>%
  rename(clean_name = feature_id)


# Create bar plots for each comparisons top ten significant metabolites 

pdf("Log2Fold_Bar_chart_NL_vs_Y.pdf", width = 8, height = 6)
p1 <- plot_fc_barplot(res_NL_Y_clean, "NL vs Y", n_metabolites = 50)
p1
dev.off()


pdf("Bar_chart_WT_vs_WRN.pdf", width = 8, height = 6)
p1 <- plot_fc_barplot(res_WT_WRN_clean, "WT vs WRN", n_metabolites = 10)
p1
dev.off()

pdf("Bar_chart_WT_vs_shPGC.pdf", width = 8, height = 6)
p2 <- plot_fc_barplot(res_WT_shPGC_clean, "WT vs shPGC", n_metabolites = 10)
p2
dev.off()

pdf("Bar_chart_WRN_vs_shPGC.pdf", width = 8, height = 6)
p3 <- plot_fc_barplot(res_WRN_shPGC_clean, "WRN vs shPGC", n_metabolites = 10)
p3
dev.off()


## Now we will extract the significant metabolites from the cleaned results dataframe for each comparison

sig_metab_NL_Y <- res_NL_Y_clean %>%
  filter(padj < 0.05) %>%
  pull(clean_name)


sig_metab_WT_WRN <- res_WT_WRN_clean %>%
  filter(padj < 0.05) %>%
  pull(clean_name)

sig_metab_WT_shPGC <- res_WT_shPGC_clean %>%
  filter(padj < 0.05) %>%
  pull(clean_name)

sig_metab_WRN_shPGC <- res_WRN_shPGC_clean %>%
  filter(padj < 0.05) %>%
  pull(clean_name)


## prepare the input for pathway analysis for each comparison
## we need one column in plain text 

NL_Y_sig_df <- data.frame(metabolite = sig_metab_NL_Y)


Wt_wrn_sig_df <- data.frame(metabolite = sig_metab_WT_WRN)
Wt_shPGC_sig_df <- data.frame(metabolite = sig_metab_WT_shPGC)
Wrn_shPGC_sig_df <- data.frame(metabolite = sig_metab_WRN_shPGC)


## Export for metaboanlyst evaluation for pathway analysis
write.csv(NL_Y_sig_df, "NL_vs_y_sig_metabolites_pathway_enrichment.csv", row.names = FALSE)

write.csv(Wt_wrn_sig_df, "wt_vs_wrn_sig_metabolites.csv", row.names = FALSE)
write.csv(Wt_shPGC_sig_df, "wt_vs_shPGC_sig_metabolites.csv", row.names = FALSE)
write.csv(Wrn_shPGC_sig_df, "wrn_vs_shPGC_sig_metabolites.csv", row.names = FALSE)


## Create an abundance heatmap using our normalized matrix and each comparison 
## Top 20 metabolites per comparison 

top_metabs_heatmap_NLvsY<- res_NL_Y_clean %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  slice_head(n=50) %>%
  pull(clean_name)



## WT vs WRN

top_metabs_heatmap_WTvsWRN<- res_WT_WRN_clean %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  slice_head(n=20) %>%
  pull(clean_name)

top_metabs_heatmap_WTvsshPGC<- res_WT_shPGC_clean %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  slice_head(n=20) %>%
  pull(clean_name)

top_metabs_heatmap_WRNvsshPGC<- res_WRN_shPGC_clean %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  slice_head(n=20) %>%
  pull(clean_name)


# Create the heatmap data for each of the 3 comparisons 
NL_vs_Y_heatData <- x_norm[top_metabs_heatmap_NLvsY,]


WT_vs_WRN_heatData <- x_norm[top_metabs_heatmap_WTvsWRN,]
WT_vs_shPGC_heatData <- x_norm[top_metabs_heatmap_WTvsshPGC, ]
WRN_vs_shPGC_heatData <- x_norm[top_metabs_heatmap_WRNvsshPGC, ]


#Annotation 
annotation_col <- meta %>% column_to_rownames("Sample") %>%
  select(Group)


# Abdundance Heatmap for each of the three comparisons
# Include Dendrograms for hierarchal clustering of samples and metabolites 

##PHeatmap is using "euclidean" distance metric with "complete" linkage method by default
## We will switch our method to "correlation" distance so we can see how metabolites co vary (go up or down together)
## We will also siwtch to "average" linkage 
## scaling is z score normalization by row (feature)

pdf("heatmap_NL_vs_Y_pairwise.pdf", width = 14, height = 10)

pheatmap(
  NL_vs_Y_heatData,
  scale = "row",
  annotation_col = annotation_col,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "average",
  fontsize_row = 7,
  fontsize_col = 7,
  treeheight_row = 100,
  treeheight_col = 100,
  main = paste("Top", nrow(NL_vs_Y_heatData), "Significant Metabolites (NL vs Y)")
)
dev.off()

pdf("heatmap_Wt_vs_shPGC_pairwise.pdf", width = 12, height = 12)
pheatmap(
  WT_vs_shPGC_heatData,
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
  main = paste("Top", nrow(WT_vs_shPGC_heatData) , "Significant Metabolites (WT vs shPGC)")
)
dev.off()

pdf("heatmap_Wrn_vs_shPGC_pairwise.pdf", width = 12, height = 12)
pheatmap(
  WRN_vs_shPGC_heatData,
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
  main = paste("Top", nrow(WRN_vs_shPGC_heatData) ,"Significant Metabolites (WRN vs shPGC)")
)
dev.off()

## =====================================
## Export Results to CSV files 
## ====================================

## Anova all and significant results

write.csv(anova_res, "anova_all_results.csv", row.names = FALSE)
write.csv(anova_res %>% filter(padj < 0.05), "anova_significant_results.csv", row.names = FALSE)

## Tukey all and significant results
write.csv(tukey_res, "tukey_all_results.csv", row.names = FALSE)
write.csv(tukey_sig, "tukey_significant_results.csv", row.names = FALSE)

## Pairwise T test Results full and cleaned
write.csv(res_NL_Y, "pairwise_NL_vs_Y_full.csv", row.names = FALSE)
write.csv(res_NL_Y_clean, "pairwise_NL_vs_Y_cleaned.csv", row.names = FALSE)


write.csv(res_WT_WRN_clean, "pairwise_WT_vs_WRN.csv", row.names = FALSE)
write.csv(res_WT_shPGC_clean, "pairwise_WT_vs_shPGC.csv", row.names = FALSE)
write.csv(res_WRN_shPGC_clean, "pairwise_WRN_vs_shPGC.csv", row.names = FALSE)


## Normalized Data Matrix
x_norm_export <- x_norm %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Metabolite")
write.csv(x_norm_export, "normalized_data_matrix.csv", row.names = FALSE)

## Metadata
write.csv(meta, "sample_metadata.csv", row.names = FALSE)

## Summary Statistics for ANOVA

#summary_stats <- data.frame(
  #Analysis = c(#"Total ANOVA Metabolites Tested",
               #"ANOVA Significant (FDR < 0.05)",
               #"Tukey Significant comparisons",
               "NL vs Y significant"
               #"WT vs shPGC significant",
               #"WRN vs shPGC significant"
              # ),
  #Count = c(
    #nrow(anova_res),
    #sum(anova_res$padj < 0.05, na.rm = TRUE),
    #nrow(tukey_sig),
    #sum(res_NL_Y_clean$padj < 0.05, na.rm = TRUE)
    #sum(res_WT_shPGC_clean$padj < 0.05, na.rm = TRUE),
    #sum(res_WRN_shPGC_clean$padj < 0.05, na.rm = TRUE)
  #)
#)


## Summary statistics for Pairwise analysis
summary_stats <- data.frame(
  Analysis = c("Total Metabolites Tested",
               "Metabolites after filtering (>70% present)",
               "Signficant Metabolites (FDR < 0.05)",
               " -Upregulated in NL",
               " -Downregulated in NL",
               "Metabolites with |LogFC| > 1 And FDR < 0.05",
               " -High FC Upregulated (LogFC > 1)",
               " -High FC Downregulated (LogFC < -1)",
               "Metabolites with |LogFC| > 2",
               "Identifiable metabolites (after cleaning)"
               ),
  Count = c(
    nrow(x), # Total before filtering
    nrow(x_norm), # After filtering missing data
    sum(res_NL_Y$padj < 0.05, na.rm = TRUE), # All significant
    sum(res_NL_Y$padj < 0.05 & res_NL_Y$logFC > 0, na.rm = TRUE), #Up
    sum(res_NL_Y$padj < 0.05 & res_NL_Y$logFC < 0, na.rm = TRUE), #Down
    sum(res_NL_Y$significance != "Not Significant"), #High FC + significant
    sum(res_NL_Y$significance == "Up"), # High FC up
    sum(res_NL_Y$significance == "Down"), # High FC Down
    sum(res_NL_Y$padj < 0.05 & abs(res_NL_Y$logFC > 2), na.rm = TRUE), #Very high FC
    nrow(res_NL_Y_clean)
  )
)

write.csv(summary_stats, "summary_statistics.csv", row.names = FALSE)

























