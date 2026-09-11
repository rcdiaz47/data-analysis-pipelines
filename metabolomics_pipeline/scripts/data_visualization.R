# ------ Visualizations ------ #
library(tidyverse)
library(yaml)
library(pheatmap)
library(ggplot2)
library(ggrepel)

# read the pipeline configurations
config <- yaml::read_yaml("C:\\Users\Proteomics\\Documents\\data-analysis-pipelines\\metabolomics_pipeline\\config.yaml")

output_dir <- config$output_dir
n_groups <- config$n_groups

intermediate_dir <- file.path(output_dir, "intermediate")

output_path <- function(filename){
  file.path(output_dir, filename)
}

intermediate_path <- function(filename){
  file.path(intermediate_dir, filename)
}

# fixed thresholds
fdr_threshold <- 0.05
logfc_threshold <- 1

# ------ Read inputs from earlier stages ------ #

comparison_names_path <- intermediate_path("comparison_names.txt")

if(!file.exists(comparison_names_path)){
  stop(
    paste0(
      "comparison_names.txt not found:", comparison_names_path,
      "\n Run statistical_analysis.R first."
    )
  )
}

comparison_names <- readLines(comparison_names_path)

x_norm_df <- read.csv(intermediate_path("normalized_matrix.csv"), stringsAsFactors = FALSE)
meta <- read.csv(intermediate_path("sample_metadata.csv"), stringsAsFactors = FALSE)

x_norm <- x_norm_df %>%
  column_to_rownames("feature_id") %>%
  as.matrix()

# Read each comparisons pairwise results into a named list 
res_pairwise_list <- lappy(comparison_names, function(comparison_name){
  
  clean_comparison <- gsub(" ", "_", comparison_name)
  read.csv(
    intermediate_path(paste0("pairwise_", clean_comparison, "_full.csv")),
    stringsAsFactors = FALSE
  )
  
})

names(res_pairwise_list) <- comparison_names

# Annotation for heatmaps
annotation_col <- meta %>% select(Group)
rownames(annotation_col) <- meta$Sample

# ------ Volcano plot function ------#

plot_volcano <- function(res, title, top_n = 35){
  
  top_up <- res %>%
    filter(significance ==  "Up") %>%
    arrange(padj) %>%
    slice_head(n = top_n)
  
  top_down <- res %>%
    filter(significance == "Down") %>%
    arrange(padj) %>%
    slice_head(n = top_n)
  
  top_labels <- bind_rows(top_up, top_down)
  
  n_up <- sum(res$significance == "Up")
  n_down <- sum(res$significance == "Down")
  n_total <- n_up + n_down
  
  ggplot(res, aes(x = logFC, y = -log10(padj))) +
    
    geom_point(
      data = res %>% filter(significance == "Not Significant"),
      color = "grey70", alpha = 0.3, size = 1.8
    ) +
    
    geom_point(
      data = res %>% filter(significance != "Not Significant"),
      aes(color = significance), alpha = 0.8, size = 2.5
    ) +
    
    geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "grey50", linewidth = 0.4) +
    
    geom_text_repel(
      data = top_labels,
      aes(label = feature_id, color = significance),
      size = 2.5, max.overlaps = 40, box.padding = 0.4, point.padding = 0.2,
      segment.color = "grey40", segment.size = 0.3, min.segment.length = 0,
      force = 3, show.legend = FALSE
    ) +
    
    scale_color_manual(
      values = c("Up" = "#D62728", "Down" = "#1F77B4", "Not Significant" = "grey70"),
      labels = c("Up" = "Upregulated", "Down" = "Downregulated", "Not Significant" = "Not Significant"),
      name = ''
    ) +
    labs(
      title = title,
      subtitle = paste0(n_total, "Significant metabolites (", n_up, "up, ",  n_down, "down)"),
      x = "Log2 Fold Change",
      y = "-log10(FDR)"
    ) +
    
    theme_minimal() +
    
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey30"),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text = element_text(size = 11, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

# ------ Fold change bar plots function -----#

plot_fc_barplot <- function(res, comparison_name, n_metanolites = 20){
  
  top_fc <- res %>% 
    filter(padj < fdr_threshold) %>%
    arrange(desc(abs(logFC))) %>% 
    slice_head(n = n_metanolites)
  
  groups <- strsplit(comparison_name, "vs")[[1]]
  group1 <- groups[1]
  group2 <- groups[2]
  
  ggplot(top_fc, aes(x = reorder(feature_id, logFC), y = logFC, fill = logFC > 0)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    coord_flip() +
    scale_fill_manual(
      values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"),
      labels = c("TRUE" = paste("Higher in", group1), "FALSE" = paste("Higher in", group2)),
      name = ""
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(size = 6),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
    
  
}

# ------ Generate volcano plots and fold-change bar plots per comparison ----- # 
for(comparison_name in names(res_pairwise_list)){
  
  res <- res_pairwise_list[[comparison_name]]
  clean_comparison <- gsub(" ", "_", comparison_name)
  
  # Volcano plots 
  pdf(output_path(paste0("volcano_plot_", clean_comparison, ".pdf")), width = 8, height = 6)
  print(plot_volcano(res, comparison_name))
  dev.off()
  
  # Fold-change bar plot
  pdf(output_path(paste0("barplot_", clean_comparison, ".pdf")), width = 8, height = 6)
  print(plot_fc_barplot(res, comparison_name, n_metabolites = 20))
  dev.off()
  
}

# ------ Generate heatmap per comparison ------ #

for(comparison_name in names(res_pairwise_list)){
  
  res <- res_pairwise_list[[comparison_name]]
  clean_comparison <- gsub(" ", "_", comparison_name)
  
  top_metabs <- res %>%
    filter(padj < fdr_threshold) %>%
    arrange(padj) %>%
    slice_head(n = 50) %>%
    pull(feature_id)
  
  if(length(top_metabs) == 0){
    
    cat("No significant metabolites for", comparison_name, "- skipping heatmap.\n")
  }
  
  heat_data <- x_norm[top_metabs, , drop = FALSE]
  
  pdf(output_path(paste0("heatmap_", comparison_name, ".pdf")), width = 12, height = 12)
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

# ------ ANOVA-specific visualizations (only when n_groups > 2) ------#

if(n_groups > 2){
  
  anova_all_path <- intermediate_path("anova_all_results.csv")
  anova_sig_path <- intermediate_path("anova_significant_results.csv")
  
  if(!file.path(anova_all_path) || !file.exists(anova_sig_path)){
    
    stop(
      paste0(
        "ANOVA result files not found in ", intermediate_dir,
        "\nRun statistical_analysis.R first."
      )
    )
    
  }
  
  anova_res <- read.csv(anova_all_path, stringsAsFactors = FALSE)
  anova_sig_res <- read.csv(anova_sig_path, stringsAsFactors = FALSE)
  
  cat("\n===== ANOVA Summary =====\n")
  cat("Total metabolites tested: ", nrow(anova_res), "\n")
  cat("Significant metabolites (FDR <", fdr_threshold, "):", nrow(anova_sig_res), "\n")
  cat("Percentage signficant: ",
      round(nrow(anova_sig_res)/nrow(anova_res) * 100, 1), "%\n")
  
  # Heatmap of top 30 ANOVA-significant metabolites
  top_anova_metabolites <- anova_sig_res %>%
    arrange(padj) %>%
    slice_head(n = 30) %>%
    pull(feature_id)
  
  pdf(output_path("heatmap_anova_top30.pdf"), width =14, height=10)
  print(pheatmap(
    x_norm[top_anova_metabolites, , drop = FALSE],
    scale = "row",
    annotation_col = annotation_col,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "average",
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 8,
    main = "Top 30 ANOVA Significant Metabolites"
  ))
  
  dev.off()
  
  # Boxplots of top 6 ANOVA-significant metabolites
  top6_anova <- anova_sig_res %>%
    arrange(padj) %>%
    slice_head(n=6) %>%
    pull(feature_id)
  
  plot_data_anova <- x_norm[top6_anova, , drop=FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature_id") %>%
    pivot_longer(-feature_id, names_to = "sample", values_to = "abundance") %>%
    left_join(meta, by = c("sample" = "Sample"))
  
  pdf(output_path("anova_boxplots_top6.pdf"), width = 8, height = 6)
  print(ggplot(plot_data_anova, aes(x = Group, y = abundance, fill = Group)) +
          geom_boxplot(outlier.shape = NA) +
          geom_point(position = position_jitter(width = 0.2), alpha = 0.5, size = 2) +
          facet_wrap(~ feature_id, scales = "free_y, ncol = 3") +
          labs(
            title = "Top 6 metabolites by ANOVA",
            y = "Log2 Normalized Abundance",
            x = ""
          ) +
          
          theme_bw() +
          theme(
            legend.position = "none",
            strip.text = element_text(size = 10, face = "bold"),
            plot.title = element_text(hjust = 0.5)
          ))
  
  dev.off()
  
}

cat("\nVisualization complete. Output written to:", output_dir, "\n")












