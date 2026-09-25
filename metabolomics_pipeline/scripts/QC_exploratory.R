# ------ Quality Control / Exploratory Analysis ----- #
library(tidyverse)
library(yaml)
library(ggplot2)
library(ggrepel)

# read the pipeline configuration
config <- yaml::read_yaml("C:\\Users\\Proteomics\\Documents\\data-analysis-pipelines\\metabolomics_pipeline\\config.yaml")

output_dir <- config$output_dir

intermediate_dir <- file.path(output_dir, "intermediate")

# Helper function for output paths
output_path <- function(filename){
  file.path(output_dir, filename)
}

# ----- Read preprocessed data ----- #
norm_matrix_path <- file.path(intermediate_dir, "normalized_matrix.csv")
metadata_path <- file.path(intermediate_dir, "sample_metadata.csv")

if(!file.exists(norm_matrix_path)){
  stop(
    paste0(
      "Normalized matrix not found: ", norm_matrix_path,
      "\nRun preprocessing.R first."
    )
  )
}

if(!file.exists(metadata_path)){
  stop(
    paste0(
      "Sample metadata not found: ", metadata_path,
      "\nRun preprocessing.R first."
    )
  )
}

x_norm_df <- read.csv(norm_matrix_path, stringsAsFactors = FALSE, check.names = FALSE)
meta <- read.csv(metadata_path, stringsAsFactors = )

# Rebuild the matrix with feature_id as rownames
x_norm <- x_norm_df %>% column_to_rownames("feature_id") %>% as.matrix()

# ------ Density plot of normalized data ----- #
norm_df <- x_norm %>% as.data.frame() %>% tibble::rownames_to_column("feature_id") %>%
  pivot_longer(-feature_id, names_to = "sample", values_to = "value")

density_plot <- ggplot(norm_df, aes(x=value)) +
  geom_density(fill = "blue", alpha = 0.4) +
  labs(
    title = "Normalized Data Distribution",
    x = 'Normalized Intensity (log2)',
    y = "Density"
  ) +
  theme_minimal()
pdf(output_path("normalized_density_plot.pdf"), width = 8, height = 6)
print(density_plot)
dev.off()

# ------ PCA ------#

# Transpose so samples are rows and features are columns; drops rows with any NA
x_pca <- t(na.omit(x_norm))

pca_res <- prcomp(x_pca, scale. = TRUE)

pca_scores <- as.data.frame(pca_res$x)

# Match samples to their group using metadata (order-safe)
pca_scores$Sample <- rownames(pca_scores)

pca_scores <- pca_scores %>% left_join(meta, by = "Sample")

pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample), size = 2.5, show.legend = FALSE) +
  stat_ellipse(
    aes(fill = Group),
    geom = "polygon",
    alpha = 0.2,
    color = NA
  ) +
  labs(
    title = "PCA Scores Plot",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 1), "%)")
  ) + theme_minimal()

pdf(output_path("pca_plot.pdf"), width = 8, height = 6)
print(pca_plot)
dev.off()

# ------ PC1 by batch ------ #
# Direct check: does PC1 differ between the two batched?

sample_number <- as.numeric(gsub("Area:\\s*(\\d+)_.*", "\\1", meta$Sample))
meta$Batch <- ifelse(sample_number <=50, "Batch 1 (39-50)", "Batch 2 (51-60)")

pca_scores$Batch <- meta$Batch[match(pca_scores$Sample, meta$Sample)]

pca1_batch_plot <- ggplot(pca_scores, aes(x = reorder(Sample, PC1), y = PC1, fill = Batch)) +
  geom_col() +
  labs(title = "PC1 Score by sample and Batch", x = "", y = "PC1") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

pdf(output_path("pc1_by_batch.pdf"), width = 10, height = 6)
print(pca1_batch_plot)
dev.off()
  

# ------ Internal standard QC ------ #

is_qc_path <- file.path(intermediate_dir, "internal_standard_qc.csv")

if(file.exists(is_qc_path)){
  is_df <- read.csv(is_qc_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  is_long <- is_df %>%
    pivot_longer(-feature_id, names_to = "sample", values_to = "intensity") %>%
    left_join(meta, by = c("sample" = "Sample"))
  
  is_plot <- ggplot(is_long, aes(x = sample, y = intensity, color = Group)) +
    geom_point(size = 2) +
    facet_wrap(~feature_id, scales = "free_y", ncol = 2) +
    labs(title = "Internal Standard QC", x = "", y = "Raw Peak Area") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
  
  pdf(output_path("internal_standard_qc.pdf"), width = 10, height = 8)
  print(is_plot)
  dev.off()
  
  is_cv <- is_long %>%
    group_by(feature_id) %>%
    summarise(
      mean_intensity = mean(intensity, na.rm = TRUE),
      sd_intensity = sd(intensity, na.rm = TRUE),
      cv_percent = round(100*sd_intensity/mean_intensity, 1)
    )
  
  write.csv(is_cv, output_path("internal_standard_cv.csv"), row.names = FALSE)
  cat("Internal standard QC: CV summary written to internal_standard_cv.csv\n")
  print(is_cv)
  
}else {
  
  cat("No internal standard QC file founf- skipping IS QC plot\n")
  
}


cat("QC / exploratory analysis complete.\n")
cat(" - Features used in PCA (complete cases):", ncol(x_pca), "\n")
cat(" - Samples:", nrow(pca_scores), "\n")
cat(" - Output written to:", output_dir, "\n")
























