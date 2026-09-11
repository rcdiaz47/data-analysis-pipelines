# ------ Statistical Analysis ------ #
library(tidyverse)
library(yaml)

# read the pipeline configuration
config <- yaml::read_yaml("C:\\Users\\Proteomics\\Documents\\data-analysis-pipelines\\metabolomics_pipeline\\config.yaml")

output_dir <- config$output_dir
n_groups <- config$n_groups
group_names <- unlist(config$group_names)

intermediate_dir <- file.path(output_dir, "intermediate")

# Helper function for intermediate output paths
intermediate_path <- function(filename){
  file.path(intermediate_dir, filename)
}

# Fixed thresholds (not user configurable)
fdr_threshold <- 0.05
logfc_threshold <- 1

# ------ Read preprocessed data ------ #
norm_matrix_path <- intermediate_path("normalized_matrix.csv")
metadata_path <- intermediate_path("sample_metadata.csv")

if(!file.exists(norm_matrix_path)){
  stop(
    paste0(
      "Normalized matrix not found: ", norm_matrix_path,
      "\nRun preprocessing.R first"
    )
  )
}

if(!file.exists(metadata_path)){
  stop(
    paste0(
      "Sample metadata not found: ", metadata_path,
      "\nRun preprocessing.R first"
    )
  )
}

x_norm_df <- read.csv(norm_matrix_path, stringsAsFactors = FALSE)
meta <- read.csv(metadata_path, stringsAsFactors = FALSE)

x_norm <- x_norm_df %>%
  column_to_rownames("feature_id") %>%
  as.matrix()
group <- factor(meta$Group)

# ------ Reusable function for pairwise analysis of two groups ------ #
run_pairwise <- function(group1, group2, x_mat, meta){
  
  samples_keep <- meta$Sample[meta$Group %in% c(group1, group2)]
  groups_keep <- meta$Group[meta$Group %in% c(group1, group2)]
  
  x_sub <- x_mat[, samples_keep]
  
  res <- apply(x_sub, 1, function(z){
    
    df <- data.frame(value = z, group = groups_keep)
    df <- df[!is.na(df$value), ]
    
    if(length(unique(df$group)) < 2) {
      return(c(logFC =NA, pvalue = NA))
    }
    
    m1 <- mean(df$value[df$group == group1], na.rm = TRUE)
    m2 <- mean(df$value[df$group == group2], na.rm = TRUE)
    logFC <- m1 - m2
    
    p <- t.test(value ~ group, data = df)$p.value
    
    c(logFC = logFC, pvalue = p)
    
  })
  
  res <- as.data.frame(t(res))
  res$feature_id <- rownames(x_mat)
  res$padj <- p.adjust(res$pvalue, method = "BH")
  
  res
  
  
}


# ------ Run t-test (2 group) or ANOVA + Tukey + pairwise (3+ groups) ------ #

if(n_groups == 2){
  res_pairwise_list <- list(
    run_pairwise(group_names[1], group_names[2], x_norm, meta)
  )
  
  names(res_pairwise_list) <- paste(group_names[1], "vs", group_names[2])
  
} else{
  
  # ------ One Way Anova per metabolite ------ #
  anova_res <- apply(x_norm, 1, function(z){
    
    df <- data.frame(value = z, group = group)
    df <- df[!is.na(df$value), ]
    
    if(nrow(df) < 3 || length(unique(df$group)) < 2 || sum(!is.na(df$value)) < 3){
      
      return(data.frame(pvalue = NA))
      
    }
    
    fit <- aov(value ~ group, data = df)
    p <- summary(fit)[[1]][["Pr(>F)"]][[1]]
    
    data.frame(pvalue = p)
    
  })
  
  anova_res <- do.call(rbind, anova_res)
  anova_res$feature_id <- rownames(x_norm)
  rownames(anova_res) <- NULL
  
  anova_res$padj <- p.adjust(anova_res$pvalue, method = "BH")
  
  anova_sig_features <- anova_res$feature_id[anova_res$padj < fdr_threshold]
  anova_significant_features_clean <- anova_sig_features[!is.na(anova_sig_features)]
  
  # ------ Tukey HSD post hoc ------
  tukey_list <- lapply(anova_significant_features_clean, function(fid){
    z <- x_norm[fid,]
    df <- data.frame(value = z, group = group)
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
  tukey_sig <- tukey_res %>% filter('p adj' < fdr_threshold)
  
  # ------ All pairwise comparisons ------ #
  pairs <- combn(group_names, 2, simplify = FALSE)
  
  res_pairwise_list <- lapply(pairs, function(pair){
    
    run_pairwise(pair[1], pair[2], x_norm,meta)
    
  })
  
  names(res_pairwise_list) <- sapply(pairs, function(pair){
    paste(pair[1], "vs", pair[2])
  })
  
  
}

# ------ Add significance labels to each pairwise comparison ------- #
res_pairwise_list <- lappy(res_pairwise_list, function(res){
  
  res$significance  <- "Not Significant"
  res$significance[res$padj < fdr_threshold & res$logFC > logfc_threshold] <- "Up"
  res$significance[res$padj < fdr_threshold & res$logFC < -logfc_threshold] <- "Down"
  
  
})

# ------ Drop rows with a missing feature name ------ #
# data integrity check

res_pairwise_list <- lapply(res_pairwise_list, function(res){
  
  n_before <- nrow(res)
  res <- res %>% filter(!is.na(feature_id))
  n_dropped <- n_before - nrow(res)
  
  if(n_dropped > 0){
    cat(" Dropped", n_dropped, "row(s) with missing feature_id\n")
  }
  
  res
  
})

# ------ Print summary of statistical testing ------ #
cat("\n======================================================\n")
cat("Statistical Testing Complete\n")
cat("=========================================================\n")
cat("Groups tested:", paste(group_names, collapse = ", "), "\n")
cat("Number of groups:", n_groups, "\n")
cat("Statistical test used:", ifelse(n_groups ==2, "T-test (pairwise)", "Anova + Tukey + Pairwise t-tests"), "\n")

if(n_groups > 2){
  cat("\nAnova Results:\n")
  cat("Total metabolites tested:", nrow(anova_res), "\n")
  cat(" Significant (FDR <", fdr_threshold, "):", length(anova_significant_features_clean), "\n")
  cat(" Tukey significant comparisons:", nrow(tukey_sig), "\n")

}

cat("\nPairwise Comparisons Run:\n")
for(comparison_name in names(res_pairwise_list)){
  
  res <- res_pairwise_list[[comparison_name]]
  cat(" -", comparison_name, ":", sum(res$padj < fdr_threshold, na.rm = TRUE), "significant metabolites\n")
  
}
cat("=====================================================\n")

# ------ Export results to intermediate directory ------ #

for(comparison_name in names(res_pairwise_list)){
  
  clean_comparison <- gsub(" ", "_", comparison_name)
  
  write.csv(
    res_pairwise_list[[comparison_name]],
    intermediate_path(paste0("pairwise_", clean_comparison,"full.csv")),
    row.names = FALSE
  )
  
}

if(n_groups > 2){
  
  write.csv(anova_res, intermediate_path("anova_all_results.csv"), row.names = FALSE)
  write.csv(anova_res %>% filter(padj < fdr_threshold), intermediate_path("anova_significant_results.csv"), row.names = FALSE)
  write.csv(tukey_res,intermediate_path("tukey_all_results.csv"), row.names = FALSE)
  write.csv(tukey_sig, intermediate_path("tukey_significant_results.csv"), row.names = FALSE)
            
}

# save the comparison names so downstream scripts know what to loop over
writeLines(names(res_pairwise_list), intermediate_path("comparison_names.txt"))

cat("\nOutput written to:", intermediate_dir, "\n")














































