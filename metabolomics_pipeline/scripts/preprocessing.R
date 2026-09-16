library(tidyverse)
library(readxl)
library(yaml)

# read the pipeline configuration
config <- yaml::read_yaml("C:\\Users\\Proteomics\\Documents\\data-analysis-pipelines\\metabolomics_pipeline\\config.yaml")

# assign configuration values 
input_file <- config$input_file 
n_groups <- config$n_groups 
group_names <- unlist(config$group_names)
output_dir <- config$output_dir
missingness_threshold <- config$missingness_threshold
name_column <- config$name_column
sample_area_prefix <- config$sample_area_prefix 

#--------- validate configuration -----------------

if(n_groups < 2){
  stop("n_groups must be at least 2.")
}

if(length(group_names) != n_groups ){
  stop(
    paste0(
      "n_groups is", n_groups,
      ", but", length(group_names),
      " group names were provided."
    )
  )
}

if(anyDuplicated(group_names)){
  stop("group names must contain unique group names.")
}

if(missingness_threshold < 0 || missingness_threshold > 1){
  stop("missingness threshold must be between 0 and 1.")
}

if(!file.exists(input_file)){
  stop(
    paste0(
      "Input file not found:",
      input_file,
      "\nPlease check input file in config.yaml."
    )
  )
}

# Create the intermediate results directory
intermediate_dir <- file.path(output_dir, "intermediate")

if(!dir.exists(intermediate_dir)){
  dir.create(intermediate_dir, recursive = TRUE)
}

# Helper function for output paths 
intermediate_path <- function(filename){
  file.path(intermediate_dir, filename)
}


# ----- Read input data --------
cd <- read_xlsx(input_file)

# Validate the metabolite-name column
if(!name_column %in% colnames(cd)){
  stop(
    paste0(
      "The configured metabolite-name column '",
      name_column,
      "' was not found in the input file.\n",
      "avaliable columns are: ",
      paste(colnames(cd),collapse = ",")
    )
  )
}


# ----- Clean metabolite names ----- 

process_metabolite_names <- function(x, max_length = 40){
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  
  ifelse(
    nchar(x) > max_length, paste0(substr(x,1,max_length-3), "..."), x
  )
}

cd$feature_id <- cd[[name_column]]
cd$display_name <- process_metabolite_names(cd$feature_id)

# ----- Identify sample-area columns ------- 

area_cols <- grep(paste0("^", sample_area_prefix),
                  colnames(cd),
                  value = TRUE
)

# exclude blanks and qcs
area_cols <- area_cols[!grepl("blank|QC", area_cols, ignore.case = TRUE)]

# Confirm that sample columns were found
if(length(area_cols) == 0){
  stop(
    paste0(
      "No sample-area columns were found.\n",
      "Expected column names beginning with '",
      sample_area_prefix,
      "'."
    )
  )
}

# ------ Build and validate the numeric matrix ------- 

# Confirm that all the selected sample columns are numeric
non_numeric_columns <- area_cols[!vapply(cd[area_cols], is.numeric, logical(1))]

if(length(non_numeric_columns) > 0){
  stop(
    paste0(
      "The following sample-area columns are not numeric:\n",
      paste(non_numeric_columns, collapse = ", ")
    )
  )
}

# Build the numeric matrix 
x <- cd %>%
  dplyr::select(all_of(area_cols)) %>% as.matrix()

# Create unique identifiers for duplicates without changing the original names
feature_key <- make.unique(
  as.character(cd$feature_id),
  sep = "_duplicate_"
)


rownames(x) <- feature_key

# Validate the completed matrix 
if(any(colSums(!is.na(x)) == 0)){
  empty_samples <- colnames(x)[colSums(!is.na(x)) == 0]
  
  stop(
    paste0(
      "The following samples contain no numeric values:\n",
      paste(empty_samples, collapse = ", ")
    )
  )
}


# ------ Identify and export internal standard raw signal for QC ------ #

is_pattern <- "-d[0-9]+$"
is_rows <- grepl(is_pattern, trimws(cd$feature_id))

if(any(is_rows)){
  
  is_data <- x[is_rows, drop = FALSE]
  rownames(is_data) <- feature_key[is_rows]
  
  is_export <- is_data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature_id")
  
  write.csv(is_export, intermediate_path("internal_standard_qc.csv"), row.names = FALSE)
  cat("Internal standard QC: found", sum(is_rows), "IS feature(s), exported to internal_standard_qc.csv\n")
  
} else{
  cat("Internal Standard QC: no features matched pattern '", is_pattern, "' - skipping\n")
}


# ----- Load the sample -> group mapping file ----- 

group_mapping_file <- config$group_mapping_file 

if(!file.exists(group_mapping_file)){
  stop(
    paste0(
      "Group mapping file not found: ",
      group_mapping_file,
      "\nExpected a CSV with columns 'sample' and 'group'."
    )
  )
}

group_map <- read.csv(group_mapping_file, stringsAsFactors = FALSE)

if(!all(c("sample", "group") %in% colnames(group_map))){
  stop("group_mapping_file must contain 'sample' and 'group' columns.")
}

# Confirm every sample-area column has a matching entry in the mapping file
missing_from_map <- setdiff(area_cols, group_map$sample)

if(length(missing_from_map) > 0){
  stop(
    paste0(
      "The following sample columns have no entry in the group mapping file:\n",
      paste(missing_from_map, collapse = ", ")
    )
  )
}

group_map$sample <- trimws(group_map$sample)
group_map$group <- trimws(group_map$group)
group_names <- trimws(group_names)

# Confirm every group in the mapping file is one of the configured group_names
unknown_groups <- setdiff(unique(group_map$group), group_names)

if(length(unknown_groups) > 0){
  stop(
    paste0(
      "The group mapping gile referencens group(s) not listed in config.yaml group_names:\n",
      paste(unknown_groups, collapse = ", ")
    )
  )
}

# Build sample metadata in the same order as the columns of x
meta <- data.frame(
  Sample = area_cols,
  Group = group_map$group[match(area_cols, group_map$sample)],
  stringsAsFactors = FALSE
)

# ----- Log2 transform ------ #
x[x == ""] <- NA
x_log <- log2(x)

# ----- Filter metabolites by missingness ------ #
# NOTE : despite the name, missingness_threshold is used as a PRESENCE-rate cutoff here

# Keep features present in at least missingness_threshold of samples

keep <- rowMeans(!is.na(x_log)) >= missingness_threshold

x_filt <- x_log[keep, ]

# ----- Normalize using median centering ----- #
column_median <- apply(x_filt, 2, median, na.rm = TRUE)
x_norm <- sweep(x_filt, 2, column_median, FUN = "-")

# ----- Export normalized matrix and sample metadata ----- #
x_norm_export <- x_norm %>% as.data.frame() %>% tibble::rownames_to_column("feature_id")

write.csv(x_norm_export, intermediate_path("normalized_matrix.csv"), row.names = FALSE)

write.csv(meta, intermediate_path("sample_metadata.csv"), row.names = FALSE)

cat("Preprocessing complete.\n")
cat("- Metabolites before filtering:", nrow(x_log), "\n")
cat("- Metabolites after filtering (>=", missingness_threshold * 100, "% present):",
    nrow(x_norm), "\n")
cat("- Samples:", ncol(x_norm), "\n")
cat(" - Output written to:", intermediate_dir, "\n")













































