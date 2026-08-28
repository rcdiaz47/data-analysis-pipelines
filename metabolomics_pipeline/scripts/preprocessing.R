library(tidyverse)
library(readxl)
library(yaml)

# read the pipeline configuration
config <- yaml::read_yaml("config.yaml")

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




















