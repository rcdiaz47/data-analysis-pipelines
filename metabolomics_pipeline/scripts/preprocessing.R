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

process_metabolite_names <- function










