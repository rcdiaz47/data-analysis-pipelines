# ----- Expression Matrix Builder (for Omics Studio) ----- # 
# Purpose: Omics studios expression matrix wants normalized values on the
# LINEAR scale, not log transformed, with a unique molecule ID per row.
# This runs on raw peak areas: filter -> median normalization, attach InChIKeys from
# map_metabolite_ids.R output

library(readxl)
library(dplyr)


# ==========================================
# USER CONFIGURATION - Edit this section only 
# ===========================================

input_file <- "2Group_test_metabolomics.xlsx"
id_mapping_file <- "C:\\Users\\Proteomics\\Desktop\\results\\id_mapping.csv"

name_column <- "Name"
sample_area_prefix <- "Group Area:"

output_dir <- "./results"
expression_output_file <- "expression_matrix.csv"

# Drop metabolites missing in more than this fraction of samples 
max_missing_fraction <- 0.30 

# ===========================================
# END USER CONFIGURATION 
# ===========================================

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

if(!file.exists(input_file)) stop(paste0("Input file not found :", input_file))
if(!file.exists(id_mapping_file)){
  stop(paste0(
    "ID mapping file not found: ", id_mapping_file,
    "\nRun map_metabolite_ids.R first to generate it"
  ))
}

cd <- read_xlsx(input_file)
id_table <- read.csv(id_mapping_file, stringsAsFactors = FALSE)

if(nrow(cd) != nrow(id_mapping_file)){
  stop("Row count mismatch between input_file and id_mappping_file - ",
       "Make sure both were generated from the same source spreadsheet, ",
       "in the same row order.")
}

# ----- Build the numeric matrix of peak areas ----- #
area_cols <- grep(paste0("^", sample_area_prefix), colnames(cd), value = TRUE)
area_cols <- area_cols[!grep("Blank|blank|BLANK|QC|Qc|qc", area_cols)]

if(length(area_cols) == 0){
  stop(paste0("No Sample area columns found with the prefix '", sample_area_prefix, "'."))
}

x <- cd %>% dplyr::select(all_of(area_cols)) %>% as.matrix()
x[x == ""] <- NA
storage.mode(x) <- "numeric"
















