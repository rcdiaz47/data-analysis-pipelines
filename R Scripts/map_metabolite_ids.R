# ------ Metabolite ID Mapping (SMILES -> InChI) ----- #
# Purpose: Omics Studio requires molecule IDs to be Uniprot, Ensembl, CHEBI,
# or InCHI (Only one type per file), Our data has Name, Formula, CSID, and
# SMILES but no native InChI/CHEBI column, so we derive InChI locally from 
# SMILES using rcdk. This is a structural conversion, deterministic, not a name 
# Search, so theres no ambiguity the way name-matching would be 

library(rcdk)
library(dplyr)
library(readxl)
library(httr)
library(jsonlite)
#==============================================
# USER CONFIGURATION - Edit this section only
#==============================================

input_file <- "C:\\Users\\Proteomics\\Desktop\\2Group_test_metabolomics.xlsx"

name_columns <- "Name"
formula_column <- "Formula"
csid_column <- "CSID"
smiles_column <- "SMILES"

output_dir <- "./results"
mapping_output_file <- "id_mapping.csv"

#============================================
# END USER CONFIGURATION
#============================================

if(!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}

if(!file.exists(input_file)){
  stop(paste0("Input file not found: ", input_file))
}

cd <- read_xlsx(input_file)

required_cols <- c(name_columns, formula_column, csid_column, smiles_column)
missing_cols <- setdiff(required_cols, colnames(cd))

if(length(missing_cols) > 0){
  stop(paste0(
    "Missing expected column(s): ", paste(missing_cols, collapse = ", "),
    "\nAvaliable columns are:", paste(colnames(cd)), collapse = " ,"
  ))
}

# ----- Convert each SMILES string to InChI ----- #

smiles_vec <- cd[[smiles_column]]

get_inchi_from_smiles <- function(smiles){
  
  if(is.na(smiles) || trimws(smiles) == ""){
    return(NA_character_)
  }
  
  #URL encode the smiles strings so special characters dont break the request URL
  encoded_smiles <- URLencode(smiles, reserved = TRUE)
  
  url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/",
    encoded_smiles,
    "/property/InChI/JSON"
  )
  
  response <- tryCatch(GET(url), error = function(e) NULL)
  
  if(is.null(response) || status_code(response) != 200){
    
    return(NA_character_)
    
  }
  
  parsed <- tryCatch(fromJSON(content(response, as = "text", encoding = "UTF-8")),
                     error = function(e) NULL)
  
  if(is.null(parsed)) return(NA_character_)
  
  inchi_result <- parsed$PropertyTable$Properties$InChI
  
  # If PubChem returned a response but it didnt contain an InChI value,
  # for this compound, this will be length 0- so return NA instead of trying to index
  # nothing 
  
  if(length(inchi_result) == 0) return(NA_character_)
  
  inchi_result[1]
  
}

inchi_vec <- vapply(smiles_vec, get_inchi_from_smiles, character(1))

id_table <- cd %>%
  dplyr::select(all_of(required_cols)) %>%
  dplyr::mutate(InChI = inchi_vec)

# ----- QC Checks ----- #

n_total <- nrow(id_table)
n_failed <- sum(is.na(id_table$InChI))
n_duplicated <- sum(duplicated(id_table$InChI[!is.na(id_table$InChI)]))

cat("Total metabolites:", n_total, "\n")
cat("Failed SMILES to InChI conversion :", n_failed, "\n")
cat("Duplicated InChI values (same structure, different rows) :", n_duplicated, "\n")


if(n_failed > 0){
  cat("\nRows with failed conversions :\n")
  print(id_table[is.na(id_table$InChI), c(name_columns, smiles_column)])
}

if(n_duplicated > 0){
  dup_inchi <- id_table$InChI[duplicated(id_table$InChI) & !is.na(id_table$InChI)]
  cat("\nRpws sharing a duplicated InChI :\n")
  cat("e.g stereoisomers collapsing to one structure :\n")
  print(id_table[id_table$InChI %in% dup_inchi & !is.na(id_table$InChI),  ])
}

write.csv(id_table, file.path(output_dir, mapping_output_file), row.names = FALSE)









