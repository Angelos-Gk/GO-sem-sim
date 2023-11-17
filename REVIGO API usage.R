# setwd("C:/Users/angel/Desktop/Applied Bioinformatics/thesis/pleiotropy/R")

start.time <- Sys.time()

#===================== SETTING UP ===============================

# Load required libraries
library(GO.db)
library(AnnotationDbi)
library(ontologyIndex)
library(ontologySimilarity)
library(GOSemSim)
library(progress)

# Load the GO OBO file and create a GO graph (https://current.geneontology.org/ontology/go-basic.obo)
# obo file info: 
# -format-version: 1.2
# -data-version: releases/2023-06-11
go_obo_file <- "C:/Users/angel/Desktop/Applied Bioinformatics/Thesis/pleiotropy/data/go-basic.obo"
go <- get_ontology(go_obo_file)

# Load the GAF annotation file (https://current.geneontology.org/products/pages/downloads.html)
# list of gaf files:
# -Caenorhabditis elegans: wb.gaf
# -Danio rerio: zfin.gaf
# -Drosophila melanogaster: fb.gaf
# -Escherichia coli: ecocyc.gaf
# -Gallus gallus: goa-chicken.gaf
# -Homo sapiens: goa-human.gaf
# -Mus musculus: mgi.gaf
# -Saccharomyces cerevisiae: sgd.gaf

annotation_file <- "C:/Users/angel/Desktop/Applied Bioinformatics/Thesis/pleiotropy/data/goa_human.gaf"

# Read the GAF file while skipping header lines
lines <- readLines(annotation_file)
data_start <- grep("^!", lines, invert = TRUE)
annotations <- read.table(text = lines[data_start], sep = "\t", quote = "", comment.char = "", 
                          col.names = c('DB', 'DB_Object_ID', 'DB_Object_Symbol', 'Qualifier',
                                        'GO_ID', 'DB_Reference', 'Evidence_Code', 'With_From',
                                        'Aspect', 'DB_Object_Name', 'DB_Object_Synonym',
                                        'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By',
                                        'Annotation_Extension', 'Gene_Product_Form_ID'),
                          stringsAsFactors = FALSE)

#===================== FILTERING===============================

#NA ELEGXO PANTA POSA PROTEIN EXEI TO KATHE SPECIES KAI META NA DIALEGO DB Objects Types

# 1) Filter rows based on criteria: 'NOT' in Qualifier, Aspect not equal to 'C' or 'F'
filtered_annotations <- annotations[!(grepl("NOT", annotations$Qualifier) | annotations$Aspect %in% c("C", "F")), ]

# keep only the rows that in the column DB Object Type have either the value "protein","protein_coding_gene" or "gene_product"
# some organisms have only the "protein" DB Object Type
filtered_annotations <- filtered_annotations[filtered_annotations$DB_Object_Type %in% c("protein"), ]

# 2) Group annotations by DB Object ID
grouped_annotations <- split(filtered_annotations, filtered_annotations$DB_Object_ID)

# 3) check for duplicates

# Search for the annotations associated with the given ID (we see that this ID has duplicates)
grouped_annotations$A0A087WXM9

# Create an empty list to store unique annotations for each ID
clean_data <- list()

# Iterate through each ID and remove duplicate GO terms
for (target_id in names(grouped_annotations)) {
  annotations_for_id <- grouped_annotations[[target_id]]
  
  # Remove duplicates based on 'GO_ID'
  unique_annotations_for_id <- annotations_for_id[!duplicated(annotations_for_id$GO_ID), ]
  
  clean_data[[target_id]] <- unique_annotations_for_id
}

clean_data$A0A087WXM9

#===================== REVIGO API usage ================

library(httr)
library(rvest)
library(dplyr)
library(progress)

# Function to process a dataframe and get filtered GO terms and Dispensability values
process_dataframe <- function(df, cutoff = 0.7) {
  go_ids <- df$GO_ID
  
  # Convert GO_IDs to a string
  go_list <- paste(go_ids, collapse = "\n")
  
  # Submit job to Revigo
  res <- httr::POST(
    url = "http://revigo.irb.hr/Revigo",
    body = list(
      cutoff = as.character(cutoff),
      valueType = "pvalue",
      speciesTaxon = "0",
      measure = "SIMREL",
      goList = go_list
    ),
    encode = "form"
  )
  
  # Extract Dispensability values from the resulting HTML
  dat <- httr::content(res, as = "text", encoding = "UTF-8")
  
  # Parse HTML using rvest
  html <- read_html(dat)
  
  # Extract Dispensability values and GO terms from the table
  rows <- html_nodes(html, xpath = "//table[@id='BiologicalProcess']/tr[position()>1]") 
  
  data <- lapply(rows, function(row) {
    row_values <- html_nodes(row, xpath = "./td") %>% html_text()
    as.list(row_values)
  }) %>% 
    do.call(rbind, .) %>%
    as.data.frame()
  
  # Handle cases where column names are not recognized
  colnames(data) <- c("Term ID", "Name", "Frequency", "Value", "Uniqueness", "Dispensability", "Eliminated", "Representative")
  
  # Filter based on cutoff
  filtered_data <- data[data$Dispensability < cutoff, ]
  
  # Return a data frame with GO terms and Dispensability values
  return(filtered_data)
}

# Initialize the progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent :current/:total dataframes processed, Elapsed: :elapsedfull ETA: :eta",
  total = length(clean_data)
)

# Process all dataframes in clean_data
results <- lapply(names(clean_data), function(df_name) {
  df <- clean_data[[df_name]]
  result <- process_dataframe(df)
  pb$tick()  # Increment progress bar
  return(result)
})

# Terminate the progress bar
pb$terminate()

# Print the resulting list of data frames
print(results)
