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

# Load the GO OBO file and create a GO graph
go_obo_file <- "C:/Users/angel/Desktop/Applied Bioinformatics/Thesis/pleiotropy/python/go-basic.obo"
go <- get_ontology(go_obo_file)

# Load the GAF annotation file
annotation_file <- "C:/Users/angel/Desktop/Applied Bioinformatics/Thesis/pleiotropy/python/goa_human.gaf"

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

# 1) Filter rows based on criteria: 'NOT' in Qualifier, Aspect not equal to 'C' or 'F'
filtered_annotations <- annotations[!(grepl("NOT", annotations$Qualifier) | annotations$Aspect %in% c("C", "F")), ]

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

#===================== GET GO TERMS LEVEL INFO ===============================

# Read GO terms with levels from the file
go_terms_with_level <- read.table("go_terms_with_levels.txt", sep="\t", col.names=c("GO.Term", "Level"), stringsAsFactors=FALSE)

# Convert GO terms to a named vector of levels
go_levels <- setNames(go_terms_with_level$Level, go_terms_with_level$GO.Term)

# Create a function to get the level of a GO term
get_level <- function(go_term) {
  if (go_term %in% names(go_levels)) {
    return(go_levels[[go_term]])
  } else {
    return(NA)  # Return NA if level is not found
  }
}

#===================== SEMANTIC SIMILARITY ===============================

# Create GOSemSimData object 
go_sem_data <- godata("org.Hs.eg.db", ont = "BP", computeIC = TRUE)

# Create a new object called clean_data_sub and subset only the desired dataframes (quick testing)
#clean_data_sub <- clean_data[c("A0A087WXM9", "A0A024RBG1", "A0A075B6H5")]

# Create a new list to store filtered GO terms for each ID
filtered_terms_list <- list()

# Initialize the progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent :current/:total ids processed, Elapsed: :elapsedfull ETA: :eta",
  total = length(names(clean_data))
)

# Get the start time
start_time <- Sys.time()

# Iterate through each ID and filter GO terms
for (target_id in names(clean_data)) {
  annotations_for_id <- clean_data[[target_id]]
  
  # If there's only one GO term, keep it
  if (nrow(annotations_for_id) == 1) {
    filtered_terms_list[[target_id]] <- annotations_for_id$GO_ID
    pb$tick()  # Update the progress bar
    next
  }
  
  # Calculate semantic similarity between terms of the same ID
  go_terms <- annotations_for_id$GO_ID
  n <- length(go_terms)
  sim_matrix <- matrix(0, nrow = n, ncol = n)
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      sim <- goSim(go_terms[i], go_terms[j], semData = go_sem_data, measure = "Wang")
      sim_matrix[i, j] <- sim
      #sim_matrix[j, i] <- sim
    }
  }
  
  # Apply the reduction criteria
  kept_terms <- character(0)
  used_indices <- c()
  for (i in 1:n) {
    if (i %in% used_indices) {
      next
    }
    
    sim_row <- sim_matrix[i, ]
    similar_indices <- which(sim_row >= 0.4)
    similar_levels <- sapply(similar_indices, function(index) go_levels[go_terms[index]])
    similar_and_same_level <- similar_indices[similar_levels == go_levels[go_terms[i]]]
    
    if (length(similar_and_same_level) > 0) {
      # Set a seed for reproducibility
      #set.seed(124)
      selected_index <- sample(c(i, similar_and_same_level), 1)
      used_indices <- c(used_indices, similar_and_same_level[similar_and_same_level != selected_index])
      kept_terms <- c(kept_terms, go_terms[selected_index])
    } else {
      kept_terms <- c(kept_terms, go_terms[i])
    }
  }
  
  # Remove duplicates and store the kept terms
  filtered_terms_list[[target_id]] <- unique(kept_terms)
  
  pb$tick()  # Update the progress bar
}

# Terminate the progress bar
pb$terminate()

#===================== TESTING =============================== 

# In order to test, paste the clean data to revigo (http://revigo.irb.hr/) and compare results with 
# the reduced_grouped_annotations results for each ID

# Testing
clean_data$A0A087WXM9$GO_ID
filtered_terms_list$A0A087WXM9

# Testing 
clean_data$A0A024RBG1$GO_ID
filtered_terms_list$A0A024RBG1

# Testing
clean_data$A0A075B6H5$GO_ID
filtered_terms_list$A0A075B6H5

sessionInfo()

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
