
#================= SAVE THE OUTPUT =========

# Initialize an empty character vector to store the results
results <- character()

# Loop through each element (character vector) in the list and extract information
for (i in seq_along(run1)) {
  go_ids <- run1[[i]]
  df_name <- names(run1)[i]
  go_ids_combined <- paste(go_ids, collapse = ";")
  go_id_count <- length(go_ids)
  result <- paste(df_name, go_id_count, go_ids_combined, sep = "\t")
  results <- c(results, result)
}

# Specify the path and filename for the text file ("organism_threshold_method.txt")
output_file <- "human_7_wang.txt"

# Write the results to the text file
writeLines(results, output_file)

#================= CALCULATE NUMBER OF TERMS BEFORE AND AFTER REDUCTION =========

# BEFORE REDUCTION

# Initialize a list to store the counts
go_id_counts <- list()

# Calculate the number of GO IDs for each sub-dataframe and store in the list
for (id in names(clean_data)) {
  go_id_counts[[id]] <- nrow(clean_data[[id]])
}

# Calculate the total number of GO IDs in the clean_data object
total_go_ids <- sum(sapply(go_id_counts, sum))

# Print the counts and total
cat("Total number of GO IDs (before correction):", total_go_ids)

# AFTER REDUCTION

# Initialize a list to store the counts
filtered_go_id_counts <- list()

# Calculate the number of GO IDs for each sub-list and store in the list
for (id in names(run2)) {
  filtered_go_id_counts[[id]] <- length(run1[[id]])
}

# Calculate the total number of GO IDs in the filtered_terms_list object
total_filtered_go_ids <- sum(sapply(filtered_go_id_counts, sum))

# Print the counts and total
cat("Total number of GO IDs (after correction):", total_filtered_go_ids)
