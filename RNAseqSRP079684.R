#setwd("C:/Users/angel/Desktop/Applied Bioinformatics/Thesis/RNAseq analysis")
# data download from refinebio (no transformation but quantile normalized)

# Differential Expression Analysis 

# ==== Set the  file paths ====
# Define the file path to the data directory
data_dir <- file.path("data", "SRP079684")

# Define the file path to the data directory
results_dir <- file.path("data", "SRP079684","results")

# Declare the file path to the gene expression matrix file inside directory saved as `data_dir`
data_file <- file.path(data_dir, "SRP079684.tsv")

# Declare the file path to the metadata file inside the directory saved as `data_dir`
metadata_file <- file.path(data_dir, "metadata_SRP079684.tsv")

# ==== Libraries ====
# Attach the DESeq2 library
library(DESeq2)
# Attach the ggplot2 library for plotting
library(ggplot2)
# We will need this so we can use the pipe: %>%
library(magrittr)
# For good quality volcano plots
library (EnhancedVolcano)

# The jitter plot we make later on with the DESeq2::plotCounts() function involves some randomness. 
# As is good practice when our analysis involves randomness, we will set the seed.
set.seed(12345)

# ==== Import data and metadata ====
# Read in metadata TSV file
metadata <- readr::read_tsv(metadata_file)

# Read in data TSV file
expression_df <- readr::read_tsv(data_file) %>%
  tibble::column_to_rownames("Gene")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

head(metadata$refinebio_title)

# ==== Setup metadata ====
# create a new column with cancer/control status of each sample
metadata <- metadata %>%
  # Let's get the RPL10 mutation status from this variable
  dplyr::mutate(status = dplyr::case_when(
    stringr::str_detect(refinebio_title, "tumor") ~ "cancer",
    stringr::str_detect(refinebio_title, "normal") ~ "control"
  ))

dplyr::select(metadata, refinebio_title, status)

# Print out a preview of `mutation_status`
str(metadata$status)
table(metadata$status)
# Make mutation_status a factor and set the levels appropriately
metadata <- metadata %>%
  dplyr::mutate(
    # Here we define the values our factor variable can have and their order.
    status = factor(status, levels = c("control", "cancer"))
  )

levels(metadata$status)

# ==== Define a minimum counts cutoff ====
# Define a minimum counts cutoff and filter the data to include (dont do it if z-score)
# only rows (genes) that have total counts above the cutoff
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 200)

# ==== Create a DESeq2Dataset ====
# round all expression counts
gene_matrix <- round(filtered_expression_df)

ddset <- DESeqDataSetFromMatrix(
  # Here we supply non-normalized count data
  countData = gene_matrix,
  # Supply the `colData` with our metadata data frame
  colData = metadata,
  # Supply our experimental variable to `design`
  design = ~status
)

# ==== Run differential expression analysis ====
deseq_object <- DESeq(ddset)
deseq_results <- results(deseq_object)

deseq_results <- lfcShrink(
  deseq_object, # The original DESeq2 object after running DESeq()
  coef = 2, # The log fold change coefficient used in DESeq(); the default is 2.
  res = deseq_results # The original DESeq2 results table
)

head(deseq_results)
summary(deseq_results)

# this is of class DESeqResults -- we want a data frame
deseq_df <- deseq_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  # higher expression in RPL10 mutated samples
  dplyr::arrange(dplyr::desc(log2FoldChange))

head(deseq_df)

# ==== Count number of significant and non significant ====
# Count the number of rows with TRUE and FALSE in the threshold column
threshold_count <- table(deseq_df$threshold)

# Print the result
print(threshold_count)

# ==== Check results by plotting one gene ====
plotCounts(ddset, gene = "ENSG00000141367", intgroup = "status")

# ==== Create a Volcano plot # ==== 
# We'll assign this as `volcano_plot`
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  title = "Dataset: SRP079684",
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05 # Loosen the cutoff since we supplied corrected p-values
)

# Print out plot here
volcano_plot

# ==== Save results to TSV ====
# Filter for differentially expressed genes based on both threshold and log2FoldChange
DE_gene_ids <- deseq_df$Gene[!is.na(deseq_df$padj) & deseq_df$padj < 0.05 & abs(deseq_df$log2FoldChange) > 1]

# Save the gene IDs to a text file
write.table(DE_gene_ids, file = "differential_expression_gene_ids_SRP079684.txt", col.names = FALSE, quote = FALSE, row.names = FALSE)

# ====== top 20 genes =========
# Assuming your data frame is named 'deseq_df'
filtered_df <- deseq_df[!(deseq_df$threshold == FALSE | is.na(deseq_df$threshold)), ]

# Sort by absolute fold change
sorted_df <- filtered_df[order(abs(filtered_df$log2FoldChange), decreasing = TRUE), ]

# Pick the top 20 genes
top_20_genes <- sorted_df[1:20, ]

# Print or use the 'top_20_genes' data frame as needed
print(top_20_genes)

write.table(top_20_genes, file = "differential_expression_top_20_gene_ids_SRP079684.txt", col.names = FALSE, quote = FALSE, row.names = FALSE)





