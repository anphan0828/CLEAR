library(EnrichmentBrowser)
library(KEGGdzPathwaysGEO)
library(limma)

#' Process a single GEO microarray dataset for differential expression analysis
#'
#' @param files List of file paths to raw SummarizedExperiment files
#' @param analysis_func Function to use for analysis (run_deseq2_analysis)
#' @param output_dir Directory to save results (default: "results")
#' @return A list of result dataframes
perform_de_analysis <- function(eset, dataset_name, file_path) {
  cat("\n=== Analyzing dataset:", dataset_name, "===\n")
  
  # Get expression and phenotype data
  expr_data <- exprs(eset)
  pheno_data <- pData(eset)
  
  print(paste("Dataset dimensions:", nrow(expr_data), "features x", ncol(expr_data), "samples"))
  print(paste("Available phenotype variables:", colnames(pheno_data)))
  
  if(!"Group" %in% colnames(pheno_data)) {
    cat("No 'Group' variable found. Skipping", dataset_name, "\n")
    return(NULL)
  }
  
  # Check groups
  group_factor <- factor(pheno_data$Group)
  groups <- levels(group_factor)
  print(table(pheno_data$Group))
  
  if(length(groups) != 2) {
    cat("Dataset does not have exactly 2 groups. Skipping", dataset_name, "\n")
    return(NULL)
  }
  
  # Create design matrix
  design <- model.matrix(~0 + group_factor)
  colnames(design) <- levels(group_factor)
  
  # Fit linear model
  fit <- lmFit(expr_data, design)
  
  # Create contrast
  contrast_matrix <- makeContrasts(
    contrasts = paste(groups[2], groups[1], sep=" - "),
    levels = design
  )
  
  # Apply contrasts and compute statistics
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  # Extract results
  results <- topTable(fit2, number = Inf, sort.by = "none")
  results$ProbeID <- rownames(results)
  
  # Get annotation package based on platform to map probe IDs to gene symbols
  annotation_pkg <- annotation(eset)
  
  if(annotation_pkg == "hgu133a") {
    library(hgu133a.db)
    gene_symbols <- mapIds(hgu133a.db, keys = rownames(results), 
                          column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
  } else if(annotation_pkg == "hgu133plus2") {
    library(hgu133plus2.db)
    gene_symbols <- mapIds(hgu133plus2.db, keys = rownames(results),
                          column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
  } else {
    cat("Annotation package", annotation_pkg, "not supported. Using probe IDs only.\n")
    gene_symbols <- rep(NA, nrow(results))
  }
  
  # Add annotations
  results$GeneSymbol <- gene_symbols
  
  # Reorder columns
  results <- results[, c("ProbeID", "GeneSymbol", "logFC", "t", "P.Value", "adj.P.Val", "AveExpr", "B")]
  colnames(results) <- c("ProbeID", "GeneSymbol", "log2FoldChange", "stat", "pvalue", "padj", "baseMean", "B-stat")
  
  # Remove rows with NA gene symbols first
  results_with_genes <- results[!is.na(results$GeneSymbol), ]
  cat("After removing NA gene symbols:", nrow(results_with_genes), "probes\n")
  
  if(nrow(results_with_genes) > 0) {
    # Order by highest average expression across all samples (baseMean)
    results_with_genes <- results_with_genes[order(-results_with_genes$baseMean), ]
    # Remove duplicated gene symbols, keeping the first occurrence (lowest p-value)
    results_unique <- results_with_genes[!duplicated(results_with_genes$GeneSymbol), ]
    
    cat("After removing duplicate gene symbols:", nrow(results_unique), "genes\n")
    cat("Number of duplicated gene symbols removed:", nrow(results_with_genes) - nrow(results_unique), "\n")
    
    results <- results_unique
  } else {
    cat("Warning: No genes with valid symbols found!\n")
    results <- results_with_genes
  }
  # Save results
  write.table(results, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Summary
  n_up <- sum(results$pvalue < 0.05 & results$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(results$pvalue < 0.05 & results$log2FoldChange < 0, na.rm = TRUE)
  cat("Results saved to:", file_path, "\n")
  return(results)
}

### Main commands ####
data_names <- data(package="KEGGdzPathwaysGEO")$results[,"Item"]
print(data_names)
# data(GSE1297)

# Access data for target pathway and other details
# str(GSE1297)
# GSE1297@experimentData@other

results_list <- list()
all_datasets <- data_names
print("Available datasets:")
print(all_datasets)

# For de-duplicating by average expression
if(!dir.exists("real_data/geo_microarray_data")) {
  dir.create("real_data/geo_microarray_data")    
}

# Process each dataset
for (dataset_name in all_datasets) {
  tryCatch({
    # Load dataset
    data(list = dataset_name, package = "KEGGdzPathwaysGEO")
    eset <- get(dataset_name)
    
    # Perform DE analysis
    results <- perform_de_analysis(eset, dataset_name, file_path = file.path("real_data/geo_microarray_data", paste0(dataset_name, "_de_limma.tsv")))
    if(!is.null(results)) {
      results_list[[dataset_name]] <- results
    }   
  }, error = function(e) {
    cat("Error processing", dataset_name, ":", e$message, "\n")
  })
}

# Summary of all analyses (write to file)
summary_file <- "real_data/geo_microarray_data_summary.txt"

cat("=== SUMMARY OF ALL ANALYSES ===\n", file = summary_file, sep = "\n")
cat("Successfully processed datasets:", length(results_list), "out of", length(all_datasets), "\n", 
    file = summary_file, append = TRUE)
for(dataset_name in names(results_list)) {
  results <- results_list[[dataset_name]]
  n_up <- sum(results$pvalue < 0.05 & results$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(results$pvalue < 0.05 & results$log2FoldChange < 0, na.rm = TRUE)
  cat(sprintf("%-12s: %5d up, %5d down, %5d total genes\n", 
              dataset_name, n_up, n_down, nrow(results)), file = summary_file, append = TRUE)
}
