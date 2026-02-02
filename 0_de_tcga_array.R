library(BiocManager)
library(GSEABenchmarkeR)
library(reticulate)

#' Extract metadata from TCGA raw RData files
#'
#' @param tcga_raw_dir Directory containing TCGA raw RData files  
#' @param output_file Path to save the metadata TSV file
#' @return A data frame containing TCGA dataset metadata
extract_tcga_metadata <- function(tcga_raw_dir = "tcga_raw", output_file = "tcga_datasets_metadata.tsv") {
  # TCGA study abbreviations lookup
  tcga_studies <- list(
    'ACC' = 'Adrenocortical carcinoma',
    'BLCA' = 'Bladder Urothelial Carcinoma',
    'BRCA' = 'Breast invasive carcinoma',
    'CESC' = 'Cervical squamous cell carcinoma and endocervical adenocarcinoma',
    'CHOL' = 'Cholangiocarcinoma',
    'COAD' = 'Colon adenocarcinoma',
    'ESCA' = 'Esophageal carcinoma',
    'GBM' = 'Glioblastoma multiforme',
    'HNSC' = 'Head and Neck squamous cell carcinoma',
    'KICH' = 'Kidney Chromophobe',
    'KIRC' = 'Kidney renal clear cell carcinoma',
    'KIRP' = 'Kidney renal papillary cell carcinoma',
    'LAML' = 'Acute Myeloid Leukemia',
    'LGG' = 'Brain Lower Grade Glioma',
    'LIHC' = 'Liver hepatocellular carcinoma',
    'LUAD' = 'Lung adenocarcinoma',
    'LUSC' = 'Lung squamous cell carcinoma',
    'MESO' = 'Mesothelioma',
    'OV' = 'Ovarian serous cystadenocarcinoma',
    'PAAD' = 'Pancreatic adenocarcinoma',
    'PCPG' = 'Pheochromocytoma and Paraganglioma',
    'PRAD' = 'Prostate adenocarcinoma',
    'READ' = 'Rectum adenocarcinoma',
    'SARC' = 'Sarcoma',
    'SKCM' = 'Skin Cutaneous Melanoma',
    'STAD' = 'Stomach adenocarcinoma',
    'TGCT' = 'Testicular Germ Cell Tumors',
    'THCA' = 'Thyroid carcinoma',
    'THYM' = 'Thymoma',
    'UCEC' = 'Uterine Corpus Endometrial Carcinoma',
    'UCS' = 'Uterine Carcinosarcoma',
    'UVM' = 'Uveal Melanoma'
  )
  
  # Get list of RData files
  files <- list.files(tcga_raw_dir, pattern = "*.Rdata", full.names = TRUE)
  
  # Initialize results list
  metadata_list <- list()
  
  cat("Extracting metadata from", length(files), "TCGA raw files...\n")
  
  for(i in 1:length(files)) {
    rdata_file <- files[i]
    dataset_code <- gsub(".Rdata", "", basename(rdata_file))
    
    cat("Processing:", dataset_code, "\n")
    
    tryCatch({
      # Load RData file
      load(rdata_file)
      
      # Extract metadata from SummarizedExperiment object
      se <- summarized_exp  # assuming the object is named summarized_exp
      
      # Extract sample information
      col_data <- colData(se)
      group_info <- col_data$GROUP
      
      # Count samples by group
      group_counts <- table(group_info)
      reference_size <- as.numeric(group_counts["0"])  # assuming 0 is reference
      test_size <- as.numeric(group_counts["1"])       # assuming 1 is test
      
      # Handle NAs
      if(is.na(reference_size)) reference_size <- 0
      if(is.na(test_size)) test_size <- 0
      
      # Create metadata entry
      metadata_entry <- data.frame(
        dataset = dataset_code,
        cancer_name = ifelse(dataset_code %in% names(tcga_studies), 
                           tcga_studies[[dataset_code]], 
                           "Unknown"),
        reference_group = "0",
        test_group = "1",
        reference_size = reference_size,
        test_size = test_size,
        total_samples = ncol(se),
        total_genes = nrow(se),
        stringsAsFactors = FALSE
      )
      
      metadata_list[[i]] <- metadata_entry
      
    }, error = function(e) {
      cat("Error processing", dataset_code, ":", e$message, "\n")
      
      # Add placeholder for failed extractions
      metadata_entry <- data.frame(
        dataset = dataset_code,
        cancer_name = ifelse(dataset_code %in% names(tcga_studies), 
                           tcga_studies[[dataset_code]], 
                           "Unknown"),
        reference_group = "0",
        test_group = "1",
        reference_size = NA,
        test_size = NA,
        total_samples = NA,
        total_genes = NA,
        extraction_error = TRUE,
        stringsAsFactors = FALSE
      )
      
      metadata_list[[i]] <- metadata_entry
    })
  }
  
  # Combine all metadata
  metadata_df <- do.call(rbind, metadata_list)
  
  # Save to file
  write.table(metadata_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("TCGA metadata saved to:", output_file, "\n")
  
  # Print summary
  cat("\n=== TCGA Metadata Summary ===\n")
  print(metadata_df)
  
  return(metadata_df)
}

#' Load the MalaCards phenotype rankings from RDS and save as pickle
#'
#' @param rds_path Path to the input RDS file (downloaded from GSEABenchmarkeR github)
#' @param pickle_path Path to save the output pickle file
#' @return None
rds_to_pickle <- function(rds_path, pickle_path) {
  # Import required Python modules through reticulate
  pickle <- import("pickle")
  pd <- import("pandas")

  # Read the RDS file
  r_data <- readRDS(rds_path)

  # Convert R data frames to Python pandas DataFrames
  if (is.data.frame(r_data)) {
    # If single dataframe, convert directly
    py_data <- r_to_py(r_data)
  } else if (is.list(r_data)) {
    # If list of dataframes, convert each element
    py_data <- list()
    for (name in names(r_data)) {
      if (class(r_data[[name]])[1] == "DFrame") {
        # Convert DFrame to regular R data.frame first
        go_ids <- rownames(r_data[[name]])
	df <- data.frame(
	  GO.ID = go_ids,
          TITLE = as.character(r_data[[name]]$TITLE),
          REL.SCORE = as.numeric(r_data[[name]]$REL.SCORE),
          MATCHED.GENES = as.integer(r_data[[name]]$MATCHED.GENES),
          TOTAL.GENES = as.integer(r_data[[name]]$TOTAL.GENES)
        )
        py_data[[name]] <- r_to_py(df)
      } else {
        warning(paste("Skipping non-dataframe element:", name))
      }
    }
  } else {
    stop("Input must be either a dataframe or a list of dataframes")
  }

  # Save as pickle file using Python's built-in open function
  py_save_object(py_data,filename=pickle_path,pickle="pickle")
}



#' Run DESeq2 analysis on RNA-seq data
#'
#' @param se The experiment summary object loaded from RData file
#' @param return_all_objects Logical indicating whether to return all objects (default: FALSE)
#' @return A data frame containing differential expression results or a list of objects
library(DESeq2)
run_deseq2_analysis <- function(se, return_all_objects = FALSE) {
  
  isSE <- is(se, "SummarizedExperiment")
  if(isSE) 
  { 
    expr <- assay(se)
    grp <- colData(se)[, "GROUP"]
    blk <- colData(se)[, "BLOCK"]
  }
  group <- factor(grp)
  paired <- !is.null(blk)
  f <- "~"
  if(paired) 
  {
      block <- factor(blk)
      f <- paste0(f, "block + ") 
  }
  f <- formula(paste0(f, "group"))
    
  # de.tbl <- .deseq(expr, group, paired, block, f, stat.only)
  colData2 <- data.frame(group=group)
  if(paired) colData2$block <- block
  suppressMessages({
      dds <- DESeq2::DESeqDataSetFromMatrix(
          countData=expr, colData=colData2, design=f)
      # Add pre-filtering step
      smallestGroupSize <- 3
      keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
      dds <- dds[keep,]
      # Do DESeq2 analysis
      dds <- DESeq2::DESeq(dds)
  })
  res <- DESeq2::results(dds, pAdjustMethod="BH")

  # Convert to data frame and add gene IDs
  res_df <- as.data.frame(res)
  res_df <- res_df[order(res_df$pvalue), ]
  
  # Add gene IDs as a column if they're only in rownames; for TCGA data it is entrez_gene_id
  if(!("entrez_gene_id" %in% colnames(res_df))) {
    res_df$entrez_gene_id <- rownames(res_df)
    res_df <- map_gene_name(res_df)
  }
  
  if(return_all_objects) {
    return(list(
      dds = dds,
      results_obj = res,
      results = res_df
    ))
  } else {
    res_df <- res_df[order(res_df$pvalue), c("entrez_gene_id", "external_gene_name", "baseMean","log2FoldChange","lfcSE","stat","pvalue", "padj")]
    return(res_df)
  }
}

map_gene_name <- function(df){
    library(org.Hs.eg.db)
    symbols <- mapIds(org.Hs.eg.db,
                        keys = df$entrez_gene_id,
                        column = "SYMBOL",
                        keytype = "ENTREZID",
                        multiVals = "first")
    df$external_gene_name <- symbols
    return(df)
}

#' Process raw SummarizedExperiment file and run analysis
#'
#' @param files List of file paths to raw SummarizedExperiment files
#' @param analysis_func Function to use for analysis (run_deseq2_analysis)
#' @param output_dir Directory to save results (default: "results")
#' @return A list of result dataframes
process_raw_file <- function(files, analysis_func, output_dir = "results",id=0) {
  
  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  n = length(files)
  # Process each contrast
  if (id==0) {
    print(sprintf("Processing all %d datasets",n))
  }
  else {
    print(sprintf("Processing contrasts number %d",id))
  }
  for(i in 1:n) {
    if (id == 0 | id == i){
      # File name for RData file
      rdata_file <- files[i]
      
      # Skip if file doesn't exist
      if(!file.exists(rdata_file)) {
        warning(paste("File not found:", rdata_file))
        next
      }
      
      # Load RData file
      load(rdata_file)
      
      # Run analysis
      message(paste("Processing", files[i]))
      output_file <- file.path(output_dir, paste0(gsub(".Rdata", ".tsv", basename(rdata_file))))
      if (!file.exists(output_file)){
        results <- analysis_func(summarized_exp)
        write.table(results, output_file, row.names = FALSE,sep="\t",quote=FALSE)
      } else {
        message(paste(basename(rdata_file),"already processed."))
      }
    }
  }
}

### Main commands ####
root_folder = "real_data/tcga_raw"

# Retrieve RData files of TCGA RNA-seq datasets from GSEABenchmarkeR::loadEData ####
tcga <- loadEData("tcga", nr.datasets=15)
output_dir <- root_folder
if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

for (i in 1:length(tcga)){
  dataset=names(tcga)[i]
  summarized_exp = tcga[[i]]
  save(summarized_exp,file=paste0(output_dir,"/",dataset,".Rdata"))
}

# Extract metadata from TCGA raw files first
cat("=== Extracting TCGA Metadata ===\n")
tcga_metadata <- extract_tcga_metadata(tcga_raw_dir = root_folder, 
                                      output_file = "real_data/tcga_datasets_metadata.tsv")

# Get MalaCards phenotype rankings and save as pickle
mala_go_path = "GO_BP.rds"
pickle_path = sub(".rds", ".cp", mala_go_path)
if (!file.exists(pickle_path)){
  rds_to_pickle(mala_go_path,pickle_path)
}

# List files in root_folder
files <- list.files(root_folder, pattern = "*.Rdata", full.names = TRUE)

# Set up BiocParallel for parallel processing
BiocParallel::register(BiocParallel::MulticoreParam(workers=16))
bp.param <- BiocParallel::registered()[[1]]
BiocParallel::bpprogressbar(bp.param) <- TRUE

# Run as SLURM array job
args = commandArgs(trailingOnly=TRUE)
SLURM_ARRAY_TASK_ID = as.integer(args[1])
ARRAY_SIZE = as.integer(args[2])
selected_list = c()
for (id in 1:length(files)){
  if (id %% ARRAY_SIZE == SLURM_ARRAY_TASK_ID){
    process_raw_file(files, run_deseq2_analysis, output_dir = "real_data/tcga_data",id)
  }
}