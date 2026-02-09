library(mgsa)
library(dplyr)
library(tidyr)
library(stringr)
library(optparse)

# Parse arguments
option_list=list(
  make_option(c("-g", "--gene_file"), type="character", help="Expression data (csv or txt). Each row contains the gene name, optionally followed by tab plus expression level. No header!"),
  make_option(c("-p", "--preprocess_data"), type="character", help="Path to script to preprocess data (only if gene_file has more than 2 columns and with header)"),
  make_option(c("-a", "--annotation"), type="character", help="Annotation file path"),
  make_option(c("-o", "--output_filename"), type="character", help="Name for output file"),
  make_option(c("-c", "--lower_cutoff"),type="character",help="Lower cutoff for number of genes in term"), 
  make_option(c("-C", "--upper_cutoff"),type="character",help="Upper cutoff for number of genes in term"),
  make_option(c("-n", "--repeats"),type="character",help="Number of repeats")
)
parser = OptionParser(option_list=option_list)
args = parse_args(parser)
if(args$help) {
  print_help(parser)
  quit(save="no", status=0)
}

output_filename <- args$output_filename
annt <- read.csv(file = args$annotation, header = TRUE, sep= "\t")
if(!is.null(args$lower_cutoff)){
  lower_cutoff <- as.integer(args$lower_cutoff)
} else {
  lower_cutoff = 0
}
if(!is.null(args$upper_cutoff)){
  upper_cutoff <- as.integer(args$upper_cutoff)
} else {
  upper_cutoff = 100000
}
if(!is.null(args$repeats)){
  repeats = as.integer(args$repeats)
} else {
  repeats = 1
}
# Parsing annotations
colnames(annt) <- c("term","T")
T <- apply(annt['T'], 1, function(x) {
  str <- gsub("\\[|\\]|\'","",x)
  vec <- strsplit(str, ", ")[[1]]
  return(vec)
})
annt2 <- data.frame(term=rep(annt$term, sapply(T, length)), gene=unlist(T))
print(length(unique(unlist(annt2$term))))
# Implement upper and lower cutoff
if (lower_cutoff > 0 | upper_cutoff < 100000){
  annt2 <- annt2%>%
    group_by(term)%>%
    filter(n() >= lower_cutoff & n() <= upper_cutoff)%>%
    ungroup()
}
print(head(annt2))
print(length(unique(unlist(annt2$term))))
                
# Turn 2-column df of annt2 to list by term name and genes
GO <- split(annt2$gene, annt2$term)

# Preprocess gene_file if not 2 columns
if(is.null(args$preprocess_data)){
  gene_file <- read.table(file = args$gene_file, header = F, sep = "\t")
  if (ncol(gene_file) > 2){
    print("Not default gene file format. Please provide a script to preprocess the data.")
    quit(save="no", status=0)
  } 
} else {
  source(args$preprocess_data)
  gene_file <- preprocess_data(args$gene_file)
  gene_file$V2 <- as.numeric(gene_file$V2) # force numeric
  gene_file$V3 <- as.numeric(gene_file$V3) # force numeric
  gene_file$V4 <- as.numeric(gene_file$V4) # force numeric
}

# Only gene names as input for mgsa
o <- gene_file%>%
  filter(V3 < 0.05)

# genes <- gene_file$V1
# log2FC <- gene_file$V2
# p_values <- gene_file$V3
# stat <- gene_file$V4

# model p-values
result <- mgsa(o$V1, GO, alpha=seq(0.05,0.95, length.out=19), beta=seq(0.05,0.95, length.out=19))

df = as.data.frame((result@setsResults))
df$ID = rownames(df)
                        
write.table(df[, c("ID", "estimate")],
            file = output_filename,
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t")