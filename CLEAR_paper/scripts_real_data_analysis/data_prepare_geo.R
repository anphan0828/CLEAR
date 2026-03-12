library(dplyr)

preprocess_data <- function(gene_file){
    data <- read.table(gene_file, header = TRUE, sep = "\t")
    data <- data %>%
        select(c("GeneSymbol", "log2FoldChange", "pvalue", "stat"))
    data$log2FoldChange <- as.numeric(data$log2FoldChange)
    data$pvalue <- as.numeric(data$pvalue)
    data$stat <- as.numeric(data$stat)
    colnames(data) <- c("V1", "V2", "V3", "V4")
    data <- data %>%
        filter(!is.na(V1) & !is.na(V2) & !is.na(V3) & !is.na(V4))
    return(data)
}