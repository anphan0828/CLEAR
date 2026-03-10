# Generate in-silico datasets and save metadata

rm(list = ls())

## =========================
## Paths
## =========================
go_path <- "Your path to GO_ecoli.rds"
data_dir <- "Your path to save generated datasets"
result_dir <- "Your path to result folder"

## =========================
## Helper functions
## =========================
is_fully_contained <- function(vec1, vec2) {
  all(vec1 %in% vec2)
}

get_not_hierarchical_terms <- function(GO, n) {
  unrelated_terms <- list()
  names_GO <- sample(names(GO))
  
  for (term_name in names_GO) {
    term <- GO[[term_name]]
    
    if (!any(sapply(unrelated_terms, is_fully_contained, term)) &&
        !any(sapply(unrelated_terms, function(x) is_fully_contained(term, x)))) {
      unrelated_terms[[term_name]] <- term
    }
    
    if (length(unrelated_terms) == n) {
      break
    }
  }
  
  unrelated_terms
}

generate_data <- function(GO, n_on_terms, alt_para, up_ratio = 0.5, df = 50) {
  repeat {
    genes_all <- unique(unlist(GO))
    term_on <- get_not_hierarchical_terms(GO, n_on_terms)
    
    gene_on <- sample(unique(unlist(term_on)))
    gene_off <- sample(setdiff(genes_all, gene_on))
    
    n_up <- ceiling(length(gene_on) * up_ratio)
    n_down <- length(gene_on) - n_up
    
    gene_up <- gene_on[seq_len(n_up)]
    gene_down <- gene_on[(n_up + 1):length(gene_on)]
    
    t_up <- rt(length(gene_up), df = df, ncp = alt_para)
    t_down <- rt(length(gene_down), df = df, ncp = -alt_para)
    t_off <- rt(length(gene_off), df = df, ncp = 0)
    
    genes <- c(gene_up, gene_down, gene_off)
    t_stat <- c(t_up, t_down, t_off)
    p_values <- 2 * pt(-abs(t_stat), df = df)
    
    valid <- all(is.finite(t_stat)) &&
      all(!is.na(t_stat)) &&
      all(is.finite(p_values)) &&
      all(!is.na(p_values)) &&
      all(p_values >= 0 & p_values <= 1)
    
    if (valid) {
      return(list(
        term_on = names(term_on),
        genes = genes,
        t_stat = t_stat,
        p_values = p_values,
        alt_para = alt_para
      ))
    }
    
    cat("Regenerating due to invalid t statistics or p-values...\n")
  }
}

## =========================
## Design settings
## =========================
ncp_vals <- c(1, 2, 3, 5, 10)
nterm_vals <- c(5, 10, 20, 50)

design_grid <- expand.grid(
  ncp = ncp_vals,
  n_term_on = nterm_vals,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

design_list <- lapply(seq_len(nrow(design_grid)), function(i) {
  list(
    control_para = c("ncp", "n_term_on"),
    control_para_thres = list(
      c(design_grid$ncp[i], design_grid$ncp[i]),
      c(design_grid$n_term_on[i], design_grid$n_term_on[i])
    ),
    plot_main = paste0(
      "ncp=", design_grid$ncp[i],
      ", ", design_grid$n_term_on[i], " terms on"
    )
  )
})

n_datasets <- 100

## =========================
## Load GO and create folders
## =========================
GO_ecoli <- readRDS(go_path)

if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

if (!dir.exists(result_dir)) {
  dir.create(result_dir, recursive = TRUE)
}

## =========================
## Generate datasets
## =========================
gen_id <- 1
design_record <- integer()

for (i in seq_along(design_list)) {
  design <- design_list[[i]]
  para_ncp <- design$control_para_thres[[1]][1]
  n_term_on_val <- design$control_para_thres[[2]][1]
  
  cat("Generating data for design:", design$plot_main, "\n")
  
  for (j in seq_len(n_datasets)) {
    gen <- generate_data(
      GO = GO_ecoli,
      n_on_terms = n_term_on_val,
      alt_para = para_ncp
    )
    
    saveRDS(gen, file = file.path(data_dir, paste0("gen", gen_id, ".rds")))
    design_record[gen_id] <- i
    gen_id <- gen_id + 1
  }
}

## =========================
## Build metadata
## =========================
terms_on <- list()
para_ncp <- numeric()
n_term_on <- numeric()

n_total <- n_datasets * length(design_list)

for (i in seq_len(n_total)) {
  gen <- readRDS(file.path(data_dir, paste0("gen", i, ".rds")))
  
  terms_on[[i]] <- gen$term_on
  para_ncp[i] <- gen$alt_para[1]
  n_term_on[i] <- length(gen$term_on)
}

meta <- data.frame(
  para_ncp = para_ncp,
  n_term_on = n_term_on
)

saveRDS(
  list(meta = meta, terms_on = terms_on),
  file = file.path(result_dir, "meta.rds")
)

## =========================
## Create result subfolders
## =========================
output_dirs <- c(
  "MGSA_result_fpr.05",
  "MGSA_result_fpr.1",
  "MGSA_result_fpr.2",
  "CLEAR3_s_tnormal_result",
  "CLEAR3_tnormal_result",
  "CLEAR3_gamma_result",
  "CLEAR3_beta_result",
  "GSEA_result",
  "hyper_result"
)

for (dir_name in output_dirs) {
  dir_to_create <- file.path(result_dir, dir_name)
  
  if (!dir.exists(dir_to_create)) {
    dir.create(dir_to_create, recursive = TRUE)
    cat(sprintf("Directory '%s' created.\n", dir_name))
  } else {
    cat(sprintf("Directory '%s' already exists.\n", dir_name))
  }
}