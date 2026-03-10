suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(stringr)
  library(readr)
})

# =========================================================
# USER INPUT
# =========================================================
dataset_num <- 389
model_name  <- "beta"
top_n       <- 20


res_dir  <- "Your path to /res"
base_out_dir <- "Your path to /inspect_chains"
dataset_dir <- file.path(
  base_out_dir,
  paste0(model_name, "_gen", dataset_num)
)

dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)


dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)

# Optional burn-in to drop from trace plots
burnin <- 0

# =========================================================
# CONFIG: candidate variable names inside each RDS
# edit here only if your saved names differ
# =========================================================

posterior_trace_candidates <- c(
  "log_likelihoods"
)

para1_candidates <- c(
  "p1_trace"
)

para2_candidates <- c(
  "p2_trace"
)

gene_prob_candidates <- c(
  "on_frequency"
)

# =========================================================
# HELPERS
# =========================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

find_first_existing <- function(obj, candidates) {
  nm <- candidates[candidates %in% names(obj)]
  if (length(nm) == 0) return(NULL)
  nm[1]
}

as_clean_numeric <- function(x, nm = "object", burnin = 0) {
  if (is.null(x)) return(NULL)
  
  if (is.data.frame(x)) {
    if (ncol(x) == 1) {
      x <- x[[1]]
    } else {
      stop(nm, " is a data.frame with ", ncol(x), " columns; expected one numeric vector.")
    }
  }
  
  if (is.matrix(x)) {
    if (ncol(x) == 1) {
      x <- x[, 1]
    } else if (nrow(x) == 1) {
      x <- x[1, ]
    } else {
      stop(nm, " is a matrix with dimension ", paste(dim(x), collapse = "x"),
           "; expected one chain vector.")
    }
  }
  
  if (is.list(x) && !is.atomic(x)) {
    stop(nm, " is a list; expected one numeric vector.")
  }
  
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if (burnin > 0) {
    if (length(x) <= burnin) stop(nm, " length <= burnin.")
    x <- x[(burnin + 1):length(x)]
  }
  
  x
}

as_named_numeric <- function(x, nm = "object") {
  if (is.null(x)) return(NULL)
  
  if (is.data.frame(x)) {
    if (ncol(x) == 1) {
      vals <- as.numeric(x[[1]])
      names(vals) <- rownames(x)
      x <- vals
    } else {
      stop(nm, " is a data.frame with ", ncol(x), " columns; expected named numeric vector.")
    }
  }
  
  if (is.matrix(x)) {
    if (ncol(x) == 1) {
      vals <- as.numeric(x[, 1])
      names(vals) <- rownames(x)
      x <- vals
    } else if (nrow(x) == 1) {
      vals <- as.numeric(x[1, ])
      names(vals) <- colnames(x)
      x <- vals
    } else {
      stop(nm, " is a matrix with dimension ", paste(dim(x), collapse = "x"),
           "; expected 1-column or 1-row object.")
    }
  }
  
  x <- as.numeric(x) %>% setNames(names(x))
  if (is.null(names(x))) stop(nm, " has no names.")
  x <- x[is.finite(x)]
  x
}

extract_one_chain <- function(path, burnin = 0) {
  obj <- readRDS(path)
  
  posterior_nm <- find_first_existing(obj, posterior_trace_candidates)
  para1_nm     <- find_first_existing(obj, para1_candidates)
  para2_nm     <- find_first_existing(obj, para2_candidates)
  geneprob_nm  <- find_first_existing(obj, gene_prob_candidates)
  
  posterior_trace <- if (!is.null(posterior_nm)) {
    as_clean_numeric(obj[[posterior_nm]], posterior_nm, burnin = burnin)
  } else NULL
  
  para1_trace <- if (!is.null(para1_nm)) {
    as_clean_numeric(obj[[para1_nm]], para1_nm, burnin = burnin)
  } else NULL
  
  para2_trace <- if (!is.null(para2_nm)) {
    as_clean_numeric(obj[[para2_nm]], para2_nm, burnin = burnin)
  } else NULL
  
  gene_prob <- if (!is.null(geneprob_nm)) {
    as_named_numeric(obj[[geneprob_nm]], geneprob_nm)
  } else NULL
  
  list(
    path = path,
    posterior_name = posterior_nm,
    para1_name = para1_nm,
    para2_name = para2_nm,
    geneprob_name = geneprob_nm,
    posterior_trace = posterior_trace,
    para1_trace = para1_trace,
    para2_trace = para2_trace,
    gene_prob = gene_prob
  )
}

build_trace_df <- function(chain_list, field) {
  rows <- lapply(seq_along(chain_list), function(i) {
    x <- chain_list[[i]][[field]]
    if (is.null(x)) return(NULL)
    tibble(
      chain = paste0("rep", i),
      iter = seq_along(x),
      value = x
    )
  })
  bind_rows(rows)
}

plot_trace <- function(df, title, ylab) {
  ggplot(df, aes(x = iter, y = value, color = chain)) +
    geom_line(linewidth = 0.35, alpha = 0.9) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = title,
      x = "Iteration",
      y = ylab,
      color = "Chain"
    )
}

# =========================================================
# FIND FILES
# =========================================================

pattern <- paste0("^result_", model_name, "_gen", dataset_num, "_rep[1-9][0-9]*\\.rds$")
files <- list.files(res_dir, pattern = pattern, full.names = TRUE)

if (length(files) == 0) {
  stop("No files found for model=", model_name, ", dataset gen", dataset_num)
}

meta <- tibble(
  file = files,
  file_base = basename(files),
  rep = as.integer(str_extract(file_base, "rep\\d+") %>% str_remove("rep"))
) %>%
  arrange(rep)

if (length(unique(meta$rep)) != nrow(meta)) {
  warning("Duplicate rep IDs found.")
}

message("Found files:")
print(meta)

# =========================================================
# LOAD 5 CHAINS
# =========================================================

chains <- lapply(meta$file, extract_one_chain, burnin = burnin)

# show detected variable names
var_detect_tbl <- tibble(
  rep = meta$rep,
  posterior_trace_var = sapply(chains, `[[`, "posterior_name"),
  para1_trace_var     = sapply(chains, `[[`, "para1_name"),
  para2_trace_var     = sapply(chains, `[[`, "para2_name"),
  gene_prob_var       = sapply(chains, `[[`, "geneprob_name")
)

message("Detected variable names:")
print(var_detect_tbl)

# =========================================================
# TRACE PLOTS
# =========================================================

posterior_df <- build_trace_df(chains, "posterior_trace")
para1_df     <- build_trace_df(chains, "para1_trace")
para2_df     <- build_trace_df(chains, "para2_trace")

if (nrow(posterior_df) > 0) {
  p_post <- plot_trace(
    posterior_df,
    title = paste0("Posterior trace: ", model_name, ", gen", dataset_num),
    ylab = "Posterior"
  )
  ggsave(
    file.path(dataset_dir, paste0("trace_posterior_", model_name, "_gen", dataset_num, ".png")),
    p_post, width = 5, height = 4, units = "in", bg = "white"
  )
} else {
  message("No posterior trace variable found.")
}

if (nrow(para1_df) > 0) {
  p_para1 <- plot_trace(
    para1_df,
    title = paste0("Para1 trace: ", model_name, ", gen", dataset_num),
    ylab = "para1"
  )
  ggsave(
    file.path(dataset_dir, paste0("trace_para1_", model_name, "_gen", dataset_num, ".png")),
    p_para1, width = 5, height = 4, units = "in", bg = "white"
  )
} else {
  message("No para1 trace variable found.")
}

if (nrow(para2_df) > 0) {
  p_para2 <- plot_trace(
    para2_df,
    title = paste0("Para2 trace: ", model_name, ", gen", dataset_num),
    ylab = "para2"
  )
  ggsave(
    file.path(dataset_dir, paste0("trace_para2_", model_name, "_gen", dataset_num, ".png")),
    p_para2, width = 5, height = 4, units = "in", bg = "white"
  )
} else {
  message("No para2 trace variable found.")
}

# =========================================================
# TOP N GENE SET TABLE BY AVG POSTERIOR PROB
# =========================================================

gene_prob_list <- lapply(chains, `[[`, "gene_prob")

if (all(vapply(gene_prob_list, is.null, logical(1)))) {
  message("No gene posterior probability object found, so top gene table was not created.")
} else {
  common_genes <- Reduce(intersect, lapply(gene_prob_list, names))
  
  if (length(common_genes) == 0) {
    stop("No overlapping named gene sets across chains for posterior probabilities.")
  }
  
  gene_prob_mat <- sapply(gene_prob_list, function(x) x[common_genes])
  colnames(gene_prob_mat) <- paste0("rep", meta$rep)
  
  avg_prob <- rowMeans(gene_prob_mat, na.rm = TRUE)
  
  top_genes <- names(sort(avg_prob, decreasing = TRUE))[seq_len(min(top_n, length(avg_prob)))]
  
  top_tbl <- data.frame(
    gene_set = top_genes,
    avg_posterior_prob = avg_prob[top_genes],
    gene_prob_mat[top_genes, , drop = FALSE],
    row.names = NULL,
    check.names = FALSE
  )
  
  top_tbl <- top_tbl %>%
    arrange(desc(avg_posterior_prob))
  
  print(top_tbl)
  
  write_csv(
    top_tbl,
    file.path(dataset_dir, paste0("top_gene_probs_", model_name, "_gen", dataset_num, ".csv"))
  )
  
  saveRDS(
    top_tbl,
    file.path(dataset_dir, paste0("top_gene_probs_", model_name, "_gen", dataset_num, ".rds"))
  )
}

message("Done.")
