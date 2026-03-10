suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(readr)
})

# =========================================================
# USER SETTINGS
# =========================================================

res_dir  <- "Your path to /res"
conv_dir <- "Your path to /conv_res"
indiv_dir <- file.path(conv_dir, "indiv")
plot_dir <- file.path(dirname(res_dir), "plots")

# Optional GO filter:
# If NULL, keep all terms found in files.
# If you want to restrict to a known GO universe, set GO_names <- names(readRDS("..."))
GO_names <- NULL

# If you want to drop burn-in for continuous traces, set this to an integer.
# Example: burnin <- 5000
burnin <- 500

# Thresholds used in summaries / plot annotations
rhat_threshold <- 1.05
weiss_alpha <- 0.05
bill_alpha  <- 0.05

dir.create(conv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(indiv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# HELPERS
# =========================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

.clip <- function(x, lo = -0.999, hi = 0.999) {
  pmin(pmax(x, lo), hi)
}

.arrow_order <- c("1->0", "1->1", "0->0", "0->1")

safe_numeric_vector <- function(x, nm = "trace") {
  if (is.null(x)) return(NULL)
  
  if (is.data.frame(x)) {
    if (ncol(x) != 1) stop(nm, " is a data.frame with ", ncol(x), " columns; expected one numeric chain vector.")
    x <- x[[1]]
  }
  if (is.matrix(x)) {
    if (ncol(x) != 1) stop(nm, " is a matrix with ", ncol(x), " columns; expected one numeric chain vector.")
    x <- as.numeric(x[, 1])
  }
  if (is.list(x) && !is.atomic(x)) {
    stop(nm, " is a list; expected one numeric chain vector.")
  }
  
  x <- as.numeric(x)
  
  # drop burn-in by iteration index first
  if (burnin > 0) {
    if (length(x) <= burnin) stop(nm, " raw length <= burnin.")
    x <- x[(burnin + 1):length(x)]
  }
  
  # now check whether post-burnin part is usable
  if (any(!is.finite(x))) {
    warning(nm, " contains non-finite values after burn-in; removing them.")
    x <- x[is.finite(x)]
  }
  
  if (length(x) == 0) stop(nm, " has no finite values after burn-in.")
  
  x
}

rename_trans_cols <- function(df) {
  stopifnot(is.matrix(df) || is.data.frame(df))
  df <- as.data.frame(df)
  cn <- colnames(df)
  
  if (!all(.arrow_order %in% cn)) {
    stop(
      "Transition table must contain columns: ",
      paste(.arrow_order, collapse = ", "),
      ". Found: ", paste(cn, collapse = ", ")
    )
  }
  
  df[, .arrow_order, drop = FALSE]
}

normalize_transition_object <- function(obj) {
  if (!is.null(obj$transition_counts)) {
    trans_counts <- rename_trans_cols(as.data.frame(obj$transition_counts))
    if (is.null(rownames(trans_counts))) {
      stop("transition_counts must have rownames = gene sets.")
    }
    rs <- rowSums(trans_counts, na.rm = TRUE)
    trans_freq <- trans_counts
    nz <- rs > 0
    if (any(nz)) {
      trans_freq[nz, ] <- sweep(trans_counts[nz, , drop = FALSE], 1, rs[nz], "/")
    }
    list(type = "counts", transition_freq = trans_freq, transition_counts = trans_counts)
  } else if (!is.null(obj$transition_freq)) {
    trans_freq <- rename_trans_cols(as.data.frame(obj$transition_freq))
    if (is.null(rownames(trans_freq))) {
      stop("transition_freq must have rownames = gene sets.")
    }
    list(type = "freq", transition_freq = trans_freq, transition_counts = NULL)
  } else {
    stop("RDS must contain either transition_counts or transition_freq.")
  }
}

infer_n_iter <- function(obj, trans_info) {
  if (!is.null(trans_info$transition_counts) && nrow(trans_info$transition_counts) > 0) {
    rs <- rowSums(trans_info$transition_counts, na.rm = TRUE)
    rs <- rs[rs > 0]
    if (length(rs) > 0) {
      T_est <- rs + 1L
      return(as.integer(names(sort(table(T_est), decreasing = TRUE)[1])))
    }
  }
  
  p1 <- obj$p1_trace %||% NULL
  if (!is.null(p1)) {
    if (is.data.frame(p1)) {
      if (ncol(p1) != 1) return(NA_integer_)
      p1 <- p1[[1]]
    }
    if (is.matrix(p1)) {
      if (ncol(p1) != 1) return(NA_integer_)
      p1 <- p1[, 1]
    }
    if (is.list(p1) && !is.atomic(p1)) return(NA_integer_)
    
    p1 <- as.numeric(p1)
    return(length(p1))
  }
  
  NA_integer_
}

# =========================================================
# CONVERGENCE METRICS
# =========================================================

rhat_basic <- function(chains_list) {
  stopifnot(length(chains_list) >= 2)
  
  n <- min(vapply(chains_list, length, integer(1)))
  
  if (n < 10) warning("Very short chains for Rhat; results may be unstable.")
  
  X <- sapply(chains_list, function(v) v[seq_len(n)])
  chain_means <- colMeans(X)
  
  B <- n * var(chain_means)
  W <- mean(apply(X, 2, var))
  
  if (!is.finite(W) || W <= 0) return(NA_real_)
  
  var_hat <- ((n - 1) / n) * W + (B / n)
  sqrt(var_hat / W)
}

ess_continuous <- function(x, max_lag = NULL) {
  x <- as.numeric(x)
  n <- length(x)
  
  if (n < 10) return(NA_real_)
  
  x <- x - mean(x)
  
  ac <- acf(
    x,
    plot = FALSE,
    lag.max = if (is.null(max_lag)) floor(n / 2) else max_lag,
    type = "correlation"
  )$acf[-1]
  
  pos <- ac > 0
  pos[is.na(pos)] <- FALSE
  
  if (!any(pos)) return(n)
  
  k <- which(!pos)
  if (length(k)) ac <- ac[seq_len(min(k) - 1)]
  
  tau <- 1 + 2 * sum(ac)
  ess <- n / tau
  
  max(1, min(ess, n))
}

ess_continuous_chains <- function(chains_list) {
  sapply(chains_list, ess_continuous)
}

ess_binary_from_phi <- function(n_iter, phi) {
  N <- pmax(as.numeric(n_iter), 2)
  ph <- .clip(as.numeric(phi))
  ess <- N * (1 - ph) / (1 + ph)
  pmax(1, ess)
}

weiss_between_chains <- function(on_freq_vec, trans_freq_mat, n_iter) {
  S <- length(on_freq_vec)
  stopifnot(nrow(trans_freq_mat) == S)
  
  if (length(n_iter) == 1L) n_iter <- rep(n_iter, S)
  
  p01 <- trans_freq_mat[, "0->1"]
  p10 <- trans_freq_mat[, "1->0"]
  phi <- .clip(1 - p01 - p10)
  
  c_hat_i <- (1 + phi) / (1 - phi)
  w <- n_iter / sum(n_iter)
  c_hat <- sum(w * c_hat_i, na.rm = TRUE)
  
  p_i    <- on_freq_vec
  p_pool <- sum(n_iter * p_i) / sum(n_iter)
  
  X2 <- sum(n_iter * (p_i - p_pool)^2 / p_pool) +
    sum(n_iter * ((1 - p_i) - (1 - p_pool))^2 / (1 - p_pool))
  
  X2_adj <- X2 / c_hat
  df <- (2 - 1) * (S - 1)
  pval <- pchisq(X2_adj, df = df, lower.tail = FALSE)
  
  list(stat = X2_adj, df = df, pval = pval, phi_per_chain = phi, c_hat = c_hat)
}

billingsley_between_chains <- function(trans_freq_mat, n_iter) {
  stopifnot(all(colnames(trans_freq_mat) %in% c("1->0", "1->1", "0->0", "0->1")))
  
  S <- nrow(trans_freq_mat)
  if (length(n_iter) == 1L) n_iter <- rep(n_iter, S)
  
  T_i <- pmax(n_iter - 1L, 1L)
  counts <- trans_freq_mat * T_i
  
  f10 <- counts[, "1->0"]
  f11 <- counts[, "1->1"]
  f00 <- counts[, "0->0"]
  f01 <- counts[, "0->1"]
  
  f1  <- f10 + f11
  f0  <- f00 + f01
  
  p10_i <- ifelse(f1 > 0, f10 / f1, NA_real_)
  p11_i <- ifelse(f1 > 0, f11 / f1, NA_real_)
  p00_i <- ifelse(f0 > 0, f00 / f0, NA_real_)
  p01_i <- ifelse(f0 > 0, f01 / f0, NA_real_)
  
  p10 <- sum(f10, na.rm = TRUE) / sum(f1, na.rm = TRUE)
  p11 <- sum(f11, na.rm = TRUE) / sum(f1, na.rm = TRUE)
  p00 <- sum(f00, na.rm = TRUE) / sum(f0, na.rm = TRUE)
  p01 <- sum(f01, na.rm = TRUE) / sum(f0, na.rm = TRUE)
  
  term1 <- ifelse(is.finite(p10) & p10 > 0, sum(f1 * (p10_i - p10)^2 / p10, na.rm = TRUE), 0)
  term2 <- ifelse(is.finite(p11) & p11 > 0, sum(f1 * (p11_i - p11)^2 / p11, na.rm = TRUE), 0)
  term3 <- ifelse(is.finite(p00) & p00 > 0, sum(f0 * (p00_i - p00)^2 / p00, na.rm = TRUE), 0)
  term4 <- ifelse(is.finite(p01) & p01 > 0, sum(f0 * (p01_i - p01)^2 / p01, na.rm = TRUE), 0)
  
  X2_f <- term1 + term2 + term3 + term4
  
  a0 <- sum(f0 > 0, na.rm = TRUE)
  a1 <- sum(f1 > 0, na.rm = TRUE)
  b0 <- sum(c(p00, p01) > 0, na.rm = TRUE)
  b1 <- sum(c(p10, p11) > 0, na.rm = TRUE)
  
  df <- max(0, (a0 - 1) * (b0 - 1)) + max(0, (a1 - 1) * (b1 - 1))
  pval <- if (df > 0) pchisq(X2_f, df = df, lower.tail = FALSE) else NA_real_
  
  list(stat = X2_f, df = df, pval = pval)
}

convergence_for_dataset <- function(transition_freq_list, on_frequency_list, n_iter_per_chain) {
  S <- length(transition_freq_list)
  stopifnot(S >= 2, S == length(on_frequency_list))
  
  genes_intersect <- Reduce(intersect, lapply(transition_freq_list, rownames))
  genes_intersect <- intersect(genes_intersect, Reduce(intersect, lapply(on_frequency_list, names)))
  
  if (!is.null(GO_names)) {
    genes_intersect <- intersect(genes_intersect, GO_names)
  }
  
  if (length(genes_intersect) == 0) {
    stop("No overlapping gene-set names across reps.")
  }
  
  if (length(n_iter_per_chain) == 1L) n_iter_per_chain <- rep(n_iter_per_chain, S)
  
  res <- lapply(genes_intersect, function(g) {
    tf_mat <- do.call(rbind, lapply(transition_freq_list, function(df) {
      as.numeric(df[g, c("1->0", "1->1", "0->0", "0->1")])
    }))
    colnames(tf_mat) <- c("1->0", "1->1", "0->0", "0->1")
    
    of_vec <- sapply(on_frequency_list, function(v) as.numeric(v[g]))
    
    w_out <- weiss_between_chains(of_vec, tf_mat, n_iter_per_chain)
    b_out <- billingsley_between_chains(tf_mat, n_iter_per_chain)
    
    p01 <- tf_mat[, "0->1"]
    p10 <- tf_mat[, "1->0"]
    phi <- .clip(1 - p01 - p10)
    ess_bin <- ess_binary_from_phi(n_iter_per_chain, phi)
    
    data.frame(
      gene_set      = g,
      weiss_stat    = w_out$stat,
      weiss_df      = w_out$df,
      weiss_p       = w_out$pval,
      bill_stat     = b_out$stat,
      bill_df       = b_out$df,
      bill_p        = b_out$pval,
      phi_mean      = mean(phi, na.rm = TRUE),
      phi_min       = min(phi, na.rm = TRUE),
      phi_max       = max(phi, na.rm = TRUE),
      ess_bin_mean  = mean(ess_bin, na.rm = TRUE),
      ess_bin_min   = min(ess_bin, na.rm = TRUE),
      ess_bin_max   = max(ess_bin, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  
  out_df <- bind_rows(res)
  out_df$weiss_p_bh <- p.adjust(out_df$weiss_p, method = "BH")
  out_df$bill_p_bh  <- p.adjust(out_df$bill_p, method = "BH")
  
  out_df
}

rhat_ess_for_two_params <- function(param1_chains, param2_chains,
                                    param1_name = "para1",
                                    param2_name = "para2") {
  rhat1 <- rhat_basic(param1_chains)
  rhat2 <- rhat_basic(param2_chains)
  
  ess1 <- ess_continuous_chains(param1_chains)
  ess2 <- ess_continuous_chains(param2_chains)
  
  bind_rows(
    data.frame(
      param = param1_name,
      chain = seq_along(ess1),
      ESS = ess1,
      Rhat = rhat1,
      stringsAsFactors = FALSE
    ),
    data.frame(
      param = param2_name,
      chain = seq_along(ess2),
      ESS = ess2,
      Rhat = rhat2,
      stringsAsFactors = FALSE
    )
  )
}

# =========================================================
# FILE PARSING
# =========================================================

parse_result_files <- function(res_dir) {
  files <- list.files(res_dir, pattern = "^result_.*\\.rds$", full.names = TRUE)
  
  if (length(files) == 0) {
    stop("No files found in res_dir with pattern ^result_.*\\.rds$")
  }
  
  meta <- tibble(
    file = files,
    file_base = basename(files),
    stem = str_remove(basename(files), "\\.rds$"),
    model = str_match(stem, "^result_(.+)_gen\\d+_rep\\d+$")[, 2],
    dataset = str_extract(stem, "gen\\d+"),
    dataset_num = as.integer(str_remove(str_extract(stem, "gen\\d+"), "gen")),
    rep = as.integer(str_remove(str_extract(stem, "rep\\d+"), "rep"))
  )
  
  bad <- meta %>% filter(is.na(model) | is.na(dataset) | is.na(rep))
  if (nrow(bad) > 0) {
    stop("Some filenames do not match expected pattern result_<model>_gen<id>_rep<k>.rds:\n",
         paste(bad$file_base, collapse = "\n"))
  }
  
  # Preserve your old dataset-to-(ncp,n_on) mapping
  ncp_values <- c(1, 2, 3, 5, 10)
  n_on_values <- c(5, 10, 20, 50)
  
  meta <- meta %>%
    mutate(
      category_index = (dataset_num - 1) %/% 100 + 1,
      ncp = ncp_values[(category_index - 1) %% 5 + 1],
      n_on = n_on_values[(category_index - 1) %/% 5 + 1]
    ) %>%
    arrange(model, dataset_num, rep)
  
  meta
}

load_one_result <- function(path) {
  obj <- readRDS(path)
  
  trans_info <- normalize_transition_object(obj)
  trans_freq <- trans_info$transition_freq
  
  on_freq <- obj$on_frequency %||% NULL
  if (is.null(on_freq)) stop("Missing on_frequency in: ", basename(path))
  on_freq <- on_freq[!is.na(on_freq)]
  
  keep <- intersect(rownames(trans_freq), names(on_freq))
  if (!is.null(GO_names)) {
    keep <- intersect(keep, GO_names)
  }
  if (length(keep) == 0) {
    stop("No overlapping gene sets between transitions and on_frequency in: ", basename(path))
  }
  
  trans_freq <- trans_freq[keep, , drop = FALSE]
  on_freq <- on_freq[keep]
  
  n_iter <- infer_n_iter(obj, trans_info)
  
  p1 <- obj$p1_trace %||% NULL
  p2 <- obj$p2_trace %||% NULL
  
  if (!is.null(p1)) p1 <- safe_numeric_vector(p1, paste0(basename(path), ": p1_trace"))
  if (!is.null(p2)) p2 <- safe_numeric_vector(p2, paste0(basename(path), ": p2_trace"))
  
  list(
    transition_freq = trans_freq,
    on_frequency = on_freq,
    n_iter = n_iter,
    p1_trace = p1,
    p2_trace = p2
  )
}

# =========================================================
# DATASET-LEVEL ANALYSIS
# =========================================================

analyze_one_dataset <- function(group_df) {
  stopifnot(nrow(group_df) >= 2)
  
  objs <- lapply(group_df$file, load_one_result)
  
  # Basic rep check
  rep_ids <- sort(group_df$rep)
  if (!identical(rep_ids, sort(unique(rep_ids)))) {
    warning("Duplicate reps found for ", unique(group_df$model), " ", unique(group_df$dataset))
  }
  
  tf_list <- lapply(objs, `[[`, "transition_freq")
  of_list <- lapply(objs, `[[`, "on_frequency")
  
  n_iter_vec <- vapply(objs, function(x) x$n_iter %||% NA_integer_, integer(1))
  if (all(is.na(n_iter_vec))) {
    warning("Could not infer n_iter for ", unique(group_df$model), " ", unique(group_df$dataset),
            "; using 5e5 as fallback for binary ESS / Weiss / Billingsley weighting.")
    n_iter_vec <- rep(5e5L, nrow(group_df))
  } else {
    fallback <- max(n_iter_vec, na.rm = TRUE)
    n_iter_vec[is.na(n_iter_vec)] <- fallback
  }
  
  binary_df <- convergence_for_dataset(
    transition_freq_list = tf_list,
    on_frequency_list = of_list,
    n_iter_per_chain = n_iter_vec
  )
  
  genes <- binary_df$gene_set
  onfreq_mat <- do.call(cbind, lapply(of_list, function(v) v[genes]))
  colnames(onfreq_mat) <- paste0("rep", group_df$rep)
  
  onfreq_avg <- data.frame(
    gene_set = genes,
    on_freq_avg = rowMeans(onfreq_mat, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  cont_ok <- all(vapply(objs, function(x) !is.null(x$p1_trace) && !is.null(x$p2_trace), logical(1)))
  
  cont_df <- if (cont_ok) {
    rhat_ess_for_two_params(
      param1_chains = lapply(objs, `[[`, "p1_trace"),
      param2_chains = lapply(objs, `[[`, "p2_trace"),
      param1_name = paste0(unique(group_df$model), "_para1"),
      param2_name = paste0(unique(group_df$model), "_para2")
    )
  } else {
    NULL
  }
  
  summary_df <- data.frame(
    model = unique(group_df$model),
    dataset = unique(group_df$dataset),
    dataset_num = unique(group_df$dataset_num),
    ncp = unique(group_df$ncp),
    n_on = unique(group_df$n_on),
    n_reps = nrow(group_df),
    n_genes = nrow(binary_df),
    n_iter_min = min(n_iter_vec, na.rm = TRUE),
    n_iter_max = max(n_iter_vec, na.rm = TRUE),
    weiss_sig_prop = mean(binary_df$weiss_p_bh < weiss_alpha, na.rm = TRUE),
    bill_sig_prop = mean(binary_df$bill_p_bh < bill_alpha, na.rm = TRUE),
    phi_mean_overall = mean(binary_df$phi_mean, na.rm = TRUE),
    ess_bin_mean_overall = mean(binary_df$ess_bin_mean, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(cont_df)) {
    rhat_summary <- cont_df %>%
      group_by(param) %>%
      summarise(
        Rhat = unique(Rhat)[1],
        ESS_mean = mean(ESS, na.rm = TRUE),
        ESS_min = min(ESS, na.rm = TRUE),
        .groups = "drop"
      )
    
    summary_df <- summary_df %>%
      mutate(
        rhat_para1 = rhat_summary$Rhat[1] %||% NA_real_,
        rhat_para2 = rhat_summary$Rhat[2] %||% NA_real_,
        ess_para1_mean = rhat_summary$ESS_mean[1] %||% NA_real_,
        ess_para2_mean = rhat_summary$ESS_mean[2] %||% NA_real_
      )
  } else {
    summary_df <- summary_df %>%
      mutate(
        rhat_para1 = NA_real_,
        rhat_para2 = NA_real_,
        ess_para1_mean = NA_real_,
        ess_para2_mean = NA_real_
      )
  }
  
  list(
    meta = group_df,
    summary = summary_df,
    binary = binary_df,
    onfreq_avg = onfreq_avg,
    continuous = cont_df
  )
}

# =========================================================
# RUN ALL DATASETS
# =========================================================

meta <- parse_result_files(res_dir)

# write meta
saveRDS(meta, file.path(conv_dir, "meta.rds"))
write_csv(meta, file.path(conv_dir, "meta.csv"))

groups <- meta %>%
  group_by(model, dataset, dataset_num, ncp, n_on) %>%
  group_split()

group_keys <- meta %>%
  group_by(model, dataset, dataset_num, ncp, n_on) %>%
  group_keys()

all_results <- vector("list", length(groups))

for (i in seq_along(groups)) {
  g <- groups[[i]]
  key <- group_keys[i, ]
  
  message("Analyzing: ", key$model, " / ", key$dataset)
  
  res <- analyze_one_dataset(g)
  all_results[[i]] <- res
  
  out_name <- paste0("conv_", key$model, "_", key$dataset, ".rds")
  saveRDS(res, file.path(indiv_dir, out_name))
}

names(all_results) <- paste(group_keys$model, group_keys$dataset, sep = "__")

combined <- list(
  meta = meta,
  group_keys = group_keys,
  results = all_results,
  summary = bind_rows(lapply(all_results, `[[`, "summary"))
)

saveRDS(combined, file.path(conv_dir, "combined_convergence.rds"))
write_csv(combined$summary, file.path(conv_dir, "combined_summary.csv"))

# =========================================================
# TIDY TABLES FOR PLOTTING
# =========================================================

summary_df <- combined$summary

binary_all <- bind_rows(lapply(all_results, function(x) {
  x$binary %>%
    mutate(
      model = x$summary$model,
      dataset = x$summary$dataset,
      dataset_num = x$summary$dataset_num,
      ncp = x$summary$ncp,
      n_on = x$summary$n_on
    )
}))

continuous_all <- bind_rows(lapply(all_results, function(x) {
  if (is.null(x$continuous)) return(NULL)
  x$continuous %>%
    mutate(
      model = x$summary$model,
      dataset = x$summary$dataset,
      dataset_num = x$summary$dataset_num,
      ncp = x$summary$ncp,
      n_on = x$summary$n_on
    )
}))


#################### Ploting

# helper
`%||%` <- function(a, b) if (!is.null(a)) a else b

## =========================================================
## 2) Helpers
## =========================================================
genes_over_0p5 <- function(onfreq_df) {
  if (is.null(onfreq_df) || nrow(onfreq_df) == 0) return(character(0))
  onfreq_df %>%
    filter(on_freq_avg > 0.5) %>%
    pull(gene_set) %>%
    unique()
}

fix_binary_p <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  df$weiss_p_bh <- ifelse(is.na(df$weiss_p_bh) | !is.finite(df$weiss_p_bh), 1, df$weiss_p_bh)
  df$bill_p_bh  <- ifelse(is.na(df$bill_p_bh)  | !is.finite(df$bill_p_bh),  1, df$bill_p_bh)
  df
}

.first_non_na <- function(x) {
  y <- x[!is.na(x)]
  if (length(y)) y[1] else NA_real_
}

extract_rhat_two <- function(model, dataset, ncp, n_on, cont_df) {
  if (is.null(cont_df) || !is.data.frame(cont_df) || nrow(cont_df) == 0) {
    return(tibble())
  }
  
  cont_df %>%
    group_by(param) %>%
    summarise(Rhat = .first_non_na(Rhat), .groups = "drop") %>%
    mutate(
      param_label = case_when(
        str_detect(param, "para1$") ~ "para1",
        str_detect(param, "para2$") ~ "para2",
        param == "mean" ~ "para1",
        param == "sd"   ~ "para2",
        TRUE ~ paste0("param_", param)
      ),
      model   = as.character(model),
      dataset = as.character(dataset),
      ncp     = as.integer(ncp),
      n_on    = as.integer(n_on),
      xlab    = paste(model, param_label, sep = "_")
    ) %>%
    select(model, dataset, ncp, n_on, Rhat, xlab)
}

## =========================================================
## 3) Flatten combined$result list to dataset table
## =========================================================
datasets_tbl <- purrr::imap_dfr(combined$results, function(g, nm) {
  g$binary <- fix_binary_p(g$binary)
  
  tibble(
    model   = as.character(g$summary$model[1]),
    dataset = as.character(g$summary$dataset[1]),
    ncp     = as.integer(g$summary$ncp[1]),
    n_on    = as.integer(g$summary$n_on[1]),
    binary  = list(g$binary),
    onfreq  = list(g$onfreq_avg),
    cont    = list(g$continuous)
  )
})

## =========================================================
## 4) Adjusted p-value aggregation
## same logic as old code
## =========================================================
pval_long <- datasets_tbl %>%
  mutate(genes_gt = purrr::map(onfreq, genes_over_0p5)) %>%
  mutate(
    weiss_vec = purrr::map2(binary, genes_gt, ~{
      if (is.null(.x) || length(.y) == 0) return(numeric(0))
      .x %>% filter(gene_set %in% .y) %>% pull(weiss_p_bh) %>% as.numeric()
    }),
    bill_vec = purrr::map2(binary, genes_gt, ~{
      if (is.null(.x) || length(.y) == 0) return(numeric(0))
      .x %>% filter(gene_set %in% .y) %>% pull(bill_p_bh) %>% as.numeric()
    })
  ) %>%
  select(model, dataset, ncp, n_on, weiss_vec, bill_vec)

pval_pooled_by_model <- pval_long %>%
  group_by(ncp, n_on, model) %>%
  summarise(
    weiss_vals = list(unlist(weiss_vec, use.names = FALSE)),
    bill_vals  = list(unlist(bill_vec,  use.names = FALSE)),
    .groups = "drop"
  )

pval_plot_df <- pval_pooled_by_model %>%
  pivot_longer(
    cols = c(weiss_vals, bill_vals),
    names_to = "metric",
    values_to = "pvals"
  ) %>%
  mutate(
    metric_label = ifelse(metric == "weiss_vals", "weiss", "bill"),
    xlab = paste(model, metric_label, sep = "_"),
    pretty_metric = ifelse(metric_label == "weiss", "Weiss (BH)", "Billingsley (BH)")
  ) %>%
  unnest(pvals) %>%
  mutate(
    pvals = as.numeric(pvals),
    pvals = pmin(pmax(pvals, 0), 1)
  ) %>%
  filter(is.finite(pvals)) %>%
  mutate(
    n_on = factor(n_on, levels = c(5, 10, 20, 50)),
    ncp  = factor(ncp,  levels = c(1, 2, 3, 5, 10))
  )

## proportion below 0.05 for each violin
pval_prop_df <- pval_plot_df %>%
  group_by(ncp, n_on, xlab) %>%
  summarise(
    prop_below_005 = mean(pvals < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("%.2f", prop_below_005),
    y = 0.03
  )

## p-value violin plot
p_violin <- ggplot(pval_plot_df, aes(x = xlab, y = pvals)) +
  geom_violin(trim = FALSE, na.rm = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.6, na.rm = TRUE) +
  geom_hline(yintercept = 0.05, linetype = 2) +
  geom_text(
    data = pval_prop_df,
    aes(x = xlab, y = y, label = label),
    inherit.aes = FALSE,
    color = "red",
    size = 3
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = NULL,
    y = "Adjusted p-value"
  ) +
  facet_grid(
    rows = vars(n_on),
    cols = vars(ncp),
    labeller = labeller(
      n_on = function(x) paste0("m_on=", x),
      ncp  = function(x) paste0("NCP=", x)
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

## =========================================================
## 5) R-hat aggregation
## same logic as old code
## =========================================================
rhat_long <- datasets_tbl %>%
  filter(!purrr::map_lgl(cont, is.null)) %>%
  mutate(
    rhat_rows = purrr::pmap(
      list(model, dataset, ncp, n_on, cont),
      extract_rhat_two
    )
  ) %>%
  select(rhat_rows) %>%
  unnest(rhat_rows) %>%
  mutate(
    n_on = factor(n_on, levels = c(5, 10, 20, 50)),
    ncp  = factor(ncp,  levels = c(1, 2, 3, 5, 10))
  )

## proportion of R-hat > 1.1
rhat_prop_df <- rhat_long %>%
  group_by(ncp, n_on, xlab) %>%
  summarise(
    prop_above_105 = mean(Rhat > 1.1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("%.2f", prop_above_105),
    y = 2.35
  )

## R-hat boxplot
p_rhat <- ggplot(rhat_long, aes(x = xlab, y = Rhat)) +
  geom_boxplot(outlier.size = 0.6) +
  geom_hline(yintercept = 1.1, linetype = 2) +
  geom_text(
    data = rhat_prop_df,
    aes(x = xlab, y = y, label = label),
    inherit.aes = FALSE,
    color = "red",
    size = 3
  ) +
  labs(
    x = NULL,
    y = "R-hat"
  ) +
  facet_grid(
    rows = vars(n_on),
    cols = vars(ncp),
    labeller = labeller(
      n_on = function(x) paste0("m_on=", x),
      ncp  = function(x) paste0("NCP=", x)
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  coord_cartesian(ylim = c(0.9, 2.5))

## =========================================================
## 6) Save
## =========================================================
ggsave(
  file.path(plot_dir, "violin_by_category.png"),
  p_violin, width = 14, height = 4, units = "in", bg = "white"
)

ggsave(
  file.path(plot_dir, "rhat_by_category.png"),
  p_rhat, width = 14, height = 6, units = "in", bg = "white"
)

