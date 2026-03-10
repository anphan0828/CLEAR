# Analyze simulation results and generate PR plot
rm(list = ls())

## =========================
## Paths
## =========================
meta_path <- "Your path to meta.rds"
go_path <- "Your path to GO_ecoli.rds"
plot_dir <- "Your path to plot folder"

result_dirs <- list(
  clear_s_tnormal = "Your path to CLEAR_s_tnormal_result/",
  clear_tnormal   = "Your path to CLEAR_tnormal_result/",
  clear_gamma     = "Your path to CLEAR_gamma_result/",
  clear_beta      = "Your path to CLEAR_beta_result/",
  gsea            = "Your path to GSEA_result/",
  hyper           = "Your path to hyper_result/",
  mgsa_fpr05      = "Your path to MGSA_result_fpr.05/",
  mgsa_fpr1       = "Your path to MGSA_result_fpr.1/",
  mgsa_fpr2       = "Your path to MGSA_result_fpr.2/"
)

## =========================
## Load data
## =========================
meta <- readRDS(meta_path)
GO <- readRDS(go_path)

## =========================
## Helper functions
## =========================
extract_clear_res <- function(res_dir, i) {
  res <- readRDS(file.path(res_dir, paste0("gen", i, "_res.rds")))$on_frequency
  sort(res, decreasing = TRUE)
}

extract_mgsa_res <- function(res_dir, i) {
  res <- readRDS(file.path(res_dir, paste0("gen", i, "_res.rds")))@setsResults
  res <- res[order(res$estimate, decreasing = TRUE), ]
  res[is.na(res)] <- 0
  setNames(res$estimate, rownames(res))
}

extract_gsea_res <- function(res_dir, i) {
  res <- readRDS(file.path(res_dir, paste0("gen", i, "_res.rds")))@result
  res <- res[order(res$pvalue, decreasing = FALSE), ]
  res[is.na(res)] <- 0
  setNames(res$pvalue, rownames(res))
}

extract_hyper_res <- function(res_dir, i) {
  res <- readRDS(file.path(res_dir, paste0("gen", i, "_res.rds")))
  res <- res[order(res$P_Value, decreasing = FALSE), ]
  res[is.na(res)] <- 0
  setNames(res$P_Value, res$Term)
}

safe_extract <- function(fun, res_dir, i) {
  tryCatch(fun(res_dir, i), error = function(e) NULL)
}

precision_recall <- function(term_truth, term_pred) {
  precision <- numeric()
  recall <- numeric()
  
  true_positive <- 0
  total_positive <- 0
  
  for (element in term_pred) {
    total_positive <- total_positive + 1
    if (element %in% term_truth) {
      true_positive <- true_positive + 1
    }
    precision <- c(precision, true_positive / total_positive)
    recall <- c(recall, true_positive / length(term_truth))
  }
  
  data.frame(precision = precision, recall = recall)
}

calc_auc <- function(pr_list) {
  unlist(lapply(pr_list, function(pr) {
    auc <- 0
    for (i in 2:nrow(pr)) {
      auc <- auc + (pr$precision[i] + pr$precision[i - 1]) *
        (pr$recall[i] - pr$recall[i - 1]) / 2
    }
    auc
  }))
}

avg_precision_recall <- function(pr_list, recall_points = seq(0, 1, by = 0.01)) {
  interpolated_precisions <- lapply(pr_list, function(pr) {
    if (sum(!is.na(pr$recall)) >= 2 && sum(!is.na(pr$precision)) >= 2) {
      approx(x = pr$recall, y = pr$precision, xout = recall_points, rule = 2)$y
    } else {
      rep(0, length(recall_points))
    }
  })
  
  precision_matrix <- do.call(cbind, interpolated_precisions)
  data.frame(
    precision = rowMeans(precision_matrix, na.rm = TRUE),
    recall = recall_points
  )
}

## =========================
## Find valid indices
## =========================
valid_indices <- c()

for (i in 1:2000) {
  files_exist <- all(
    file.exists(file.path(result_dirs$clear_s_tnormal, paste0("gen", i, "_res.rds"))),
    file.exists(file.path(result_dirs$clear_tnormal,   paste0("gen", i, "_res.rds"))),
    file.exists(file.path(result_dirs$clear_gamma,     paste0("gen", i, "_res.rds"))),
    file.exists(file.path(result_dirs$clear_beta,      paste0("gen", i, "_res.rds"))),
    file.exists(file.path(result_dirs$gsea,            paste0("gen", i, "_res.rds"))),
    file.exists(file.path(result_dirs$hyper,           paste0("gen", i, "_res.rds"))),
    file.exists(file.path(result_dirs$mgsa_fpr05,      paste0("gen", i, "_res.rds"))),
    file.exists(file.path(result_dirs$mgsa_fpr1,       paste0("gen", i, "_res.rds"))),
    file.exists(file.path(result_dirs$mgsa_fpr2,       paste0("gen", i, "_res.rds")))
  )
  
  if (files_exist) {
    valid_indices <- c(valid_indices, i)
  }
}

cat("Number of valid indices:", length(valid_indices), "\n")

## =========================
## Load results
## =========================
res_list <- list()

for (i in valid_indices) {
  res_list[[as.character(i)]] <- list(
    term_truth          = meta$meta$ncps[i],
    clear_s_tnormal_res = safe_extract(extract_clear_res, result_dirs$clear_s_tnormal, i),
    clear_tnormal_res   = safe_extract(extract_clear_res, result_dirs$clear_tnormal, i),
    clear_gamma_res     = safe_extract(extract_clear_res, result_dirs$clear_gamma, i),
    clear_beta_res      = safe_extract(extract_clear_res, result_dirs$clear_beta, i),
    gsea_res            = safe_extract(extract_gsea_res,  result_dirs$gsea, i),
    hyper_res           = safe_extract(extract_hyper_res, result_dirs$hyper, i),
    mgsa_fpr05_res      = safe_extract(extract_mgsa_res,  result_dirs$mgsa_fpr05, i),
    mgsa_fpr1_res       = safe_extract(extract_mgsa_res,  result_dirs$mgsa_fpr1, i),
    mgsa_fpr2_res       = safe_extract(extract_mgsa_res,  result_dirs$mgsa_fpr2, i)
  )
}

## =========================
## PR/AUC analysis
## =========================
pr_meta_function <- function(meta_res_list,
                             gen_meta,
                             term_truths,
                             control_para,
                             control_para_thres,
                             plot_main,
                             plot_dir) {
  library(ggplot2)
  
  true_index <- rep(TRUE, nrow(gen_meta))
  for (k in seq_along(control_para)) {
    true_index <- true_index &
      abs(gen_meta[[control_para[k]]]) >= control_para_thres[[k]][1] &
      abs(gen_meta[[control_para[k]]]) <= control_para_thres[[k]][2]
  }
  
  meta_res_list <- meta_res_list[true_index]
  term_truths <- term_truths[true_index]
  
  method_names <- c(
    "clear_s_tnormal_res",
    "clear_tnormal_res",
    "clear_gamma_res",
    "clear_beta_res",
    "gsea_res",
    "hyper_res",
    "mgsa_fpr05_res",
    "mgsa_fpr1_res",
    "mgsa_fpr2_res"
  )
  
  pr_lists <- setNames(vector("list", length(method_names)), method_names)
  
  for (j in seq_along(method_names)) {
    pr_lists[[j]] <- vector("list", length(meta_res_list))
  }
  
  for (i in seq_along(meta_res_list)) {
    term_truth <- term_truths[[i]]
    for (method in method_names) {
      pr_lists[[method]][[i]] <- precision_recall(
        term_truth,
        names(meta_res_list[[i]][[method]])
      )
    }
  }
  
  auc_list <- lapply(pr_lists, calc_auc)
  avg_pr_list <- lapply(pr_lists, avg_precision_recall)
  
  method_labels <- c(
    clear_s_tnormal_res = "CLEAR (s_tnormal)",
    clear_tnormal_res   = "CLEAR (Tnormal)",
    clear_gamma_res     = "CLEAR (Gamma)",
    clear_beta_res      = "CLEAR (Beta)",
    gsea_res            = "GSEA",
    hyper_res           = "ORA (FPR = 0.05)",
    mgsa_fpr05_res      = "MGSA (FPR = 0.05)",
    mgsa_fpr1_res       = "MGSA (FPR = 0.1)",
    mgsa_fpr2_res       = "MGSA (FPR = 0.2)"
  )
  
  for (method in names(avg_pr_list)) {
    avg_pr_list[[method]]$method <- method_labels[method]
  }
  
  all_data <- do.call(rbind, avg_pr_list)
  
  method_colors <- c(
    "CLEAR (s_tnormal)" = "purple",
    "CLEAR (Tnormal)"   = "#33a02c",
    "CLEAR (Gamma)"     = "#1f78b4",
    "CLEAR (Beta)"      = "blue",
    "GSEA"              = "black",
    "ORA (FPR = 0.05)"  = "grey",
    "MGSA (FPR = 0.05)" = "red",
    "MGSA (FPR = 0.1)"  = "orange",
    "MGSA (FPR = 0.2)"  = "yellow"
  )
  
  method_linetypes <- c(
    "CLEAR (s_tnormal)" = "solid",
    "CLEAR (Tnormal)"   = "solid",
    "CLEAR (Gamma)"     = "solid",
    "CLEAR (Beta)"      = "solid",
    "GSEA"              = "solid",
    "ORA (FPR = 0.05)"  = "solid",
    "MGSA (FPR = 0.05)" = "dashed",
    "MGSA (FPR = 0.1)"  = "dashed",
    "MGSA (FPR = 0.2)"  = "dashed"
  )
  
  p <- ggplot(all_data, aes(x = recall, y = precision, group = method)) +
    geom_line(aes(color = method, linetype = method), linewidth = 1.2) +
    scale_color_manual(values = method_colors) +
    scale_linetype_manual(values = method_linetypes) +
    labs(
      title = plot_main,
      x = "Recall",
      y = "Precision",
      color = "Method",
      linetype = "Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  ggsave(
    filename = file.path(plot_dir, paste0(plot_main, ".png")),
    plot = p,
    width = 8,
    height = 6,
    bg = "white"
  )
  
  data.frame(
    clear_s_tnormal_auc = auc_list$clear_s_tnormal_res,
    clear_tnormal_auc   = auc_list$clear_tnormal_res,
    clear_gamma_auc     = auc_list$clear_gamma_res,
    clear_beta_auc      = auc_list$clear_beta_res,
    gsea_auc            = auc_list$gsea_res,
    hyper_auc           = auc_list$hyper_res,
    mgsa_fpr05_auc      = auc_list$mgsa_fpr05_res,
    mgsa_fpr1_auc       = auc_list$mgsa_fpr1_res,
    mgsa_fpr2_auc       = auc_list$mgsa_fpr2_res
  )
}

auc_df <- pr_meta_function(
  meta_res_list = res_list,
  gen_meta = meta$meta[valid_indices, ],
  term_truths = meta$terms_on,
  control_para = "n_term_on",
  control_para_thres = list(c(5, 50)),
  plot_main = "all data sets",
  plot_dir = plot_dir
)

# Additional analysis for simulation results

## =========================
## Output path
## =========================
plot_path <- "Your path to save avg_auc_heat.png"

## =========================
## Helper function
## =========================
filter_index <- function(gen_meta, control_para, control_para_thres) {
  true_index <- rep(TRUE, nrow(gen_meta))
  
  for (k in seq_along(control_para)) {
    true_index <- true_index &
      abs(gen_meta[[control_para[k]]]) >= control_para_thres[[k]][1] &
      abs(gen_meta[[control_para[k]]]) <= control_para_thres[[k]][2]
  }
  
  which(true_index)
}

## =========================
## Design settings
## =========================
ncp_vals <- c(1, 2, 3, 5, 10)
nterm_vals <- c(5, 10, 20, 50)

design_grid <- expand.grid(
  para_ncp = ncp_vals,
  n_term_on = nterm_vals,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

design_list <- lapply(seq_len(nrow(design_grid)), function(i) {
  ncp <- design_grid$para_ncp[i]
  nterm <- design_grid$n_term_on[i]
  
  list(
    control_para = c("para_ncp", "n_term_on"),
    control_para_thres = list(c(ncp, ncp), c(nterm, nterm)),
    plot_label = paste0("mu_t=", ncp, ", m_on=", nterm)
  )
})

## =========================
## Filter indices
## =========================
index_list <- list()

for (i in seq_along(design_list)) {
  sublist <- design_list[[i]]
  
  index_list[[i]] <- filter_index(
    gen_meta = meta$meta[1:2000, ],
    control_para = sublist$control_para,
    control_para_thres = sublist$control_para_thres
  )
  
  names(index_list)[i] <- sublist$plot_label
}

## =========================
## Map indices to valid datasets
## =========================
index_map <- setNames(seq_along(valid_indices), valid_indices)

transform_vector <- function(vec, map) {
  unname(map[as.character(vec)])
}

transformed_list <- lapply(index_list, transform_vector, map = index_map)

## =========================
## Average AUC by condition
## =========================
avg_auc <- list()

for (i in seq_along(transformed_list)) {
  subdf <- auc_df[transformed_list[[i]], , drop = FALSE]
  avg_auc[[i]] <- colMeans(subdf, na.rm = TRUE)
  names(avg_auc)[i] <- names(transformed_list)[i]
}

avg_auc_df <- as.data.frame(do.call(rbind, avg_auc))

avg_auc_df <- avg_auc_df[, c(
  "clear_tnormal_auc",
  "clear_beta_auc",
  "clear_gamma_auc",
  "clear_s_tnormal_auc",
  "mgsa_fpr05_auc",
  "mgsa_fpr1_auc",
  "mgsa_fpr2_auc",
  "gsea_auc",
  "hyper_auc"
)]

colnames(avg_auc_df) <- c(
  "CLEAR_tnormal",
  "CLEAR_beta",
  "CLEAR_gamma",
  "CLEAR_tnormal_stat",
  "MGSA FPR 0.05",
  "MGSA FPR 0.1",
  "MGSA FPR 0.2",
  "GSEA",
  "ORA"
)

## =========================
## Heatmap plot
## =========================
library(ggplot2)
library(reshape2)
library(viridis)

avg_auc_df_long <- melt(as.matrix(avg_auc_df))
colnames(avg_auc_df_long)[1:2] <- c("Condition", "Method")
avg_auc_df_long$Condition <- factor(
  avg_auc_df_long$Condition,
  labels = names(transformed_list)
)

condition_expr <- as.expression(lapply(design_list, function(x) {
  ncp <- x$control_para_thres[[1]][1]
  nterm <- x$control_para_thres[[2]][1]
  bquote(mu[t] == .(ncp) ~ "," ~ m[plain(on)] == .(nterm))
}))

avg_auc_plot <- ggplot(avg_auc_df_long, aes(Method, Condition, fill = value)) +
  geom_tile() +
  geom_text(aes(
    label = sprintf("%.2f", value),
    color = ifelse(value < 0.5, "white", "black")
  )) +
  scale_color_identity() +
  scale_fill_viridis(option = "A", name = "Avg. PR-AUC") +
  scale_y_discrete(labels = condition_expr) +
  labs(
    x = "Method",
    y = "Data-generating setting"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = plot_path,
  plot = avg_auc_plot,
  width = 6.5,
  height = 6,
  units = "in"
)