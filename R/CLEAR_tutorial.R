rm(list = ls())

## =========================================================
## CLEAR tutorial script
## =========================================================

## This script demonstrates how to run the CLEAR method on example data.
## If running R from the command line, make sure you are in the root directory of the CLEAR repository

## 1) Set working directory
## Change this to your local CLEAR example folder
work_dir <- "example/" # Change this to your data folder

## Make sure the path ends with "/"
if (!grepl("/$", work_dir)) {
  work_dir <- paste0(work_dir, "/")
}

## 2) Load CLEAR function
source(file.path("R/CLEAR.R"))

## 3) Create results folder (subfolder "results" in the working directory)
result_dir <- file.path(work_dir, "results")
if (!dir.exists(result_dir)) {
  dir.create(result_dir, recursive = TRUE)
}

## 4) Load example data
GO <- readRDS(file.path(work_dir, "GO_ecoli.rds"))
gen <- readRDS(file.path(work_dir, "example_data.rds"))

genes <- gen$genes
p_values <- gen$p_values
t_stat <- gen$t_stat

## 5) Run CLEAR with test statistics
## Normal model for test statistics
cat("Running CLEAR with test statistics (normal model)...\n")
result_s_tnormal <- CLEAR(
  genes,
  t_stat,
  GO,
  n_iterations = 1000000,
  burn_in = 500000
)

saveRDS(
  result_s_tnormal,
  file = file.path(result_dir, "result_stat_normal.rds")
)

## 6) Run CLEAR with p-values
## Truncated normal model
cat("Running CLEAR with p-values (truncated normal model)...\n")
result_tnormal <- CLEAR(
  genes,
  p_values,
  GO,
  stat_type = "p-value"
)

saveRDS(
  result_tnormal,
  file = file.path(result_dir, "result_pvalue_truncated_normal.rds")
)

## Beta model
cat("Running CLEAR with p-values (beta model)...\n")
result_beta <- CLEAR(
  genes,
  p_values,
  GO,
  stat_type = "p-value",
  model_dist = "beta"
)

saveRDS(
  result_beta,
  file = file.path(result_dir, "result_pvalue_beta.rds")
)

## Gamma model
cat("Running CLEAR with p-values (gamma model)...\n")
result_gamma <- CLEAR(
  genes,
  p_values,
  GO,
  stat_type = "p-value",
  model_dist = "gamma"
)

saveRDS(
  result_gamma,
  file = file.path(result_dir, "result_pvalue_gamma.rds")
)

cat("All analyses finished.\n")
cat("Results saved in: ", result_dir, "\n")