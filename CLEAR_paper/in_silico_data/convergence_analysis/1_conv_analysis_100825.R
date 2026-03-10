# convergence analysis of CLEARv3.2
# Xinglin
# 10/08/2025

set.seed(1)
library(foreach)
library(doParallel)
source("Your path to CLEAR.R")
GO <- readRDS("Your path to GO_ecoli.rds")

nCores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
myCluster <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(myCluster)


ite <- 1000000 # number of iterations



data_dir <- "Your path to data/"


for (batch in 1:5){

  # randomly select one data set from each of the ncp and n_on combinations
  replicate <- 100
  n_blocks <- 20
  gen_num <- sapply(0:(n_blocks - 1), function(i) {
    sample((i * replicate + 1):(i * replicate + replicate), 1)
  })

  # run methods and record the configurations
  save_dir <- paste0("Your path", batch, "/")
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  foreach(i = gen_num, .packages = c("mgsa")) %dopar% {
    gen <- readRDS(paste0(data_dir, "gen", i, ".rds"))
    genes <- gen$genes
    stats <- gen$p_values
    t_stat <- gen$t_stat

    for (j in 1:5) {
      result_tnormal <- CLEAR(genes, stats, GO, n_iterations = ite, burn_in = ite/2,
                              stat_type = "p-value", record_conf = T)
      saveRDS(result_tnormal, paste0(save_dir, "result_tnormal_gen", i, "_rep", j, ".rds"))

      result_gamma <- CLEAR(genes, stats, GO, n_iterations = ite, burn_in = ite/2,
                            stat_type = "p-value", model_dist = "gamma", record_conf = T)
      saveRDS(result_gamma, paste0(save_dir, "result_gamma_gen", i, "_rep", j, ".rds"))

      result_s_tnormal <- CLEAR(genes, t_stat, GO, n_iterations = ite, burn_in = ite/2, record_conf = T)
      saveRDS(result_s_tnormal, paste0(save_dir, "result_s_tnormal_gen", i, "_rep", j, ".rds"))

      result_beta <- CLEAR(genes, stats, GO, n_iterations = ite, burn_in = ite/2,
                           stat_type = "p-value", model_dist = "beta", record_conf = T)
      saveRDS(result_beta, paste0(save_dir, "result_beta_gen", i, "_rep", j, ".rds"))
    }
  }
}

stopCluster(myCluster)



