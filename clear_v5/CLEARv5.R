# CLEAR v5.0 - Optimized
# - Implements vectorization caching for massive speedup in Set Update phase
# - Uses functional closures to remove if/else branching from MCMC loop
# - Retains specific priors:
#   - TN: Cauchy(mu), Half-Cauchy(sigma)
#   - Gamma: Exp(shape), Half-Cauchy(rate)
#   - Beta: Beta(a,1), prior Uniform(a)
# 07/21/2025
# Xinglin Jia

CLEAR <- function(genes, stats, GO,
                  stat_type = "stat",
                  model_dist = "truncated_normal",
                  n_iterations = 1e6,
                  burn_in = 5e5,
                  param_update_ratio = 0.2,
                  record_likelihoods = 1000,
                  random_initial = FALSE,
                  penalty_max_bound = 20,
                  record_conf = FALSE) {

  # --- Setup & Validation ---
  stat_type  <- match.arg(stat_type,  choices = c("stat", "p-value"))
  model_dist <- match.arg(model_dist, choices = c("normal", "truncated_normal", "gamma", "beta"))

  dist_fns <- get_dist_handler(stat_type, model_dist)

  # Preprocess
  gene_index <- setNames(seq_along(genes), genes)
  go_indices <- lapply(GO, function(x) unname(gene_index[as.character(x)]))
  n_genes <- length(genes)
  n_sets  <- length(go_indices)

  # Initialize
  init_data <- dist_fns$init(stats, random_initial)
  stats_vec <- init_data$stats
  p1 <- init_data$p1; p2 <- init_data$p2
  refs <- init_data$refs

  current_set_states <- if (random_initial) sample(c(0L,1L), n_sets, replace=TRUE) else integer(n_sets)
  current_gene_states <- get_gene_states(current_set_states, go_indices, n_genes)

  # Pre-calc Likelihoods
  ll_vec_on  <- dist_fns$log_lik_on(stats_vec, p1, p2)
  ll_vec_off <- dist_fns$log_lik_off(stats_vec)
  ll_diff    <- ll_vec_on - ll_vec_off
  base_ll    <- sum(ll_vec_off)

  current_ll_data <- base_ll + sum(ll_diff[current_gene_states == 1])
  current_ll_pen  <- calculate_penalty(current_set_states, penalty_max_bound)
  current_ll_prior <- dist_fns$log_prior(p1, p2, refs)
  current_total_score <- current_ll_data + current_ll_pen + current_ll_prior

  # Recording
  record_idx <- 1L
  rec_interval <- floor(n_iterations / record_likelihoods)

  rec_ll_total <- numeric(record_likelihoods)
  rec_p1 <- numeric(record_likelihoods)
  rec_p2 <- numeric(record_likelihoods)

  # Tuning & Tallying
  proposal_sd <- 0.1
  n_accept <- 0; n_attempt <- 0
  tune_interval <- floor(n_iterations / 1000)

  # Tallies
  on_freq <- numeric(n_sets) # Total iterations spent in state 1
  flip01 <- integer(n_sets)  # Transitions 0 -> 1
  flip10 <- integer(n_sets)  # Transitions 1 -> 0

  # --- MCMC Loop ---
  for (iter in 1:n_iterations) {

    # 1. Update Sets
    if (runif(1) > param_update_ratio) {
      set_id <- sample.int(n_sets, 1)
      old_val <- current_set_states[set_id]
      new_val <- 1L - old_val

      new_set_states <- current_set_states
      new_set_states[set_id] <- new_val
      new_gene_states <- get_gene_states(new_set_states, go_indices, n_genes)

      new_ll_data <- base_ll + sum(ll_diff[new_gene_states == 1])
      new_ll_pen  <- calculate_penalty(new_set_states, penalty_max_bound)

      delta <- (new_ll_data + new_ll_pen) - (current_ll_data + current_ll_pen)

      if (log(runif(1)) < delta) {
        current_set_states  <- new_set_states
        current_gene_states <- new_gene_states
        current_ll_data     <- new_ll_data
        current_ll_pen      <- new_ll_pen
        current_total_score <- new_ll_data + new_ll_pen + current_ll_prior

        # Track Flips
        if (iter > burn_in) {
          if (old_val == 0L) flip01[set_id] <- flip01[set_id] + 1L
          else               flip10[set_id] <- flip10[set_id] + 1L
        }
      }

      # 2. Update Params
    } else {
      n_attempt <- n_attempt + 1
      prop <- dist_fns$propose(p1, p2, proposal_sd)
      p1_new <- prop$p1; p2_new <- prop$p2

      new_ll_prior <- dist_fns$log_prior(p1_new, p2_new, refs)

      if (is.finite(new_ll_prior)) {
        ll_vec_on_new  <- dist_fns$log_lik_on(stats_vec, p1_new, p2_new)
        ll_vec_off_new <- dist_fns$log_lik_off(stats_vec)
        ll_diff_new <- ll_vec_on_new - ll_vec_off_new
        base_ll_new <- sum(ll_vec_off_new)
        new_ll_data <- base_ll_new + sum(ll_diff_new[current_gene_states == 1])

        delta <- (new_ll_data + new_ll_prior) - (current_ll_data + current_ll_prior)

        if (log(runif(1)) < delta) {
          p1 <- p1_new; p2 <- p2_new
          current_ll_prior <- new_ll_prior
          current_ll_data  <- new_ll_data
          ll_diff <- ll_diff_new
          base_ll <- base_ll_new
          current_total_score <- current_ll_data + current_ll_pen + current_ll_prior
          n_accept <- n_accept + 1
        }
      }
    }

    # Recording
    if (iter %% rec_interval == 0 && record_idx <= record_likelihoods) {
      rec_ll_total[record_idx] <- current_total_score
      rec_p1[record_idx] <- p1
      rec_p2[record_idx] <- p2
      record_idx <- record_idx + 1L
    }

    # Tallying
    if (iter > burn_in) on_freq <- on_freq + current_set_states

    # Tuning
    if (iter <= burn_in && iter %% tune_interval == 0) {
      acc_rate <- n_accept / max(n_attempt, 1)
      if (acc_rate > 0.45) proposal_sd <- proposal_sd * 1.1
      if (acc_rate < 0.25) proposal_sd <- proposal_sd * 0.9
      n_accept <- 0; n_attempt <- 0
    }
  }

  # --- Output Construction ---
  valid_idx <- 1:(record_idx - 1L)
  if (length(valid_idx) == 0) warning("No iterations were recorded.")

  T_steps <- n_iterations - burn_in
  on_probs <- if (T_steps > 0) on_freq / T_steps else on_freq
  names(on_probs) <- names(GO)

  # Calculate Transitions
  # N(1->0) and N(0->1) are raw counts from the MCMC
  N10 <- flip10
  N01 <- flip01

  # Derive Diagonals:
  # Total times in State 1 = N(1->1) + N(1->0)
  # Thus: N(1->1) = Total_ON - N(1->0)
  N11 <- pmax(0L, as.integer(round(on_freq - N10)))

  # Total times in State 0 = N(0->0) + N(0->1)
  # Thus: N(0->0) = Total_OFF - N(0->1)
  off_freq <- T_steps - on_freq
  N00 <- pmax(0L, as.integer(round(off_freq - N01)))

  transition_counts <- cbind(`1->0` = N10,
                             `1->1` = N11,
                             `0->0` = N00,
                             `0->1` = N01)
  rownames(transition_counts) <- names(GO)

  transition_freq <- if (T_steps > 0) transition_counts / T_steps else transition_counts

  list(
    on_frequency = on_probs,
    transition_counts = transition_counts,
    transition_freq = transition_freq,
    log_likelihoods = rec_ll_total[valid_idx],
    p1_trace = rec_p1[valid_idx],
    p2_trace = rec_p2[valid_idx]
  )
}

calculate_penalty <- function(set_states, max_bound) {
  m <- length(set_states); n <- sum(set_states)
  if (n == 0) return(0)
  p_min <- 1/m; p_max <- max_bound / m
  p <- max(min(n/m, p_max), p_min)
  n*log(p) + (m-n)*log(1-p)
}

get_gene_states <- function(set_states, go_indices, n_genes) {
  gs <- integer(n_genes)
  active_sets <- which(set_states == 1)
  for (s in active_sets) gs[go_indices[[s]]] <- 1L
  gs
}


get_dist_handler <- function(stat_type, model_dist) {

  # Helper to sanitize p-values
  sanitize_pvals <- function(stats) {
    stats <- as.numeric(stats)
    stats[!is.finite(stats)] <- NA_real_
    stats[stats <= 0] <- .Machine$double.xmin
    stats[stats >= 1] <- 1 - .Machine$double.eps
    return(stats)
  }

  if (model_dist == "truncated_normal") {
    is_pval <- (stat_type == "p-value")
    return(list(
      init = function(stats, random) {
        if (is_pval) {
          stats <- sanitize_pvals(stats)
          s <- -log10(stats)
        } else {
          s <- abs(stats)
        }

        # Handle NAs introduced by sanitization (treat as 0/uninformative)
        s[is.na(s) | !is.finite(s)] <- 0

        mu0 <- mean(s, na.rm=TRUE); scale_ref <- sd(s, na.rm=TRUE)
        if(is.na(scale_ref) || scale_ref <= 0) scale_ref <- 1

        if(random) {
          p1 <- mu0 * runif(1,0.5,1.5); p2 <- scale_ref * runif(1,0.5,1.5)
        } else { p1 <- mu0; p2 <- scale_ref }

        list(stats=s, p1=p1, p2=p2, refs=list(mu0=mu0, s=scale_ref))
      },
      log_lik_on = function(x, p1, p2) {
        if (p2 <= 0) return(rep(-Inf, length(x)))
        dnorm(x, mean=p1, sd=p2, log=TRUE) - log(1 - pnorm(0, mean=p1, sd=p2))
      },
      log_lik_off = function(x) {
        if(is_pval) dexp(x, rate=log(10), log=TRUE) else dnorm(x,0,1,log=TRUE) - log(0.5)
      },
      log_prior = function(p1, p2, refs) {
        if (p1 < 0 || p2 <= 0) return(-Inf)
        dt((p1 - refs$mu0)/2.5, df=1, log=TRUE) + dt(p2/2.5, df=1, log=TRUE)
      },
      propose = function(p1, p2, sd) {
        if(sample(c(TRUE,FALSE),1)) list(p1=rnorm(1,p1,sd), p2=p2)
        else list(p1=p1, p2=rnorm(1,p2,sd))
      }
    ))
  }

  if (model_dist == "gamma") {
    return(list(
      init = function(stats, random) {
        stats <- sanitize_pvals(stats)

        y <- -log(stats)
        # NAs may exist, use na.rm=TRUE
        m <- mean(y, na.rm=TRUE); v <- var(y, na.rm=TRUE)

        # Safety checks for M/V if data is sparse
        if(is.na(m) || m <= 0) m <- 1
        if(is.na(v) || v <= 0) v <- m

        k_init <- m^2/v; theta_init <- v/m

        if(random) {
          p1 <- k_init * runif(1,0.5,1.5); p2 <- theta_init * runif(1,0.5,1.5)
        } else { p1 <- k_init; p2 <- theta_init }

        # NA handling in stats vector: fill with 0 (since log(1) = 0)
        y[is.na(y)] <- 0

        list(stats=y, p1=p1, p2=p2, refs=list())
      },
      log_lik_on = function(x, p1, p2) {
        if(p1<=0 || p2<=0) return(rep(-Inf, length(x)))
        dgamma(x, shape=p1, scale=p2, log=TRUE)
      },
      log_lik_off = function(x) dexp(x, rate=1, log=TRUE),
      log_prior = function(p1, p2, refs) {
        if(p1<=0 || p2<=0) return(-Inf)
        dexp(p1, rate=1, log=TRUE) + dt((1/p2)/2.5, df=1, log=TRUE) - 2*log(p2)
      },
      propose = function(p1, p2, sd) {
        if(sample(c(TRUE,FALSE),1)) list(p1=rnorm(1,p1,sd), p2=p2)
        else list(p1=p1, p2=rnorm(1,p2,sd))
      }
    ))
  }

  if (model_dist == "beta") {
    return(list(
      init = function(stats, random) {
        p <- sanitize_pvals(stats)
        p[is.na(p)] <- 0.5 # Safe fallback for Beta

        a_init <- if(random) runif(1,0.2,0.8) else 0.5
        list(stats=p, p1=a_init, p2=1, refs=list())
      },
      log_lik_on = function(x, p1, p2) {
        if(p1<=0 || p1>=1) return(rep(-Inf, length(x)))
        dbeta(x, shape1=p1, shape2=1, log=TRUE)
      },
      log_lik_off = function(x) rep(0, length(x)),
      log_prior = function(p1, p2, refs) {
        if(p1<=0 || p1>=1) return(-Inf) else 0
      },
      propose = function(p1, p2, sd) {
        list(p1=rnorm(1, mean=p1, sd=sd), p2=1)
      }
    ))
  }
  stop("Unsupported combination")
}
