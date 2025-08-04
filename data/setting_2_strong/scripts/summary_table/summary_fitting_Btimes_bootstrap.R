source("scripts/summary_table/summary_fitting_once.R")
source("scripts/data_generation.R")
source("scripts/MSLOM.R")
source("scripts/optim/optim_no_covariates.R")
source("scripts/optim/optim_generalized_covariates.R")
source("scripts/DR.R")

library(dplyr)
library(purrr)
library(knitr)
library(kableExtra)

fit_Btimes_bootstrap <- function(beta0,
                       use     = "control",
                       n_small = 3000,
                       B,
                       seed    = 2025,
                       digits  = 3,
                       verbose = TRUE) {      # <- NEW
  
  ## master seed for reproducibility ----
  set.seed(seed)
  original_sample <- simulate_data_intercept(n_small,  # original sample
                                             beta0)
  n <- nrow(original_sample)
  
  set.seed(seed)
  sim_results <- vector("list", B)
  
  ## choose progress style only if user asked for it -------------
  if (verbose) {
    interactive_session <- interactive()
    if (interactive_session) {
      pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
      on.exit(close(pb), add = TRUE)        # close bar on exit
    }
  }
  
  ## Monte-Carlo loop --------------------------------------------
  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    dat <- original_sample[idx, , drop = FALSE]
    
    sim_results[[b]] <-
      fit_once(beta0       = beta0,
               sim_dat     = dat,
               use         = use,
               n_small     = n_small,
               seed        = NULL)$model_betas
    
    ##  progress update  ----
    if (verbose) {
      if (interactive_session) {
        utils::setTxtProgressBar(pb, b)
      } else {
        message(sprintf("Finished iteration %d / %d", b, B))
      }
    }
  }
  
  ## -------- summary (unchanged) --------------------------------
  all_models  <- names(sim_results[[1]])
  coeff_names <- paste0("b", 0:3)
  
  vec <- unlist(sim_results, use.names = FALSE)
  arr <- aperm(array(vec, dim = c(4, length(all_models), B)),
               perm = c(3, 1, 2))
  dimnames(arr) <- list(NULL, coeff_names, all_models)
  
  mean_mat <- apply(arr, c(3, 2), mean, na.rm = TRUE)
  sd_mat   <- apply(arr, c(3, 2),  sd,  na.rm = TRUE)
  
  fmt <- function(m, s)
    sprintf(paste0("%.", digits, "f(%.", digits, "f)"), m, s)
  
  summary_table <- tibble::tibble(
    Model   = all_models,
    b0      = fmt(mean_mat[, "b0"], sd_mat[, "b0"]),
    b1      = fmt(mean_mat[, "b1"], sd_mat[, "b1"]),
    b2      = fmt(mean_mat[, "b2"], sd_mat[, "b2"]),
    b3      = fmt(mean_mat[, "b3"], sd_mat[, "b3"]),
    # RERI_OR = round(population_truth[[as.character(beta0)]]$true_RERI_OR, 3),
    # RERI_RR = round(population_truth[[as.character(beta0)]]$true_RERI_RR, 3)
  )
  
  summary_table
}

# A vanilla non-parametric bootstrap is not ideal here because we need every
# replicate to preserve the target prevalence P(Y = 1).  Sampling rows with
# replacement lets the proportion of cases drift around the original value,
# so the prevalence is no longer strictly controlled.

# fit_Btimes_bootstrap(beta0 = -4.5, B = 200)
