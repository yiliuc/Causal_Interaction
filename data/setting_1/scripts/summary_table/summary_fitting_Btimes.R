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

fit_Btimes <- function(beta0,
                       use     = "control",
                       n_small = 3000,
                       B,
                       seed    = 2025,
                       digits  = 3,
                       verbose = TRUE) {      # <- NEW
  
  set.seed(seed)
  sim_results <- vector("list", B)
  
  if (verbose) {
    interactive_session <- interactive()
    if (interactive_session) {
      pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
      on.exit(close(pb), add = TRUE)
    }
  }
  
  for (b in seq_len(B)) {
    dat <- simulate_data_intercept(n_small, beta0)
    
    sim_results[[b]] <-
      fit_once(beta0       = beta0,
               sim_dat     = dat,
               use         = use,
               n_small     = n_small,
               seed        = NULL)$model_betas
    
    if (verbose) {
      if (interactive_session) {
        utils::setTxtProgressBar(pb, b)
      } else {
        message(sprintf("Finished iteration %d / %d", b, B))
      }
    }
  }
  
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

# fit_Btimes(beta0 = -4.5, B = 2)

# make_reri_table(fit_Btimes(beta0 = -4.5, B = 2000), beta0 = -4.5)

################################################################################
fit_Btimes <- function(beta0,
                       use     = "control",
                       n_small = 3000,
                       B,
                       seed    = 2025,
                       digits  = 3,
                       verbose = TRUE) {
  
  set.seed(seed)
  
  ## ---- containers -------------------------------------------------------
  all_models   <- c("Default LOM", "CLOM(all)", "MSLOM(corresp.)", "DR")
  coeff_names  <- paste0("b", 0:3)
  n_models     <- length(all_models)
  
  est_arr   <- array(NA_real_, dim = c(B, 4, n_models),
                     dimnames = list(NULL, coeff_names, all_models))
  se_arr    <- est_arr                                   # same shape
  cover_b3  <- matrix(NA_real_, nrow = B, ncol = n_models,
                      dimnames = list(NULL, all_models)) # only b3
  
  ## ---- progress bar -----------------------------------------------------
  if (verbose && interactive()) {
    pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    on.exit(close(pb), add = TRUE)
  }
  
  ## ---- true parameter vector (only need the 4th entry) ------------------
  true_b3    <- 2.64
  
  ## ---- Monte-Carlo loop -------------------------------------------------
  for (b in seq_len(B)) {
    
    sim_dat <- simulate_data_intercept(n_small, beta0)
    
    fit   <- fit_once(beta0, sim_dat = sim_dat, use = use)
    betas <- fit$model_betas
    ses   <- fit$model_ses
    
    for (m in seq_along(all_models)) {
      est_arr[b, , m] <- betas[[m]]
      se_arr[b,  , m] <- ses[[m]]
      
      ## ----- coverage for b3 only -----
      low  <- betas[[m]][4] - 1.96 * ses[[m]][4]
      high <- betas[[m]][4] + 1.96 * ses[[m]][4]
      cover_b3[b, m] <- (true_b3 >= low) & (true_b3 <= high)
    }
    
    if (verbose && interactive())
      utils::setTxtProgressBar(pb, b)
  }
  
  ## ---- summaries --------------------------------------------------------
  mean_mat <- apply(est_arr, c(3, 2), mean, na.rm = TRUE)  # model × β
  sd_mat   <- apply(est_arr, c(3, 2),  sd,  na.rm = TRUE)
  cov_vec  <- colMeans(cover_b3, na.rm = TRUE)              # model
  
  fmt_est <- function(m, s)
    sprintf(paste0("%.", digits, "f(%.", digits, "f)"), m, s)
  fmt_cov <- function(x)
    sprintf(paste0("%.", digits, "f"), x)
  
  tibble::tibble(
    Model  = all_models,
    b0     = fmt_est(mean_mat[, "b0"], sd_mat[, "b0"]),
    b1     = fmt_est(mean_mat[, "b1"], sd_mat[, "b1"]),
    b2     = fmt_est(mean_mat[, "b2"], sd_mat[, "b2"]),
    b3     = fmt_est(mean_mat[, "b3"], sd_mat[, "b3"]),
    cov_b3 = fmt_cov(cov_vec)
  )
}
