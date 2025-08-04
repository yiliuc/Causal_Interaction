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
  
  vec <- unlist(sim_results, use.names = FALSE) # Return numbers only
  arr <- aperm(array(vec, dim = c(4, length(all_models), B)),
               perm = c(3, 1, 2)) # Return number of models blocks with each block B by 4
  dimnames(arr) <- list(NULL, coeff_names, all_models) # name each block by model name and each column in the block by b1,... 
  
  mean_mat <- apply(arr, c(3, 2), mean, na.rm = TRUE)
  sd_mat   <- apply(arr, c(3, 2),  sd,  na.rm = TRUE)
  
  ci_low  <- mean_mat[, "b3"] - 1.96 * sd_mat[, "b3"]
  ci_high <- mean_mat[, "b3"] + 1.96 * sd_mat[, "b3"]
  
  coverage <- vapply(seq_along(all_models), function(m) {
    b3_vec <- arr[, "b3", m]             # length B
    mean(b3_vec >= ci_low[m] & b3_vec <= ci_high[m])
  }, numeric(1))
  
  fmt <- function(m, s)
    sprintf(paste0("%.", digits, "f(%.", digits, "f)"), m, s)
  
  summary_table <- tibble::tibble(
    Model   = all_models,
    b0      = fmt(mean_mat[, "b0"], sd_mat[, "b0"]),
    b1      = fmt(mean_mat[, "b1"], sd_mat[, "b1"]),
    b2      = fmt(mean_mat[, "b2"], sd_mat[, "b2"]),
    b3      = fmt(mean_mat[, "b3"], sd_mat[, "b3"]),
    coverage = round(coverage, 3)
    # RERI_OR = round(population_truth[[as.character(beta0)]]$true_RERI_OR, 3),
    # RERI_RR = round(population_truth[[as.character(beta0)]]$true_RERI_RR, 3)
  )
  
  summary_table
}

# fit_Btimes(beta0 = -4.5, B = 2)

# make_reri_table(fit_Btimes(beta0 = -4.5, B = 2000), beta0 = -4.5)

################################################################################
# Archived
# library(dplyr)
# library(purrr)        # map_* helpers
# library(knitr)
# library(kableExtra)
# 
# make_reri_table <- function(beta0,
#                             use = "control",
#                             n_small = 3000,
#                             B       = 2000,
#                             seed    = 2025,
#                             digits  = 3) {
#   
#   ## -----------------------------------------------------------
#   ## helper: fit all nine candidate models + DR for ONE sample
#   ## -----------------------------------------------------------
#   fit_once <- function(dat) {
#     
#     fit_nc <- fit_no_covariates_optim(dat)
#     
#     clom <- function(covs, start) fit_kcovariates_optim(dat, covs, start)
#     fit_clom_X1_X2   <- clom(c("X1", "X2"), rep(0, 6))
#     fit_clom_X3_X4   <- clom(c("X3", "X4"), rep(0, 6))
#     fit_clom_X5_X6   <- clom(c("X5", "X6"), rep(0, 6))
#     fit_clom_all_cov <- clom(paste0("X", 1:6), rep(0, 10))
#     
#     mslom <- function(covs = NULL) {
#       dat_w <- if (is.null(covs)) {
#         compute_ipw(dat, use = use)                       # weights with ALL covariates
#       } else {
#         compute_ipw(dat, covars = covs, use = use)        # user-supplied subset
#       }
#       fit_mslom_optim(dat_w)
#     }
#     fit_msl_X1_X2   <- mslom(c("X1", "X2"))
#     fit_msl_X3_X4   <- mslom(c("X3", "X4"))
#     fit_msl_X5_X6   <- mslom(c("X5", "X6"))
#     fit_msl_all_cov <- mslom()
#     
#     dat_w_full <- compute_ipw(dat, use = use)             # DR fit
#     fit_wclom_all_cov <- fit_kcov_w_optim(dat_w_full, covars = paste0("X", 1:6))
#     
#     extract_betas <- function(fit) fit$coefficients[1:4]
#     
#     list(
#       `Default LOM` = extract_betas(fit_nc),
#       `CLOM(X1X2)`  = extract_betas(fit_clom_X1_X2),
#       `CLOM(X3X4)`  = extract_betas(fit_clom_X3_X4),
#       `CLOM(X5X6)`  = extract_betas(fit_clom_X5_X6),
#       `CLOM(all)`   = extract_betas(fit_clom_all_cov),
#       `MSLOM(X1X2)` = extract_betas(fit_msl_X1_X2),
#       `MSLOM(X3X4)` = extract_betas(fit_msl_X3_X4),
#       `MSLOM(X5X6)` = extract_betas(fit_msl_X5_X6),
#       `MSLOM(all)`  = extract_betas(fit_msl_all_cov),
#       `DR`          = extract_betas(fit_wclom_all_cov)
#     )
#   }
#   
#   ## -----------------------------------------------------------
#   ## Monte-Carlo loop
#   ## -----------------------------------------------------------
#   set.seed(seed)
#   sim_results <- vector("list", B)
#   
#   for (b in seq_len(B)) {
#     dat <- simulate_data_intercept(n_small, beta0)
#     sim_results[[b]] <- fit_once(dat)
#   }
#   
#   ## -----------------------------------------------------------
#   ## reshape:  (model × coefficient) × replication  →  two tibbles
#   ## -----------------------------------------------------------
#   all_models   <- names(sim_results[[1]])
#   coeff_names  <- paste0("b", 0:3)
#   
#   ## 3-D array: rep × coef × model
#   vec <- unlist(sim_results, use.names = FALSE)
#   
#   arr <- aperm(                                    # permute dimensions
#     array(vec,
#           dim = c(4, length(all_models), B)
#     ),
#     perm = c(3, 1, 2)
#   )
#   
#   dimnames(arr) <- list(NULL,                      # (replications unnamed)
#                         coeff_names,
#                         all_models)
#   # return(list(results = sim_results, arr =arr))
#   ## empirical mean & sd
#   mean_mat <- apply(arr, c(3, 2), mean,  na.rm = TRUE)
#   sd_mat   <- apply(arr, c(3, 2), sd,    na.rm = TRUE)
#   
#   ## helper to format “mean(sd)”
#   fmt <- function(m, s)
#     sprintf(paste0("%.", digits, "f(%.", digits, "f)"), m, s)
#   
#   summary_tbl <- tibble(
#     Model = all_models,
#     b0 = fmt(mean_mat[, "b0"], sd_mat[, "b0"]),
#     b1 = fmt(mean_mat[, "b1"], sd_mat[, "b1"]),
#     b2 = fmt(mean_mat[, "b2"], sd_mat[, "b2"]),
#     b3 = fmt(mean_mat[, "b3"], sd_mat[, "b3"])
#   )
#   
#   ## -----------------------------------------------------------
#   ## pretty LaTeX table
#   ## -----------------------------------------------------------
#   summary_tbl %>%
#     rename(`$\\beta_0$` = b0,
#            `$\\beta_1$` = b1,
#            `$\\beta_2$` = b2,
#            `$\\beta_3$` = b3) %>%
#     kbl(booktabs = TRUE,
#         escape   = FALSE,
#         align    = c("l", rep("c", 4)),
#         caption  = paste0("Monte-Carlo summary (B = ", B, ")")) %>%
#     kable_styling(position      = "center",
#                   latex_options = "hold_position",
#                   full_width    = FALSE)
# }

