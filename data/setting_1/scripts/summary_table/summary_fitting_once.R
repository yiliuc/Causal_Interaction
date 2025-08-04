# load("scripts/population_truth/population_truth.Rdata")


library(dplyr)
library(knitr)
library(kableExtra)

fit_once <- function(beta0,
                     sim_dat = NULL,
                     use     = "control",
                     n_small = 3000,
                     n_large = 1e6,
                     seed    = NULL) {
  
  if (is.null(sim_dat)) {
    if (!is.null(seed)) set.seed(seed)
    sim_dat <- simulate_data_intercept(n_small, beta0)
  }
  
  fit_nc <- fit_no_covariates_optim(sim_dat)
  
  clom <- function(covs, start) fit_kcovariates_optim(sim_dat, covs, start)
  
  # fit_clom_X1_X2   <- clom(c("X1", "X2"),                rep(0, 6))
  # fit_clom_X3_X4   <- clom(c("X3", "X4"),                rep(0, 6))
  # fit_clom_X5_X6   <- clom(c("X5", "X6"),                rep(0, 6))
  fit_clom_all_cov <- clom(c("X1", "X2", "X3", "X4", "X5", "X6"), rep(0, 10))
  
  mslom <- function(covars_G = NULL, covars_E = NULL) {
    dat_w <- if (is.null(covars_G) & is.null(covars_E)) {
      compute_ipw(sim_dat, use = use)
    } else {
      compute_ipw(sim_dat, covars_G = covars_G, covars_E = covars_E, use = use)
    }
    fit_mslom_optim(dat_w)
  }
  # fit_msl_X1_X2   <- mslom(c("X1", "X2"))
  # fit_msl_X3_X4   <- mslom(c("X3", "X4"))
  # fit_msl_X5_X6   <- mslom(c("X5", "X6"))
  fit_msl_corresponding <- mslom(covars_G =  c("X1", "X2", "X5", "X6"),
                                 covars_E =  c("X3", "X4", "X5", "X6"))
  
  dat_w_full       <- compute_ipw(sim_dat, 
                                  covars_G =  c("X1", "X2", "X5", "X6"),
                                  covars_E =  c("X3", "X4", "X5", "X6"), 
                                  use = use)
  fit_wclom_all_cov <- fit_kcov_w_optim(dat_w_full,
                                        covars = paste0("X", 1:6))
  
  extract_betas <- function(fit) fit$coefficients[1:4]
  extract_ses <- function(fit) fit$se[1:4]
  
  fmt_beta <- function(fit, i, digits = 3) {
    est <- round(fit$coefficients[i], digits)
    se  <- round(fit$se[i],          digits)
    sprintf("%.*f(%.*f)", digits, est, digits, se)
  }
  
  models <- list(fit_nc,
                 # fit_clom_X1_X2, fit_clom_X3_X4, fit_clom_X5_X6, 
                 fit_clom_all_cov,
                 # fit_msl_X1_X2,  fit_msl_X3_X4,  fit_msl_X5_X6,  
                 fit_msl_corresponding,
                 fit_wclom_all_cov)
  
  return(list(summary_table = tibble(
    Model   = c("Default LOM",
                "CLOM(all)",
                "MSLOM(corresp.)",
                "DR"),
    b0  = vapply(models, fmt_beta, character(1), i = 1),
    b1  = vapply(models, fmt_beta, character(1), i = 2),
    b2  = vapply(models, fmt_beta, character(1), i = 3),
    b3  = vapply(models, fmt_beta, character(1), i = 4),
    # RERI_OR = round(population_truth[[as.character(beta0)]]$true_RERI_OR, 3),
    # RERI_RR = round(population_truth[[as.character(beta0)]]$true_RERI_RR, 3)
  ),
  model_betas = list(
    `Default LOM` = extract_betas(fit_nc),
    # `CLOM(X1X2)`  = extract_betas(fit_clom_X1_X2),
    # `CLOM(X3X4)`  = extract_betas(fit_clom_X3_X4),
    # `CLOM(X5X6)`  = extract_betas(fit_clom_X5_X6),
    `CLOM(all)`   = extract_betas(fit_clom_all_cov),
    # `MSLOM(X1X2)` = extract_betas(fit_msl_X1_X2),
    # `MSLOM(X3X4)` = extract_betas(fit_msl_X3_X4),
    # `MSLOM(X5X6)` = extract_betas(fit_msl_X5_X6),
    `MSLOM(corresp.)`  = extract_betas(fit_msl_corresponding),
    `DR`          = extract_betas(fit_wclom_all_cov)),
  model_ses = list(
    `Default LOM` = extract_ses(fit_nc),
    `CLOM(all)`   = extract_ses(fit_clom_all_cov),
    `MSLOM(corresp.)`  = extract_ses(fit_msl_corresponding),
    `DR` = extract_ses(fit_wclom_all_cov))))
}


# # ---------------------------------------------------------------------
# # Example --------------------------------------------------------------
# (1) Just the raw numbers as a tibble
# results <- fit_once(beta0 = -4.5)
# 
# # (2) Directly print the LaTeX/HTML table
# make_reri_table(results$summary_table)





################################################################################
# Archived
# 
# 
# library(dplyr)
# library(knitr)
# library(kableExtra)
# make_reri_table <- function(beta0,
#                             use = "control",
#                             n_small = 3000,
#                             n_large = 1e6,
#                             seed    = 2025) {
#   # --- simulate data ---------------------------------------------------------
#   set.seed(seed)
#   sim_dat       <- simulate_data_intercept(n_small, beta0)
#   set.seed(seed)
#   sim_dat_large <- simulate_data_intercept(n_large, beta0)
#   
#   # --- population truth ------------------------------------------------------
#   p00 <- with(sim_dat_large, mean(Y[G == 0 & E == 0]))
#   p10 <- with(sim_dat_large, mean(Y[G == 1 & E == 0]))
#   p01 <- with(sim_dat_large, mean(Y[G == 0 & E == 1]))
#   p11 <- with(sim_dat_large, mean(Y[G == 1 & E == 1]))
#   
#   true_RERI_RR <- (p11 - p10 - p01 + p00) / p00
#   
#   logit <- function(p) log(p / (1 - p))
#   OR11  <- exp(logit(p11) - logit(p00))
#   OR10  <- exp(logit(p10) - logit(p00))
#   OR01  <- exp(logit(p01) - logit(p00))
#   true_RERI_OR <- OR11 - OR10 - OR01 + 1
#   
#   # --- model fitting ---------------------------------------------------------
#   fit_nc <- fit_no_covariates_optim(sim_dat)
#   
#   clom <- function(covs, start) fit_kcovariates_optim(sim_dat, covs, start)
#   
#   fit_clom_X1_X2   <- clom(c("X1", "X2"), rep(0, 6))
#   fit_clom_X3_X4   <- clom(c("X3", "X4"), rep(0, 6))
#   fit_clom_X5_X6   <- clom(c("X5", "X6"), rep(0, 6))
#   fit_clom_all_cov <- clom(c("X1", "X2", "X3", "X4", "X5", "X6"), rep(0, 10))
#   
#   mslom <- function(covs = NULL) {
#     dat_w <- if (is.null(covs)) {
#       compute_ipw(sim_dat, use = use)                 # let compute_ipw use all covariates
#     } else {
#       compute_ipw(sim_dat, covars = covs, use = use)  # user-supplied subset
#     }
#     fit_mslom_optim(dat_w)
#   }
#   fit_msl_X1_X2   <- mslom(c("X1", "X2"))
#   fit_msl_X3_X4   <- mslom(c("X3", "X4"))
#   fit_msl_X5_X6   <- mslom(c("X5", "X6"))
#   fit_msl_all_cov <- mslom()
#   
#   dat_w_full <- compute_ipw(sim_dat, use = use)       # weights using all six X’s
#   fit_wclom_all_cov <- fit_kcov_w_optim(
#     dat_w_full,
#     covars = paste0("X", 1:6)
#   )
#   
#   # helper for “est(se)”
#   fmt_beta <- function(fit, i, digits = 3) {
#     est <- round(fit$coefficients[i], digits)
#     se  <- round(fit$se[i],          digits)
#     sprintf("%.*f(%.*f)", digits, est, digits, se)
#   }
#   
#   # --- assemble table --------------------------------------------------------
#   models <- list(fit_nc,
#                  fit_clom_X1_X2, fit_clom_X3_X4, fit_clom_X5_X6, fit_clom_all_cov,
#                  fit_msl_X1_X2,  fit_msl_X3_X4,  fit_msl_X5_X6,  fit_msl_all_cov,
#                  fit_wclom_all_cov)
#   
#   summary_tbl <- tibble(
#     Model   = c("Default LOM",
#                 "CLOM(X1X2)", "CLOM(X3X4)", "CLOM(X5X6)", "CLOM(all)",
#                 "MSLOM(X1X2)", "MSLOM(X3X4)", "MSLOM(X5X6)", "MSLOM(all)",
#                 "DR"),
#     b0  = vapply(models, fmt_beta, character(1), i = 1),
#     b1  = vapply(models, fmt_beta, character(1), i = 2),
#     b2  = vapply(models, fmt_beta, character(1), i = 3),
#     b3  = vapply(models, fmt_beta, character(1), i = 4),
#     RERI_OR = round(true_RERI_OR, 3),
#     RERI_RR = round(true_RERI_RR, 3)
#   )
#   
#   # --- kable -----------------------------------------------------------------
#   summary_tbl %>%
#     rename(`$\\beta_0$` = b0,
#            `$\\beta_1$` = b1,
#            `$\\beta_2$` = b2,
#            `$\\beta_3$` = b3,
#            `$RERI_{OR}^{true}$` = RERI_OR,
#            `$RERI_{RR}^{true}$` = RERI_RR) %>%
#     kbl(booktabs = TRUE,
#         escape   = FALSE,
#         align    = c("l", rep("c", 6))) %>%
#     kable_styling(position      = "center",
#                   latex_options = "hold_position",
#                   full_width    = FALSE)
# }
# make_reri_table(beta0 = -6)