source("scripts/data_generation.R")

compute_population_truth <- function(beta0, n_large = 1e6, seed = 2025){
  set.seed(seed)
  sim_dat_large <- simulate_data_intercept(n_large, beta0)
  
  # --- population truth ----------------------------------------------
  p00 <- with(sim_dat_large, mean(Y[G == 0 & E == 0]))
  p10 <- with(sim_dat_large, mean(Y[G == 1 & E == 0]))
  p01 <- with(sim_dat_large, mean(Y[G == 0 & E == 1]))
  p11 <- with(sim_dat_large, mean(Y[G == 1 & E == 1]))
  
  true_RERI_RR <- (p11 - p10 - p01 + p00) / p00
  
  logit <- function(p) log(p / (1 - p))
  OR11  <- exp(logit(p11) - logit(p00))
  OR10  <- exp(logit(p10) - logit(p00))
  OR01  <- exp(logit(p01) - logit(p00))
  true_RERI_OR <- OR11 - OR10 - OR01 + 1
  return(list(true_marginal_RERI_RR = true_RERI_RR,
              true_marginal_RERI_OR = true_RERI_OR))
}

beta0 <- c(-6, -5.5, -5, -4.5, -4, -3, -2, -1)
# beta0 <- seq(-6, -1, by = 0.1)
# beta0 <- -6

population_truth <- setNames(
  lapply(seq_along(beta0), function(i) {
    res <- compute_population_truth(beta0 = beta0[i])
    message(sprintf("finished b0 = %s  (%d of %d)", beta0[i], i, length(beta0)))
    res
  }),
  as.character(beta0)
)

save(population_truth, file = "scripts/population_truth/population_truth.Rdata")


# ################################################################################
# 
# # -------------------------------------------------------------
# #   Population truth by Monte‑Carlo replication (parallel)
# # -------------------------------------------------------------
# compute_population_truth_parallel <- function(beta0,
#                                               n_large = 1e6,
#                                               B       = 2000,
#                                               seed    = 2025,
#                                               verbose = TRUE) {
#   ## 1. cluster ------------------------------------------------
#   numCores <- min(parallel::detectCores() - 1, 120)
#   cl <- parallel::makeCluster(numCores, type = "PSOCK")
#   registerDoSNOW(cl)
#   on.exit(parallel::stopCluster(cl), add = TRUE)
#   
#   ## 2. progress bar (handled by doSNOW) -----------------------
#   if (verbose && interactive()) {
#     pb <- txtProgressBar(max = B, style = 3)
#     progress <- function(n) setTxtProgressBar(pb, n)
#     opts <- list(progress = progress)
#   } else {
#     opts <- NULL
#   }
#   
#   ## 3. parallel Monte‑Carlo loop ------------------------------
#   set.seed(seed, kind = "L'Ecuyer-CMRG")
#   
#   reri_list <- foreach(
#     b        = seq_len(B),
#     .combine = "rbind",
#     .export  = c("simulate_data_intercept"),
#     .options.snow = opts,
#     .options.RNG  = seed,
#     .packages     = "stats",      # simulate_data_intercept deps
#     .inorder      = TRUE,
#     .errorhandling = "remove"
#   ) %dorng% {
#     dat <- simulate_data_intercept(n_large, beta0)
#     
#     ## --- population truth for one replicate ------------------
#     p00 <- with(dat, mean(Y[G == 0 & E == 0]))
#     p10 <- with(dat, mean(Y[G == 1 & E == 0]))
#     p01 <- with(dat, mean(Y[G == 0 & E == 1]))
#     p11 <- with(dat, mean(Y[G == 1 & E == 1]))
#     
#     reri_rr <- (p11 - p10 - p01 + p00) / p00
#     
#     logit <- function(p) log(p / (1 - p))
#     OR11  <- exp(logit(p11) - logit(p00))
#     OR10  <- exp(logit(p10) - logit(p00))
#     OR01  <- exp(logit(p01) - logit(p00))
#     reri_or <- OR11 - OR10 - OR01 + 1
#     
#     c(true_marginal_RERI_RR = reri_rr, true_marginal_RERI_OR = reri_or)
#   }
#   
#   if (verbose && interactive()) close(pb)
#   
#   ## 4. summarise across B replicates --------------------------
#   colMeans(reri_list, na.rm = TRUE)
# }
# 
# beta0_vals <- c(-6, -4.5, -3.7, -3, -2, -1)
# 
# population_truth <- setNames(
#   lapply(beta0_vals, function(b0) {
#     compute_population_truth_parallel(beta0 = b0,
#                                       n_large = 1e6,
#                                       B       = 2000,
#                                       seed    = 2025,
#                                       verbose = TRUE)
#   }),
#   as.character(beta0_vals)
# )
# 
# save(population_truth,
#      file = "scripts/population_truth/population_truth.Rdata")
# 
