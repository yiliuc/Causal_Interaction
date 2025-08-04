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
  return(c(true_marginal_RERI_RR = true_RERI_RR,
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
