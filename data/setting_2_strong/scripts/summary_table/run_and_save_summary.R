source("scripts/summary_table/summary_fitting_Btimes.R")

run_many_reri_in_memory <- function(beta0_vec,
                                    B,
                                    use      = "control",
                                    n_small  = 3000,
                                    seed     = 2025,
                                    digits   = 3,
                                    verbose  = TRUE) {
  
  setNames(
    lapply(beta0_vec, function(b0) {
      if (verbose) message("Running beta0 = ", b0, ", B = ", B)
      
      fit_Btimes(beta0   = b0,
              B       = B,
              n_small = n_small,
              use     = use,
              seed    = seed,
              digits  = digits,
              verbose = verbose)
    }),
    nm <- as.character(beta0_vec)
  )
}


# beta0_vec <- c(-6, -4.5, -3.7, -3, -2, -1)
# B <- 10
# use <- "control"
# 
# reri_list <- run_many_reri_in_memory(beta0_vec = beta0_vec,
#                                      B         = B,
#                                      use       = use)
# 
# file_name <- sprintf("%d_%s_reri.Rdata", B, use)           # → "1000_control_reri.Rdata"
# file_path <- file.path("outputs", "simulation_data", file_name)
# 
# save(reri_list, file = file_path)

# my_results <- get(load("outputs/simulation_data/10_control_reri.Rdata")[1])





