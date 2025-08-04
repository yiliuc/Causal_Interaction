################################################################################
# 1. Score vector  S(beta) = (S0, S1, S2, S3)
score_no_covariates <- function(data, beta) {
  with(data, {
    b0 <- beta[1];  b1 <- beta[2];  b2 <- beta[3];  b3 <- beta[4]
    
    z  <- 1 + b1*G + b2*E + b3*G*E
    O  <- exp(b0) * z
    p  <- O / (1 + O)
    r  <- Y - p
    
    S0 <- sum(r)
    S1 <- sum(r * G / z)
    S2 <- sum(r * E / z)
    S3 <- sum(r * G * E / z)
    
    return(c(S0, S1, S2, S3))
  })
}
################################################################################
# 2.  Observed-information
info_no_covariates <- function(data, beta) {
  with(data, {
    b0 <- beta[1];  b1 <- beta[2];  b2 <- beta[3];  b3 <- beta[4]
    
    z  <- 1 + b1*G + b2*E + b3*G*E
    O  <- exp(b0) * z
    p  <- O / (1 + O)
    r  <- Y - p
    w  <- O / (1 + O)^2
    
    A  <- G / z
    B  <- E / z
    C  <- G * E / z
    
    I <- matrix(0, 4, 4)
    
    I[1,1] <- sum(w)
    I[1,2] <- I[2,1] <- sum(w * A)
    I[1,3] <- I[3,1] <- sum(w * B)
    I[1,4] <- I[4,1] <- sum(w * C)
    
    I[2,2] <- sum((w + r) * A * A)
    I[3,3] <- sum((w + r) * B * B)
    I[4,4] <- sum((w + r) * C * C)
    
    I[2,3] <- I[3,2] <- sum((w + r) * A * B)
    I[2,4] <- I[4,2] <- sum((w + r) * A * C)
    I[3,4] <- I[4,3] <- sum((w + r) * B * C)
    
    dimnames(I) <- list(paste0("b", 0:3), paste0("b", 0:3))
    return(I)
  })
}
################################################################################
# 3. Perform Newton-Raphson
nr_no_covariates <- function(data, beta_start, tol = 1e-7,
                           max_iter = 100, ridge = 1e-6){
  p <- length(beta_start)
  beta_mat  <- matrix(beta_start, nrow = 1)
  score_mat <- matrix(NA_real_, nrow = 1, ncol = p)
  info_arr  <- list()

  beta_old <- beta_start
  for (iter in 1:max_iter) {
    S <- score_no_covariates(data, beta_old)
    I <- info_no_covariates(data,  beta_old)

    delta <- tryCatch(
      solve(I, S),
      error = function(e) solve(I + diag(ridge, p), S)
    )
    beta_new <- beta_old + delta

    beta_mat  <- rbind(beta_mat,  beta_new)
    score_mat <- rbind(score_mat, S)
    info_arr[[iter]] <- I

    if (max(abs(beta_new - beta_old)) < tol) {
      message("Converged in ", iter, " iterations.")
      return(list(beta = beta_new,
                  history = list(beta = beta_mat[-1, , drop = FALSE],
                                 score = score_mat[-1, , drop = FALSE],
                                 info  = info_arr)))
    }
    beta_old <- beta_new
  }
  stop("NR did not converge in ", max_iter, " iterations.")
}

# Sample usage
# set.seed(2025)
# sim_dat <- simulate_data(3000)
# beta_start <- c(-4.5, 0, 0, 2)
# 
# fit <- nr_no_covariates(sim_dat, beta_start)
# fit$beta            # final estimate
# fit$history$beta    # iteration-by-iteration path


################################################################################
# Archived
# library(MASS)
# beta_old <- rep(0.2, 4)
# tolerance <- 1e-7
# iteration <- 0
# diff <- 1
# 
# track_beta <- c(beta_old)
# track_score <- c(1, 1, 1, 1)
# track_information <- c(1, 1, 1, 1)
# 
# 
# while (diff > tolerance) {
#   # Calculate the score vector and information matrix
#   S <- score_no_covariates(sim_dat, beta_old)
#   I <- info_no_covariates(sim_dat, beta_old)
#   
#   # Update beta using the Newton-Raphson step
#   beta_new <- beta_old + solve(I) %*% S
#   
#   diff <- sum(abs(beta_new - beta_old))
#   beta_old <- beta_new
#   
#   # Track the process
#   track_beta <- rbind(track_beta, c(beta_new))
#   track_score <- rbind(track_score, c(S))
#   track_information <- rbind(track_information, c(I))
#   
# }
# 
# beta_old
# track_beta
# track_score
# track_information