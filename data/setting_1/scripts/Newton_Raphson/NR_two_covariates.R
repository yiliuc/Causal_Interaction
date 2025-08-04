################################################################################
# 1. Score vector
score_2cov <- function(data, theta, covars){
  if (length(covars) != 2)
    stop("‘covars’ must be a character vector of length 2.")
  if (!all(covars %in% names(data)))
    stop("Not all covariate names are columns in ‘data’.")
  
  Y  <- data$Y;  G <- data$G;  E <- data$E
  X1 <- data[[covars[1]]]
  X2 <- data[[covars[2]]]
  
  b0 <- theta[1]; b1 <- theta[2]; b2 <- theta[3]; b3 <- theta[4]
  g1 <- theta[5]; g2 <- theta[6]
  
  z <- 1 + b1*G + b2*E + b3*G*E
  O <- exp(b0 + g1*X1 + g2*X2) * z
  p <- O / (1 + O)
  r <- Y - p
  
  S0 <- sum(r)
  S1 <- sum(r * G      / z)
  S2 <- sum(r * E      / z)
  S3 <- sum(r * G * E  / z)
  S4 <- sum(r * X1)
  S5 <- sum(r * X2)
  
  out <- c(S0, S1, S2, S3, S4, S5)
  names(out) <- c("b0","b1","b2","b3",
                  paste0("g_", covars[1]),
                  paste0("g_", covars[2]))
  out
}

################################################################################
# 2. Observed-information matrix
info_2cov <- function(data, theta, covars)
{
  if (length(covars) != 2)
    stop("‘covars’ must be a character vector of length 2.")
  if (!all(covars %in% names(data)))
    stop("Not all covariate names are columns in ‘data’.")
  
  Y  <- data$Y;  G <- data$G;  E <- data$E
  X1 <- data[[covars[1]]]
  X2 <- data[[covars[2]]]
  
  b0 <- theta[1]; b1 <- theta[2]; b2 <- theta[3]; b3 <- theta[4]
  g1 <- theta[5]; g2 <- theta[6]
  
  z <- 1 + b1*G + b2*E + b3*G*E
  O <- exp(b0 + g1*X1 + g2*X2) * z
  p <- O / (1 + O)
  r <- Y - p
  w <- O / (1 + O)^2
  
  A <- G / z
  B <- E / z
  C <- G * E / z
  
  I <- matrix(0, 6, 6)
  I[1,1] <- sum(w)
  I[1,2] <- I[2,1] <- sum(w * A)
  I[1,3] <- I[3,1] <- sum(w * B)
  I[1,4] <- I[4,1] <- sum(w * C)
  I[1,5] <- I[5,1] <- sum(w * X1)
  I[1,6] <- I[6,1] <- sum(w * X2)
  
  I[2,2] <- sum((w + r) * A * A)
  I[3,3] <- sum((w + r) * B * B)
  I[4,4] <- sum((w + r) * C * C)
  
  I[2,3] <- I[3,2] <- sum((w + r) * A * B)
  I[2,4] <- I[4,2] <- sum((w + r) * A * C)
  I[3,4] <- I[4,3] <- sum((w + r) * B * C)
  
  I[2,5] <- I[5,2] <- sum(w * A * X1)
  I[2,6] <- I[6,2] <- sum(w * A * X2)
  I[3,5] <- I[5,3] <- sum(w * B * X1)
  I[3,6] <- I[6,3] <- sum(w * B * X2)
  I[4,5] <- I[5,4] <- sum(w * C * X1)
  I[4,6] <- I[6,4] <- sum(w * C * X2)
  
  I[5,5] <- sum(w * X1 * X1)
  I[6,6] <- sum(w * X2 * X2)
  I[5,6] <- I[6,5] <- sum(w * X1 * X2)
  
  parNames <- c("b0","b1","b2","b3",
                paste0("g_", covars[1]),
                paste0("g_", covars[2]))
  dimnames(I) <- list(parNames, parNames)
  I
}
################################################################################
# 3. Running Newton-Raphson
nr_two_covariates <- function(data, theta_start, covars, tol = 1e-7,
                             max_iter = 100, ridge = 1e-6){
  p <- length(theta_start)
  theta_mat  <- matrix(theta_start, nrow = 1)
  score_mat <- matrix(NA_real_, nrow = 1, ncol = p)
  info_arr  <- list()
  
  theta_old <- theta_start
  for (iter in 1:max_iter) {
    S <- score_2cov(data, theta_old, covars)
    I <- info_2cov(data,  theta_old, covars)
    
    delta <- tryCatch(
      solve(I, S),
      error = function(e) solve(I + diag(ridge, p), S)
    )
    theta_new <- theta_old + delta
    
    theta_mat  <- rbind(theta_mat,  theta_new)
    score_mat <- rbind(score_mat, S)
    info_arr[[iter]] <- I
    
    if (max(abs(theta_new - theta_old)) < tol) {
      message("Converged in ", iter, " iterations.")
      return(list(theta = theta_new,
                  history = list(theta = theta_mat[-1, , drop = FALSE],
                                 score = score_mat[-1, , drop = FALSE],
                                 info  = info_arr)))
    }
    theta_old <- theta_new
  }
  stop("NR did not converge in ", max_iter, " iterations.")
}

# # Sample usage
# set.seed(2025)
# sim_dat <- simulate_data(3000)
# theta_start <- c(-4, 1.6, 0, 0, 0, 0)
# 
# covars <- c("X1", "X2")
# fit <- nr_two_covariates(sim_dat, theta_start, covars)
# fit$theta            # final estimate
# fit$history$beta    # iteration-by-iteration path