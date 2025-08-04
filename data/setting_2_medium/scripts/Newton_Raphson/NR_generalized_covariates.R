################################################################################
# 1. Check the data
linodds_check_extract <- function(data, covars) {
  if (!is.data.frame(data))
    stop("'data' must be a data.frame")
  miss <- setdiff(c("Y","G","E", covars), names(data))
  if (length(miss))
    stop("missing columns in 'data': ", paste(miss, collapse = ", "))
  if (anyDuplicated(covars))
    stop("'covars' contains duplicates")
  
  ## return a list with ready-to-use vectors / matrix
  list(
    Y  = data$Y,
    G  = data$G,
    E  = data$E,
    X  = as.matrix(data[, covars, drop = FALSE])   # n × k
  )
}

################################################################################
# 2. Score vector
score_kcov <- function(data, theta, covars = character())
{
  k <- length(covars)
  if (length(theta) != 4 + k)
    stop("'theta' must have length 4 + length(covars)")
  
  lst <- linodds_check_extract(data, covars)
  Y <- lst$Y;  G <- lst$G;  E <- lst$E;  X <- as.matrix(lst$X)
  
  b0 <- theta[1]; b1 <- theta[2]; b2 <- theta[3]; b3 <- theta[4]
  g  <- if (k) theta[5:(4+k)] else numeric(0)
  
  z  <- 1 + b1*G + b2*E + b3*G*E
  eta<- b0 + if(k) X %*% g else 0
  O  <- exp(eta) * z
  p  <- O / (1 + O)
  r  <- as.vector(Y - p)
  
  S_beta  <- c(
    sum(r),
    sum(r * G      / z),
    sum(r * E      / z),
    sum(r * G * E  / z)
  )
  
  if (k) {
    S_gamma <- as.vector( colSums(r * X) )      # length k
    S_out   <- c(S_beta, S_gamma)
    names(S_out) <- c(paste0("b", 0:3), paste0("g_", covars))
  } else {
    S_out   <- S_beta
    names(S_out) <- paste0("b", 0:3)
  }
  
  S_out
}

################################################################################
# 3. Observed information
info_kcov <- function(data, theta, covars = character())
{
  k <- length(covars)
  if (length(theta) != 4 + k)
    stop("'theta' must have length 4 + length(covars)")
  
  lst <- linodds_check_extract(data, covars)
  Y <- lst$Y;  G <- lst$G;  E <- lst$E;  X <- as.matrix(lst$X)
  
  b0 <- theta[1]; b1 <- theta[2]; b2 <- theta[3]; b3 <- theta[4]
  g  <- if (k) theta[5:(4+k)] else numeric(0)
  
  z  <- 1 + b1*G + b2*E + b3*G*E
  eta<- b0 + if(k) X %*% g else 0
  O  <- exp(eta) * z
  p  <- O / (1 + O)
  r  <- as.vector(Y - p)
  w  <- as.vector(O / (1 + O)^2)
  
  A <- G      / z
  B <- E      / z
  C <- G * E  / z
  
  I_betabeta <- matrix(0, 4, 4)
  I_betabeta[1,1] <- sum(w)
  I_betabeta[1,2:4] <- I_betabeta[2:4,1] <-
    c(sum(w * A), sum(w * B), sum(w * C))
  
  wr <- w + r
  I_betabeta[2,2] <- sum(wr * A * A)
  I_betabeta[3,3] <- sum(wr * B * B)
  I_betabeta[4,4] <- sum(wr * C * C)
  
  I_betabeta[2,3] <- I_betabeta[3,2] <- sum(wr * A * B)
  I_betabeta[2,4] <- I_betabeta[4,2] <- sum(wr * A * C)
  I_betabeta[3,4] <- I_betabeta[4,3] <- sum(wr * B * C)
  
  if (k) {
    Vbeta  <- cbind(1, A, B, C)
    I_betagamma <- t( Vbeta ) %*% (w * X)
    
    I_gammagamma <- crossprod( sqrt(w) * X )
  }
  
  d <- 4 + k
  I <- matrix(0, d, d)
  
  I[1:4, 1:4] <- I_betabeta
  
  if (k) {
    I[1:4, 5:d] <- I_betagamma
    I[5:d, 1:4] <- t(I_betagamma)
    I[5:d, 5:d] <- I_gammagamma
  }
  parNames <- c(paste0("b", 0:3),
                if (k) paste0("g_", covars))
  dimnames(I) <- list(parNames, parNames)
  I
}

################################################################################
# 4. Running Newton-Raphson
nr_k_covariates <- function(data, theta_start, covars, tol = 1e-7,
                              max_iter = 100, ridge = 1e-6){
  p <- length(theta_start)
  theta_mat  <- matrix(theta_start, nrow = 1)
  score_mat <- matrix(NA_real_, nrow = 1, ncol = p)
  info_arr  <- list()
  
  theta_old <- theta_start
  for (iter in 1:max_iter) {
    S <- score_kcov(data, theta_old, covars)
    I <- info_kcov(data,  theta_old, covars)
    
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
# sim_dat <- simulate_data_joint(3000)
# theta_start <- c(-4, 0, 0, 4.3, 0, 0, 0, 0, 0, 0)
# 
# covars <- c("X1", "X2", "X3", "X4", "X5", "X6")
# fit <- nr_k_covariates(sim_dat, theta_start, covars)
# fit$theta            # final estimate
# fit$history$beta    # iteration-by-iteration path