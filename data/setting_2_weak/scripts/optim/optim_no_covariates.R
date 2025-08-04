################################################################################
# 1. Score vector
score_no_covariates_optim <- function(data, beta) {
  with(data, {
    b0 <- beta[1]; b1 <- beta[2]; b2 <- beta[3]; b3 <- beta[4]
    z  <- 1 + b1*G + b2*E + b3*G*E
    O  <- exp(b0) * z
    p  <- O / (1 + O)
    r  <- Y - p
    c( S0 = sum(r),
       S1 = sum(r * G / z),
       S2 = sum(r * E / z),
       S3 = sum(r * G * E / z) )
  })
}

################################################################################
# 2. Negative likelihood function
neg_loglik_no_covariates_optim <- function(data, beta) {
  with(data, {
    b0 <- beta[1]; b1 <- beta[2]; b2 <- beta[3]; b3 <- beta[4]
    z  <- 1 + b1*G + b2*E + b3*G*E
    eta <- b0 + log(z)
    return(-sum( Y*eta - log1p(exp(eta)) ))
  })
}

################################################################################
# 3. Fit the optim()
fit_no_covariates_optim <- function(data,
                                    start = c(b0 = 0, b1 = 0, b2 = 0, b3 = 0),
                                    reltol = 1e-9, maxit = 500){
  opt <- optim(par     = start,
               fn      = neg_loglik_no_covariates_optim,
               gr      = function(b, data) - score_no_covariates_optim(data, b),
               data    = data,
               method  = "BFGS",
               control = list(reltol = reltol, maxit = maxit),
               hessian = TRUE)
  
  if (opt$convergence != 0)
    warning("optim() did not converge: code = ", opt$convergence)
  
  par_names <- c("b0", "b1", "b2", "b3")
  beta_hat <- opt$par
  names(beta_hat) <- par_names
  V_beta   <- solve(opt$hessian)
  se_beta  <- sqrt(diag(V_beta))
  
  list(coefficients = beta_hat,
       se           = se_beta,
       vcov         = V_beta,
       optim_details = opt)
}

# set.seed(2025)
# sim_dat <- simulate_data(3000)
# fit <- fit_no_covariates_optim(sim_dat)
# fit$coefficients
# fit$se

# glm(Y ~ G*E, family = binomial, data = sim_dat)$coef
# compute_RERI_OR(0.9616572, 0.6123888, 0.3070633)
