source("scripts/data_generation.R")

############################################################
# 1. Compute the weight w = 1/p(G|C)·1/p(E|G,C)=1/(p(G,E|C))
compute_ipw <- function(dat,
                        covars_G = paste0("X", 1:6),
                        covars_E = paste0("X", 1:6),
                        use = "control",
                        eps          = 1e-6,
                        cutoff   = 20) {
  stopifnot(all(c("Y", "G", "E", covars_G, covars_E) %in% names(dat)))
  
  # Choose the data to fit the model
  used_data <- if (use == "control"){
    dat[dat$Y == 0, ]
  } else if (use == "both") {
    dat <- dat
  }
  
  # P(G = 1 | C)
  f_G   <- reformulate(covars_G, response = "G")
  mod_G <- glm(f_G, data = used_data, family = binomial())
  
  # P(E = 1 | G, C)
  f_E   <- reformulate(c("G", covars_E), response = "E")
  mod_E <- glm(f_E, data = used_data, family = binomial())
  
  pG1 <- predict(mod_G, newdata = dat, type = "response")
  pE1 <- predict(mod_E, newdata = dat, type = "response")
  
  # pG1 <- pmin(pmax(pG1, eps), 1 - eps)
  # pE1 <- pmin(pmax(pE1, eps), 1 - eps)
  
  # P(G = g_i | C)
  pG  <- ifelse(dat$G == 1, pG1, 1 - pG1)
  
  ## P(E = e_i | G, C)
  pE  <- ifelse(dat$E == 1, pE1, 1 - pE1)
  
  w   <- 1 / (pG * pE)
  w <- pmin(w, cutoff)
  
  # if (!is.null(truncate_q)) {
  #   q_hi <- quantile(w, truncate_q, na.rm = TRUE)
  #   q_lo <- quantile(w, 1 - truncate_q, na.rm = TRUE)
  #   w    <- pmin(pmax(w, q_lo), q_hi)
  # }
  
  dat$w <- w
  dat
}

# sim_dat <- simulate_data_intercept(3000, -1, 2025)
# covars_G =  c("X1", "X2", "X5", "X6")
# covars_E =  c("X3", "X4", "X5", "X6")
# w <- compute_ipw(sim_dat, covars_G = covars_G, covars_E = covars_E, use = "control")
# w[w$w >= 20,]
# 
# sim_dat <- simulate_data_intercept(3000, -6, 2025)
# covars_G =  c("X1", "X2", "X5", "X6")
# covars_E =  c("X3", "X4", "X5", "X6")
# w2 <- compute_ipw(sim_dat, covars_G = covars_G, covars_E = covars_E, use = "control")
# w2[w2$w >= 20,]
# 
# sim_dat <- simulate_data_intercept(3000, -2, 2025)
# covars_G =  c("X1", "X2", "X5", "X6")
# covars_E =  c("X3", "X4", "X5", "X6")
# w3 <- compute_ipw(sim_dat, covars_G = covars_G, covars_E = covars_E, use = "control")
# w3[w3$w >= 20,]
############################################################
# 2. Weighted score and negative likelihood
score_mslom <- function(data, beta) {
  with(data, {
    b0 <- beta[1]; b1 <- beta[2]; b2 <- beta[3]; b3 <- beta[4]
    
    z  <- 1 + b1*G + b2*E + b3*G*E           # z_i  (must stay > 0)
    O  <- exp(b0) * z                        # odds
    p  <- O / (1 + O)                        # risk
    r  <- Y - p                              # residual
    
    c( S0 = sum(w * r),
       S1 = sum(w * r * G      / z),
       S2 = sum(w * r * E      / z),
       S3 = sum(w * r * G * E  / z) )
  })
}

neg_loglik_mslom <- function(data, beta) {
  with(data, {
    b0  <- beta[1]; b1 <- beta[2]; b2 <- beta[3]; b3 <- beta[4]
    z   <- 1 + b1*G + b2*E + b3*G*E
    eta <- b0 + log(z)
    -sum( w * (Y*eta - log1p(exp(eta))) )
  })
}

############################################################
# 3. Optimiser: the inverse of fisher information
fit_mslom_optim <- function(data,
                            start  = c(b0 = 0, b1 = 0, b2 = 0, b3 = 0),
                            reltol = 1e-9,
                            maxit  = 500) {
  
  opt <- optim(par     = start,
               fn      = neg_loglik_mslom,
               gr      = function(b, data) -score_mslom(data, b),
               data    = data,
               method  = "BFGS",
               control = list(reltol = reltol, maxit = maxit),
               hessian = TRUE)
  
  if (opt$convergence != 0)
    warning("optim() did not converge: code = ", opt$convergence)
  
  par_names <- c("b0", "b1", "b2", "b3")
  beta_hat  <- opt$par
  names(beta_hat) <- par_names
  
  V_beta    <- solve(opt$hessian)
  se_beta   <- sqrt(diag(V_beta))
  
  list(coefficients  = beta_hat,
       se            = se_beta,
       vcov          = V_beta,
       optim_details = opt)
}

# dat_w <- compute_ipw(sim_dat)                 # add column w
# fit   <- fit_mslom_optim(dat_w)           # weighted MSLOM fit
# fit$coefficients
# fit$se

##############################################################################
# sandwich-robust standard errors
fit_mslom_optim <- function(data,
                            start  = c(b0 = 0, b1 = 0, b2 = 0, b3 = 0),
                            reltol = 1e-9,
                            maxit  = 500) {

  opt <- optim(par     = start,
               fn      = neg_loglik_mslom,
               gr      = function(b, data) -score_mslom(data, b), # −score = ∇(−ℓ)
               data    = data,
               method  = "BFGS",
               control = list(reltol = reltol, maxit = maxit),
               hessian = TRUE)
  
  if (opt$convergence != 0)
    warning("optim() did not converge: code = ", opt$convergence)
  
  beta_hat <- opt$par
  names(beta_hat) <- c("b0", "b1", "b2", "b3")
  
  # Compute the sandwich variance
  y <- data$Y
  w <- data$w
  b0  <- beta_hat[1];  b1 <- beta_hat[2];  b2 <- beta_hat[3];  b3 <- beta_hat[4]
  z   <- 1 + b1 * data$G + b2 * data$E + b3 * data$G * data$E
  g   <- cbind(1,
               data$G / z,
               data$E / z,
               data$G * data$E / z)
  eta <- b0 + log(z)
  p   <- plogis(eta)
  n   <- nrow(data)
  omega <- w * p * (1 - p)
  bread  <- t(g) %*% (g * omega) / n
  meat   <- t(g) %*% (g * (w * (y - p))^2) / n
  V_sand <- solve(bread) %*% meat %*% solve(bread) / n
  se_robust <- sqrt(diag(V_sand))
  
  # Original se. Note that we used the negative likelihood here.
  V_hess <- solve(opt$hessian)
  se_naive <- sqrt(diag(V_hess))
  
  list(coefficients = beta_hat,
       se = se_robust,
       vcov_robust = V_sand,
       se_naive = se_naive,
       vcov_naive = V_hess,
       optim_details = opt)
}


################################################################################
# library(geex)
# 
# estFUN_mslom <- function(data) {           # <-- name must be "data"
#   with(data, {
#     function(theta) {
#       b0 <- theta[1]; b1 <- theta[2]; b2 <- theta[3]; b3 <- theta[4]
#       
#       z   <- 1 + b1*G + b2*E + b3*G*E
#       eta <- b0 + log(z)
#       p   <- plogis(eta)
#       S   <- w * (Y - p)
#       
#       c(
#         S * 1,            # dℓ/db0
#         S * G / z,        # dℓ/db1
#         S * E / z,        # dℓ/db2
#         S * G*E / z       # dℓ/db3
#       )
#     }
#   })
# }

# beta_hat <- fit_mslom_optim(dat_w)$coefficients
# 
# fit_geex <- m_estimate(
#   estFUN         = estFUN_mslom,
#   data           = dat_w,
#   roots          = beta_hat,      # skip root-finding …
#   compute_roots  = FALSE          # … and go straight to the sandwich
# )
# 
# coef_geex  <- coef(fit_geex)               # point estimates
# vcov_geex  <- vcov(fit_geex)               # sandwich
# se_geex    <- sqrt(diag(vcov_geex))
