################################################################################
# 1. Compute the negative log-likelihood
neg_loglik_two_covariates <- function(par, dat, covars) {
  y <- dat$Y;  g <- dat$G;  e <- dat$E
  x1 <- dat[[covars[1]]];  x2 <- dat[[covars[2]]]
  
  b0 <- par[1]; b1 <- par[2]; b2 <- par[3]; b3 <- par[4]
  g1 <- par[5]; g2 <- par[6]
  
  z <- 1 + b1*g + b2*e + b3*g*e
  if (any(z <= 0))
    return(.Machine$double.xmax)
  
  eta <- b0 + g1*x1 + g2*x2 + log(z)
  -sum(y*eta - log1p(exp(eta)))
}

################################################################################
# 2. Compute the score
score_two_covariates <- function(par, dat, covars) {
  y <- dat$Y;  g <- dat$G;  e <- dat$E
  x1 <- dat[[covars[1]]];  x2 <- dat[[covars[2]]]

  b0 <- par[1]; b1 <- par[2]; b2 <- par[3]; b3 <- par[4]
  g1 <- par[5]; g2 <- par[6]
  
  z <- 1 + b1*g + b2*e + b3*g*e
  if (any(z <= 0))
    return(rep(NA_real_, 6))
  
  Odds <- exp(b0 + g1*x1 + g2*x2) * z
  p    <- Odds / (1 + Odds)
  r    <- y - p
  
  c(sum(r),
    sum(r * g / z),
    sum(r * e / z),
    sum(r * g*e / z),
    sum(r * x1),
    sum(r * x2))
}

################################################################################
# 3. Optimizer
fit_two_covariates_optim <- function(dat, covars = c("X1","X2"),
                            start = rep(0, 6),
                            reltol = 1e-9, maxit = 500){
  ## capture data & covars in closures so optim() only sees 'par'
  nll  <- function(par) neg_loglik_two_covariates(par, dat, covars)
  grad <- function(par) -score_two_covariates  (par, dat, covars)   # minus ⇒ ∇(−ℓ)
  
  opt <- optim(par       = start,
               fn        = nll,
               gr        = grad,
               method    = "BFGS",
               control   = list(reltol = reltol, maxit = maxit),
               hessian   = TRUE)
  
  if (opt$convergence != 0)
    warning("optim() did not converge (code = ", opt$convergence, ")")
  
  par_names <- c("b0", "b1", "b2", "b3",
                 paste0("g_", covars[1]),
                 paste0("g_", covars[2]))
  
  coef <- opt$par
  names(coef) <- par_names
  vcov <- solve(opt$hessian)
  dimnames(vcov) <- list(par_names, par_names)
  se   <- sqrt(diag(vcov))
  
  structure(list(coefficients = coef,
                 se           = se,
                 vcov         = vcov,
                 optim        = opt),
            class = "linodds_fit")
}

# fit <- fit_two_covariates_optim(sim_dat, covars = c("X3","X4"),
#                        start  = c(0,0,0,0,0,0))
# 
# fit$coefficients
# fit$se

# glm(Y ~ G*E + X1 + X2, family = binomial, data = sim_dat)$coeff
# compute_RERI_OR(0.9727883, 0.6174996, 0.3237160)







