################################################################################
# 1. Negative likelihood
neg_loglik_kcov_w <- function(par, data, covars = character()) {
  
  lst <- linodds_check_extract(data, covars)
  y <- lst$Y;  g <- lst$G;  e <- lst$E;  X <- lst$X
  w <- data$w                                    # <── weights already added
  k <- ncol(X)
  
  ## unpack
  b0 <- par[1]; b1 <- par[2]; b2 <- par[3]; b3 <- par[4]
  gam <- if (k) par[5:(4 + k)] else numeric(0)
  
  z <- 1 + b1*g + b2*e + b3*g*e
  if (any(z <= 0))                   # invalid region ⇒ huge penalty
    return(.Machine$double.xmax)
  
  eta  <- b0 + if (k) drop(X %*% gam) else 0
  eta  <- eta + log(z)               # log-odds
  -sum( w * (y * eta - log1p(exp(eta))) )
}

################################################################################
# 2. Score function
score_kcov_w <- function(data, theta, covars = character()) {
  
  k <- length(covars)
  lst <- linodds_check_extract(data, covars)
  Y <- lst$Y;  G <- lst$G;  E <- lst$E;  X <- lst$X
  w <- data$w
  
  b0 <- theta[1]; b1 <- theta[2]; b2 <- theta[3]; b3 <- theta[4]
  g  <- if (k) theta[5:(4 + k)] else numeric(0)
  
  z <- 1 + b1*G + b2*E + b3*G*E
  O <- exp(b0 + if (k) X %*% g else 0) * z
  p <- O / (1 + O)
  r <- as.vector(Y - p)
  
  wr <- w * r
  
  S_beta <- c(
    sum(wr),
    sum(wr * G      / z),
    sum(wr * E      / z),
    sum(wr * G * E  / z)
  )
  
  if (k) {
    S_gamma <- as.vector(colSums(wr * X))
    names(S_gamma) <- paste0("g_", covars)
    out <- c(S_beta, S_gamma)
  } else {
    out <- S_beta
  }
  names(out)[1:4] <- paste0("b", 0:3)
  out
}

################################################################################
# 3. Optimiser
fit_kcov_w_optim <- function(data,
                             covars  = character(),
                             start   = rep(0, 4 + length(covars)),
                             method  = "BFGS",
                             reltol  = 1e-9,
                             maxit   = 1000) {
  
  k <- length(covars)
  if (length(start) != 4 + k)
    stop("'start' must have length 4 + length(covars)")
  if (!"w" %in% names(data))
    stop("column 'w' (weights) not found — run compute_ipw() first")
  
  nll  <- function(par) neg_loglik_kcov_w(par, data, covars)
  grad <- function(par) -score_kcov_w(data, par, covars)
  
  opt <- optim(par = start,
               fn = nll,
               gr = grad,
               method = method,
               control = list(reltol = reltol, maxit = maxit),
               hessian = TRUE)
  
  if (opt$convergence != 0)
    warning("optim() did not converge (code = ", opt$convergence, ")")
  
  par_names <- c(paste0("b", 0:3),
                 if (k) paste0("g_", covars))
  coef <- opt$par
  names(coef) <- par_names
  
  vcov <- tryCatch(solve(opt$hessian),
                   error = function(e)
                     matrix(NA_real_, length(coef), length(coef)))
  dimnames(vcov) <- list(par_names, par_names)
  
  structure(list(coefficients = coef,
                 se           = sqrt(diag(vcov)),
                 vcov         = vcov,
                 optim        = opt),
            class = "linodds_fit")
}

############################################################
## Example work-flow
############################################################
## 1.  Add joint-IPW weights
# set.seed(2025)
# sim_dat <- simulate_data(3000)
# dat_w <- compute_ipw(sim_dat, covars = paste0("X", 1:6))
# 
# ## 2.  Fit weighted conditional linear-odds model
# fit_wclom <- fit_kcov_w_optim(dat_w, covars = paste0("X", 1:6))
# fit_wclom$coefficients
# fit_wclom$se
