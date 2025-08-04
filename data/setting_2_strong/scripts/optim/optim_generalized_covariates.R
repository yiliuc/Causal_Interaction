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
# 2. Compute the negative log-likelihood
neg_loglik_kcov_optim <- function(par, data, covars = character()) {
  
  lst <- linodds_check_extract(data, covars)
  y <- lst$Y;  g <- lst$G;  e <- lst$E;  X <- lst$X     # X is n × k
  k <- ncol(X)
  
  ## unpack ----------------------------------------------------
  b0 <- par[1]; b1 <- par[2]; b2 <- par[3]; b3 <- par[4]
  gam <- if (k) par[5:(4 + k)] else numeric(0)
  
  z <- 1 + b1 * g + b2 * e + b3 * g * e
  if (any(z <= 0))
    return(.Machine$double.xmax)
  
  eta <- b0 + if (k) drop(X %*% gam) else 0
  eta <- eta + log(z)
  
  -sum(y * eta - log1p(exp(eta)))
}

################################################################################
# 3. Compute the score
score_kcov_optim <- function(data, theta, covars = character()){
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
    sum(r * G / z),
    sum(r * E / z),
    sum(r * G * E / z)
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

# score_kcov_optim <- function(par, data, covars = character()) {
#   
#   ## fast invalid-region check to keep optim() stable
#   lst <- linodds_check_extract(data, covars)
#   g <- lst$G;  e <- lst$E
#   b1 <- par[2]; b2 <- par[3]; b3 <- par[4]
#   z <- 1 + b1 * g + b2 * e + b3 * g * e
#   if (any(z <= 0))
#     return(rep(NA_real_, length(par)))   # forces optim() to back-track
#   
#   ## call the algebra you already coded
#   score_kcov(data, par, covars)
# }

################################################################################
# 3. Optimizer
fit_kcovariates_optim <- function(data,
                              covars  = character(),
                              start   = rep(0, 4 + length(covars)),
                              method  = "BFGS",
                              reltol  = 1e-9,
                              maxit   = 1000) {
  
  k <- length(covars)
  if (length(start) != 4 + k)
    stop("'start' must have length 4 + length(covars)")
  
  nll  <- function(par) neg_loglik_kcov_optim(par, data, covars)
  grad <- function(par) -score_kcov_optim(data, par, covars)  # minus ⇒ ∇(−ℓ)
  
  opt <- optim(par       = start,
               fn        = nll,
               gr        = grad,
               method    = method,
               control   = list(reltol = reltol, maxit = maxit),
               hessian   = TRUE)
  
  if (opt$convergence != 0)
    warning("optim() did not converge (code = ", opt$convergence, ")")
  
  par_names <- c(paste0("b", 0:3),
                 if (k) paste0("g_", covars))
  coef <- opt$par
  names(coef) <- par_names
  
  vcov <- tryCatch(solve(opt$hessian),
                   error = function(e) matrix(NA_real_, length(coef), length(coef)))
  dimnames(vcov) <- list(par_names, par_names)
  se <- sqrt(diag(vcov))
  
  structure(list(coefficients = coef,
                 se           = se,
                 vcov         = vcov,
                 optim        = opt),
            class = "linodds_fit")
}

# fit <- fit_kcovariates_optim(sim_dat,
#                          covars = c("X5","X6", "X7", "X8"),
#                          start  = rep(0, 8))
# fit$coefficients
# fit$se

