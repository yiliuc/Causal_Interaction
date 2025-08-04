################################################################################
simulate_data_intercept <- function(n, beta0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Covariates
  X <- matrix(rnorm(n * 6), ncol = 6,
              dimnames = list(NULL, paste0("X", 1:6)))
  X1 <- X[, 1]; X2 <- X[, 2]; X3 <- X[, 3]
  X4 <- X[, 4]; X5 <- X[, 5]; X6 <- X[, 6]
  
  # marginal propensities
  eta_G <-  0.3*X1 + 0.4*X2 + 0.5*X5 + 0.5*X6
  pi_G  <- plogis(eta_G)
  
  # eta_E <- -0.5*X3 + 0.2*X4 + 0.3*X5 + 0.6*X6
  eta_E <- -0.4*X3 + 0.3*X4 + 0.5*X5 + 0.5*X6
  pi_E  <- plogis(eta_E)
  
  # joint probabilities & multinomial draw
  p00 <- (1 - pi_G) * (1 - pi_E)
  p10 <-  pi_G      * (1 - pi_E)
  p01 <- (1 - pi_G) *  pi_E
  p11 <-  pi_G      *  pi_E
  
  probs <- cbind(p00, p10, p01, p11)
  draw <- apply(probs, 1, function(p) {
    sample(1:4, size = 1, prob = p)
  })
  
  G <- as.integer(draw %in% c(2, 4))  # (1,0) and (1,1)
  E <- as.integer(draw %in% c(3, 4))  # (0,1) and (1,1)
  
  # Outcome
  eta_Y <- beta0 +
    # 0.3*X1 - 0.5*X2 + 0.7*X3 + 0.9*X4 - 0.2*X5 + 0.3*X6 +
    0.4*X1 - 0.5*X2 + 0.7*X3 + 0.9*X4 - 0.2 * X5 + 0.4 * X6 +
    # 0.3*G + 0.4*E + 0.8*(G*E)
    0.3*G + 0.4*E + 0.8*(G*E)
  pi_Y <- plogis(eta_Y)
  Y    <- rbinom(n, 1, pi_Y)
  
  return(data.frame(Y, G, E, X))
}
