simulate_data_intercept <- function(n, beta0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  ## ------------------------------------------------------------------
  ## 1. Generate covariates  X1 … X8   ~  N(0,1), mutually independent
  ## ------------------------------------------------------------------
  X <- matrix(rnorm(n * 8), ncol = 8,
              dimnames = list(NULL, paste0("X", 1:8)))
  X1 <- X[, 1];  X2 <- X[, 2]        #   X_G  (instrumental for G)
  X3 <- X[, 3];  X4 <- X[, 4]        #   X_E  (instrumental for E)
  X5 <- X[, 5];  X6 <- X[, 6]        #   X_Y  (outcome-only)
  X7 <- X[, 7];  X8 <- X[, 8]        #   C    (common confounders)
  
  ## ------------------------------------------------------------------
  ## 2. Propensity scores for each treatment, given its instruments + C
  ## ------------------------------------------------------------------
  eta_G <-  0.8 * X1 + 0.5 * X2 + 0.3 * X7 + 0.5 * X8      # f(Z_G)
  pi_G  <- plogis(eta_G)
  
  eta_E <- -0.5 * X3 + 0.2 * X4 + 0.3 * X7 + 0.5 * X8      # f(Z_E)
  pi_E  <- plogis(eta_E)
  
  ## ------------------------------------------------------------------
  ## 3. Joint treatment assignment  (G,E)  ~ Multinomial(1, θ00, θ10, θ01, θ11)
  ## ------------------------------------------------------------------
  p00 <- (1 - pi_G) * (1 - pi_E)
  p10 <-  pi_G      * (1 - pi_E)
  p01 <- (1 - pi_G) *  pi_E
  p11 <-  pi_G      *  pi_E
  
  probs <- cbind(p00, p10, p01, p11)
  draw <- apply(probs, 1, function(p) {
    sample(1:4, size = 1, prob = p)
  })
  
  G <- as.integer(draw %in% c(2, 4))   # (1,0) or (1,1)
  E <- as.integer(draw %in% c(3, 4))   # (0,1) or (1,1)
  
  ## ------------------------------------------------------------------
  ## 4. Outcome model  Pr(Y = 1 | Z_Y, G, E)
  ## ------------------------------------------------------------------
  eta_Y <- beta0 +
    # 0.3 * X5  + 0.5 * X6  - 0.2 * X7 - 0.1 * X8 +   # outcome-only + C
    0.8 * X5  + 0.5 * X6  - 0.8 * X7 + 0.6 * X8 +   # outcome-only + C
    0.3 * G   + 0.4 * E  + 0.8 * (G * E)            # treatments & interaction
  pi_Y <- plogis(eta_Y)
  Y    <- rbinom(n, size = 1, prob = pi_Y)
  
  ## ------------------------------------------------------------------
  ## 5. Return tidy data frame
  ## ------------------------------------------------------------------
  return(data.frame(Y, G, E, X))
}
