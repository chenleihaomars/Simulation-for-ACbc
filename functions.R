library(MASS)   
library(FNN)
library(FOCI)

Simulate_data <- function(n, d, rho) {
  # Simulate X_1, ..., X_d ~ N(0,1), independently
  X <- matrix(rnorm(n * d), nrow = n, ncol = d)
  
  # Simulate independent Z ~ N(0,1)
  Z <- rnorm(n)
  
  # Define Y explicitly as linear combination of X_1 and Z
  Y <- rho * X[, 1] + sqrt(1 - rho^2) * Z
  X <- pnorm(X)
  Y <- pnorm(Y)
  return(list(X = X, Y = Y))
}

# Build the quadratic expansion formula for d covariates named X1...Xd
make_quad_formula <- function(d) {
  # main + interactions: .^2
  # squares: I(Xj^2)
  # No intercept
  sq_terms <- paste0("I(X", 1:d, "^2)", collapse = " + ")
  as.formula(paste("~ 0 + .^2 +", sq_terms))
}

# calculate value of T(Y, X) when Y and X are normally distributed with correlation rho
Tval <- function(rho){
  return(1 - 6 * asin(sqrt(1 - rho^2)/2)/pi)  
}

# nearest-neighbor based full-sample estimator of L
compute_Lhat <- function(X, Y, full_formula, lambda) {
  n <- nrow(X); d <- ncol(X)
  
  # 1) 1-NN on X (exclude self)
  nn_idx <- FNN::get.knn(as.matrix(X), k = 2)$nn.index[, 2]
  
  # 2) Design matrix (add intercept)
  Xdf <- as.data.frame(X); names(Xdf) <- paste0("X", 1:d)
  Z   <- model.matrix(full_formula, data = Xdf)
  Z   <- cbind(1, Z)
  
  # 3) Precompute cross-products and factorization
  XtX <- crossprod(Z)
  A <- XtX + diag(n * lambda, ncol(Z)) 
  R   <- chol(A)  # A = R'R 
  
  total <- 0
  for (j in seq_len(n)) {
    zresp <- as.integer(Y >= Y[j])    # length n
    XtY   <- crossprod(Z, zresp)      # length p
    # Solve (R'R) beta = XtY
    beta  <- backsolve(R, forwardsolve(t(R), XtY))
    G     <- as.vector(Z %*% beta)
    v     <- G * G[nn_idx] - G^2
    total <- total + (sum(v) - v[j])  # drop i=j
  }

  Lhat <- total / (n * (n - 1))
  return(Lhat)
}


