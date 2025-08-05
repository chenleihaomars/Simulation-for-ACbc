library(MASS)   


#'
#' @param n     Number of observations to draw
#' @param d     Dimension of X
#' @param rho   Correlation between each X_j and Y
#' @return A list with
#'    X: nd matrix of X observations  
#'    Y: length?\n vector of Y observations
Simulate_data <- function(n, d, rho) {
  # total dimension
  p <- d + 1

  # 1) build the AR-1 covariance matrix: Sigma[i,j] = rho^|i-j|
  idx   <- seq_len(p)
  Sigma <- rho ^ abs(outer(idx, idx, "-"))

  # 2) check PD (AR-1 is PD if |rho| < 1)
  ev <- eigen(Sigma, symmetric=TRUE, only.values=TRUE)$values
  if (min(ev) <= 0)
    warning("Covariance not PD (min eigen=", min(ev), ")")

  # 3) simulate n draws from N(0, Sigma)
  dat <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  X <- dat[, 1:d, drop = FALSE]
  X <- pnorm(X)
  Y <- dat[, p]
  # 4) split into X (first d columns) and Y (last column)
  list(
    X=X,
    Y=Y
  )
}

Simulate_data_1 <- function(n, d, rho) {
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


# Simulate_data <- function(n, d, rho, sigma_noise = 0.5) {
#   Sigma_X <- rho ^ abs(outer(1:d, 1:d, "-"))
#   X <- pnorm(mvrnorm(n, mu = rep(0, d), Sigma = Sigma_X))
# 
#   # Explicit nonlinear function, e.g. quadratic + interaction + sin terms
#   nonlinear_mean <- rowSums(X^2) + rowSums(X[,1:(d-1)] * X[,2:d]) + sin(rowSums(X))
# 
#   # Adding Gaussian noise
#   Y <- nonlinear_mean + rnorm(n, mean = 0, sd = sigma_noise)
# 
#   list(X = X, Y = Y)
# }




#' @param X A numeric vector of length n, or an n??d numeric matrix / data.frame.
#' @param Y A numeric vector of length n.
#' @return A single number: the estimate T_n.
#' Compute Chatterjee's rank correlation 

Chatterjee_coefficient <- function(X, Y) {
  # Coerce X to matrix
  Xmat <- if(is.data.frame(X)) as.matrix(X) else if(is.vector(X)) matrix(X, ncol=1) else as.matrix(X)
  n <- nrow(Xmat)
  if(length(Y)!=n) stop("Length of Y must match number of rows in X")
  
  # 1) Compute ranks of Y with ties broken at random
  R <- rank(Y, ties.method = "random")
  
  # 2) Compute nearest\neighbor indices in X
  #    (for each i, find j != i minimizing ||X[i,] - X[j,]||_2)
  #    Using 'dist' then converting to full matrix for simplicity
  D <- as.matrix(dist(Xmat, method = "euclidean"))
  diag(D) <- Inf                 # ignore self-distance
  Nn <- apply(D, 1, which.min)   # Nn[i] is index of i's nearest neighbor
  
  # 3) Compute sum of min{ R[i], R[Nn[i]] }
  sum_min <- sum(pmin(R, R[Nn]))
  
  # 4) Plug into T_n formula
  T_hat <- (6 / (n^2 - 1)) * sum_min - ((2*n + 1) / (n - 1))
  
  return(T_hat)
}

# Three data splitting schemes
make_splits <- function(m, B) {
  # how to cut into three parts
  m1 <- floor(m / 3)
  m2 <- floor(2 * m / 3)
  
  # pre-allocate list of length B
  splits <- vector("list", B)
  
  for (b in seq_len(B)) {
    # random permutation of 1:m
    perm <- sample.int(m)
    
    # assign the three folds
    splits[[b]] <- list(
      I1 = perm[1:m1],
      I2 = perm[(m1+1):m2],
      I3 = perm[(m2+1):m]
    )
  }
  
  splits
}


# Build the quadratic expansion formula for d covariates named X1...Xd
make_quad_formula <- function(d) {
  # main + interactions: .^2
  # squares: I(Xj^2)
  # No intercept
  sq_terms <- paste0("I(X", 1:d, "^2)", collapse = " + ")
  as.formula(paste("~ 0 + .^2 +", sq_terms))
}

make_cubic_formula <- function(d) {
  # explicit squares
  sq_terms   <- paste0("I(X", 1:d, "^2)", collapse = " + ")
  # explicit cubes
  cube_terms <- paste0("I(X", 1:d, "^3)", collapse = " + ")
  # .^3 gives all main effects, all two-way interactions, and all three-way interactions
  as.formula(paste(
    "~ 0 + .^3 +", 
    sq_terms, "+", 
    cube_terms
  ))
}


make_fold_df <- function(X, Y) {
  df <- as.data.frame(X)
  # Ensure standard naming convention for model.matrix
  if (ncol(df) > 0) colnames(df) <- paste0("X", 1:ncol(X))
  df$Y <- Y
  return(df)
}

# Fits ridge:  beta = (X'X + n*lambda*I)^{-1} X' y
ridge_solve <- function(X, y, lambda) {
  # X: n×p design matrix
  # y: length‐n response
  n_coef <- ncol(X)
  XtX     <- crossprod(X, X)               # p×p
  A       <- XtX + nrow(X) * lambda * diag(1, n_coef)
  XtY     <- crossprod(X, y)               # p×1
  beta    <- solve(A, XtY)
  return(beta)
}

Tval <- function(rho){
  return(1 - 6 * asin(sqrt(1 - rho^2)/2)/pi)  
}