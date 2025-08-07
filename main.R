source('functions.R')
library(parallel) 


run_sim <- function(n, d, rho, z, lambda = n^(-0.8))), Bboot = 10000, S = 1000) {
  m      <- floor(sqrt(n))
  fml    <- make_quad_formula(d)
  
  # storage
  T_n_orig <- numeric(S)
  T_n_bc   <- numeric(S)
  T_bootstrap    <- numeric(Bboot)
  cover_orig <- logical(S)
  cover_bc   <- logical(S)
  Lhat <- 0
  ## --- variance via m-out-of-n bootstrap ---
  sim0 <- Simulate_data(n = n, d = d, rho = rho)
  X0   <- sim0$X; Y0 <- sim0$Y
  T_bootstrap <- replicate(Bboot, {
    idx <- sample.int(n, m, replace = FALSE)
    codec(Y0[idx], X0[idx, , drop = FALSE])
  })
  var_bootstrap <- m * mean( (T_bootstrap - mean(T_bootstrap))^2 )
  sigma_n <- sqrt(var_bootstrap)
  
  ## --- "truth" from the closed form for the design ---
  T_true <- Tval(rho)
  
  ## --- Iteration loop ---
  for (s in seq_len(S)) {
    sim <- Simulate_data(n = n, d = d, rho = rho)
    X   <- sim$X; Y <- sim$Y
    
    # original statistic
    T_n_orig[s] <- codec(Y, X)
    
    # full-sample Lhat, then bias-correct
    Lhat <- compute_Lhat(X, Y, fml, lambda)
    T_n_bc[s] <- T_n_orig[s] - 6 * Lhat
  }
  
  # CIs & coverage
  ci_lo_orig <- T_n_orig - z * sigma_n / sqrt(n)
  ci_hi_orig <- T_n_orig + z * sigma_n / sqrt(n)
  cover_orig <- (T_true >= ci_lo_orig) & (T_true <= ci_hi_orig)
  
  ci_lo_bc <- T_n_bc - z * sigma_n / sqrt(n)
  ci_hi_bc <- T_n_bc + z * sigma_n / sqrt(n)
  cover_bc <- (T_true >= ci_lo_bc) & (T_true <= ci_hi_bc)
  
  data.frame(
    n = n, d = d, rho = rho, z = z,
    RMSE_orig = sqrt(mean((T_n_orig - T_true)^2)),
    RMSE_bc   = sqrt(mean((T_n_bc   - T_true)^2)),
    ECP_orig  = mean(cover_orig),
    ECP_bc    = mean(cover_bc),
    MAD_orig  = mad(abs(T_n_orig - T_true)),
    MAD_bc    = mad(abs(T_n_bc   - T_true))
  )
}

grid <- expand.grid(
  n   = c(300,600,900),
  rho = c(0,0.3,0.5,0.7,0.9),
  d   = c(6, 8, 10),
  z   = c(1.96),
  stringsAsFactors = FALSE
)

# 1) Decide how many workers to spawn:
ncores <- min(detectCores() - 1, nrow(grid))

# 2) Create a PSOCK cluster:
cl <- makeCluster(ncores)

# 3) Export all the objects and libraries your worker needs:
clusterExport(cl, varlist = c("run_sim", "make_quad_formula", "Simulate_data", "grid","compute_Lhat", "Tval"))
clusterEvalQ(cl, {
  library(MASS)
  library(FNN)
  library(FOCI)
})

# 4) Use parLapply:
res_list <- parLapply(cl, seq_len(nrow(grid)), function(i) {
  params <- grid[i, ]
  with(params, run_sim(n, d, rho, z))
})

# 5) Shut down the cluster:
stopCluster(cl)

# 6) Combine results:
results <- do.call(rbind, res_list)

write.csv(results, file = "results.csv", row.names = FALSE)
