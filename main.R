source('functions.R')
library(parallel)  

n<-300


# 1) Parameters 
run_sim <- function(n, d, rho, z, lambda=n^(-1/3),
                    B = 3, Bboot = 10000, S = 1000) {
  m            <- floor(sqrt(n))
  full_f       <- make_quad_formula(d)
  splits       <- make_splits(n, B)
  
  
  # 2) Storage
  T_n_orig     <- numeric(S)
  T_n_bc       <- numeric(S)
  T_bootstrap    <- numeric(Bboot)
  cover_orig <- logical(S)
  cover_bc   <- logical(S)
  
  # 3) Estimate variance using m-out-of-n bootstrap
  
  sim <- Simulate_data_1(n = n, d = d, rho = rho)
  X   <- sim$X; Y <- sim$Y
  
  T_bootstrap <- replicate(Bboot, {
    idx <- sample.int(n, m, replace = FALSE)
    Chatterjee_coefficient(X[idx, , drop = FALSE], Y[idx])
  })
  var_bootstrap <- m * mean((T_bootstrap - mean(T_bootstrap))^2) 
  sigma_n <- sqrt(var_bootstrap)
  
  # 4) True T
  T_true <- Tval(rho)
  
  # 5) Estimate T_n_orig and T_n_bc
  split_n <- make_splits(n,B)
  full_formula <- make_quad_formula(d)
  for(s in seq_len(S)){
    
    sim <- Simulate_data_1(n = n, d = d, rho = rho)
    X   <- sim$X; Y <- sim$Y
    T_n_orig[s] <- Chatterjee_coefficient(X, Y)
    L <- numeric(B)
    for(b in seq_len(B)) {
      
      # 1) Split the n-subset
      
      scheme  <- split_n[[b]] # Using the b-th scheme
      pos1    <- scheme$I1
      pos2    <- scheme$I2
      pos3    <- scheme$I3
      X1s <- X[pos1, , drop = FALSE]; Y1s <- Y[pos1]
      X2s <- X[pos2, , drop = FALSE]; Y2s <- Y[pos2]
      X3s <- X[pos3, , drop = FALSE]; Y3s <- Y[pos3]
      
      # Build data.frames for each fold
      
      
      df1 <- make_fold_df(X1s, Y1s)
      df2 <- make_fold_df(X2s, Y2s)
      df3 <- make_fold_df(X3s, Y3s)
      
      # 2) Linear ridge regression 
      
      # Find nearest neighbors within fold 2
      d2 <- as.matrix(dist(df2[, paste0("X", 1:(ncol(df2)-1))]))
      diag(d2) <- Inf
      neigh_idx <- apply(d2, 1, which.min)
      
      
      n2 <- nrow(df2)
      n3 <- nrow(df3)
      total <- 0
      
      # Build design matrices
      Z_tr   <- model.matrix(full_formula, data = df1[, 1:d, drop=FALSE])
      Z_tr_sc   <- Z_tr
      # ctr       <- attr(Z_tr_sc,"scaled:center") # Here you can choose perform scaling or not
      # scl       <- attr(Z_tr_sc,"scaled:scale")
      Z_pr   <- model.matrix(full_formula, data = df2[, 1:d, drop=FALSE])
      # Z_pr_sc   <- scale(Z_pr,    center=ctr, scale=scl)
      Z_pr_sc   <- Z_pr
      
      Z_tr_sc  <- cbind(1, Z_tr_sc)                        # now first col = intercept
      Z_pr_sc  <- cbind(1, Z_pr_sc)
      for (j in seq_len(n3)) {
        t    <- df3$Y[j]
        # Binary response indicator at threshold t
        Z_resp <- as.integer(df1$Y >= t)
        # Model fitting
        beta    <- ridge_solve(Z_tr_sc, Z_resp, lambda)
        # Return predicted conditional means
        G_i     <- as.vector(Z_pr_sc   %*% beta)
        G_nn    <- G_i[neigh_idx]
        total <- total + sum(G_i * G_nn - G_i^2)
      }
      Tb <- total / (n2 * n3)
      L[b] <- Tb
    }
    T_n_bc[s] <- T_n_orig[s]-6*mean(L)
  }
  
  # Compute RMSE and ECP for T_n_orig and T_bc
  
  ci_lo_orig <- T_n_orig - z * sigma_n / sqrt(n) 
  ci_hi_orig <- T_n_orig + z * sigma_n / sqrt(n)
  cover_orig <- (T_true >= ci_lo_orig) & (T_true <= ci_hi_orig)
  
  ci_lo_bc <- T_n_bc - z * sigma_n / sqrt(n) 
  ci_hi_bc <- T_n_bc + z * sigma_n / sqrt(n)
  cover_bc <- (T_true >= ci_lo_bc) & (T_true <= ci_hi_bc)
  
  data.frame(
    n = n, d = d, rho = rho, z = z,
    RMSE_orig = sqrt(mean((T_n_orig - T_true)^2)),
    RMSE_bc   = sqrt(mean((T_n_bc - T_true)^2)),
    ECP_orig  = mean(cover_orig),
    ECP_bc    = mean(cover_bc),
    MAD_orig = mad(abs(T_n_orig - T_true)),
    MAD_bc = mad(abs(T_n_bc - T_true))
  )
}


grid <- expand.grid(
  n   = c(300,600,900),
  rho = c(0,0.3,0.5,0.7,0.9),
  d   = c(6,8,10),
  z   = c(1.96),
  stringsAsFactors = FALSE
)

# 1) Decide how many workers to spawn:
ncores <- min(detectCores() - 1, nrow(grid))

# 2) Create a PSOCK cluster:
cl <- makeCluster(ncores)

# 3) Export all the objects and libraries your worker needs:
clusterExport(cl, varlist = c("run_sim", "make_quad_formula", "make_splits",
                              "Simulate_data_1", "Chatterjee_coefficient",
                              "grid", "make_fold_df","ridge_solve","Tval"))
clusterEvalQ(cl, {
  library(MASS)
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
print(results)

# writexl::write_xlsx(results, path = "results.xlsx")



