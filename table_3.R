#!/usr/bin/env Rscript
# ======================================================== #
#
#                  Replication of Table 3 
#
#                 Gabriel E. Cabrera-Guzmán
#                The University of Manchester
#
#                        Spring, 2026
#
#                https://gcabrerag.rbind.io
#
# ------------------------------ #
# email: gabriel.cabreraguzman@postgrad.manchester.ac.uk
# ======================================================== #

# Replication of Table 3: In-sample predictive regression results
# RRZ (2016), "Short interest and aggregate stock returns", JFE
#
# r_{t:t+h} = alpha + beta * x_t + epsilon_{t:t+h}
# where r_{t:t+h} = (1/h)(r_{t+1} + ... + r_{t+h})
# Newey-West HAC t-statistics with h lags
# Wild bootstrapped p-values (1000 iterations)

source("load.R")

# ==========================================
#      In-sample predictive regressions
# ------------------------------------------

n_pred <- ncol(GW) + 1 # 14 GW predictors + SII
beta_hat <- array(NA, dim = c(n_pred, 4, length(h_vec)))
beta_hat_PC_SII <- array(NA, dim = c(1, 4, length(h_vec)))

for (j in seq_along(h_vec)) {
  
  hh <- h_vec[j]
  n_obs <- T_full - hh  # usable observations

  y_j <- 100 * r_h[2:(T_full - (hh - 1)), j]

  # --- Individual GW predictors ---
  for (i in 1:ncol(GW)) {
    
    X_ij <- cbind(1, GW_standardize[1:(T_full - hh), i])
    res <- nwest_ols(y_j, X_ij, hh)
    beta_hat[i, , j] <- c(res$beta[2], res$tstat[2], NA, 100 * res$rsqr)
    
  }

  # --- SII (negated, so positive beta = higher SII predicts lower returns) ---
  X_sii <- cbind(1, -SII[1:(T_full - hh)])
  res_sii <- nwest_ols(y_j, X_sii, hh)
  beta_hat[n_pred, , j] <- c(res_sii$beta[2], res_sii$tstat[2], NA, 100 * res_sii$rsqr)

  # --- Multiple regression: PC(1:3) + SII ---
  # Note: MATLAB uses r_h(2:end-(h-1),j) with predictors at 1:end-1-(h-1)
  # This aligns returns at t+1:t+h with predictors at t
  n_obs_pc <- T_full - 1 - (hh - 1)
  y_pc <- 100 * r_h[2:(T_full - (hh - 1)), j]

  X_PC_j <- cbind(1, PC_GW[1:n_obs_pc,])
  X_PC_SII_j <- cbind(1, PC_GW[1:n_obs_pc,], -SII[1:n_obs_pc])

  res_PC <- nwest_ols(y_pc, X_PC_j, hh)
  res_PC_SII <- nwest_ols(y_pc, X_PC_SII_j, hh)

  SSE_reduced <- sum(res_PC$resid ^ 2)
  SSE_full <- sum(res_PC_SII$resid ^ 2)
  partial_rsqr <- (SSE_reduced - SSE_full) / SSE_reduced

  beta_hat_PC_SII[1, , j] <- c(res_PC_SII$beta[length(res_PC_SII$beta)], res_PC_SII$tstat[length(res_PC_SII$tstat)], NA, 100 * partial_rsqr)

}

# ==========================================
# Wild bootstrap p-values (fixed-regressor)
# ------------------------------------------

cat("Computing wild bootstrap p-values (B=1000)...\n")

# Kitchen-sink regression for bootstrap DGP (monthly horizon)
X_sink <- cbind(GW_standardize, -SII)
X_sink <- X_sink[, -c(4, 11)]  # remove DE and TMS
X_sink <- cbind(1, X_sink)
res_sink <- ols_fit(r_h[2:T_full, 1], X_sink[1:(T_full - 1), ])
epsilon_hat <- res_sink$resid

B <- 1000
set.seed(42)  # reproducibility (MATLAB uses rng('default'))
beta_hat_tstat_star <- array(NA, dim = c(B, n_pred, length(h_vec)))
beta_hat_PC_SII_tstat_star <- matrix(NA, B, length(h_vec))

for (b in 1:B) {
  
  if (b %% 100 == 0) cat("  Bootstrap iteration", b, "/", B, "\n")

  u_star <- rnorm(T_full - 1)
  r_star <- c(r[1], mean(r) + epsilon_hat * u_star)

  # Cumulative returns for bootstrap sample
  r_h_star <- matrix(NA, T_full, length(h_vec))
  
  for (j in seq_along(h_vec)) {
    
    hh <- h_vec[j]
    
    for (tt in 1:(T_full - (hh - 1))) {
      
      r_h_star[tt, j] <- mean(r_star[tt:(tt + hh - 1)])
      
    }
    
  }

  for (j in seq_along(h_vec)) {
    
    hh <- h_vec[j]
    n_obs_b <- T_full - 1 - (hh - 1)
    y_star_j <- 100 * r_h_star[2:(T_full - (hh - 1)), j]

    # Individual predictors
    for (i in 1:ncol(GW)) {
      
      X_ij <- cbind(1, GW_standardize[1:n_obs_b, i])
      res_b <- nwest_ols(y_star_j, X_ij, hh)
      beta_hat_tstat_star[b, i, j] <- res_b$tstat[2]
      
    }

    # SII
    X_sii <- cbind(1, -SII[1:n_obs_b])
    res_b <- nwest_ols(y_star_j, X_sii, hh)
    beta_hat_tstat_star[b, n_pred, j] <- res_b$tstat[2]

    # PC + SII
    X_PC_SII_j <- cbind(1, PC_GW[1:n_obs_b, ], -SII[1:n_obs_b])
    res_b <- nwest_ols(y_star_j, X_PC_SII_j, hh)
    beta_hat_PC_SII_tstat_star[b, j] <- res_b$tstat[length(res_b$tstat)]
    
  }
  
}

# Compute bootstrap p-values
for (j in seq_along(h_vec)) {
  
  for (i in 1:n_pred) {
    
    beta_hat[i, 3, j] <- sum(beta_hat_tstat_star[, i, j] > beta_hat[i, 2, j]) / B
    
  }
  
  beta_hat_PC_SII[1, 3, j] <- sum(beta_hat_PC_SII_tstat_star[, j] > beta_hat_PC_SII[1, 2, j]) / B
  
}

# ==========================================
#               Display results
# ------------------------------------------

cat("\n=== TABLE 3: In-sample predictive regression results ===\n")
cat("Sample: 1973:01-2014:12\n\n")

for (j in seq_along(h_vec)) {
  
  hh <- h_vec[j]
  cat(sprintf("--- h = %d ---\n", hh))
  cat(sprintf("%-6s %8s %8s %8s %8s\n", "Pred", "beta", "t-stat", "p-val", "R2(%)"))
  cat(strrep("-", 42), "\n")
  
  for (i in 1:n_pred) {
    
    pval_str <- ifelse(beta_hat[i, 3, j] < 0.005, "<0.005", sprintf("%.2f", beta_hat[i, 3, j]))
    stars <- ""
    if (!is.na(beta_hat[i, 3, j])) {
      
      if (beta_hat[i, 3, j] <= 0.01) stars <- "***"
      else if (beta_hat[i, 3, j] <= 0.05) stars <- "**"
      else if (beta_hat[i, 3, j] <= 0.10) stars <- "*"
      
    }
    
    cat(sprintf("%-6s %8.2f %8.2f%s %6s %8.2f\n",
                pred_names[i], beta_hat[i, 1, j],
                beta_hat[i, 2, j], stars, pval_str, beta_hat[i, 4, j]))
    
  }
  
  cat(sprintf("%-6s %8.2f %8.2f%s %6s %8.2f  (partial R2)\n",
              "SII(-)|PC",
              beta_hat_PC_SII[1, 1, j], beta_hat_PC_SII[1, 2, j],
              ifelse(!is.na(beta_hat_PC_SII[1, 3, j]) &&
                       beta_hat_PC_SII[1, 3, j] <= 0.01, "***",
                     ifelse(!is.na(beta_hat_PC_SII[1, 3, j]) &&
                              beta_hat_PC_SII[1, 3, j] <= 0.05, "**",
                            ifelse(!is.na(beta_hat_PC_SII[1, 3, j]) &&
                                     beta_hat_PC_SII[1, 3, j] <= 0.10, "*", ""))),
              ifelse(beta_hat_PC_SII[1, 3, j] < 0.005, "<0.005",
                     sprintf("%.2f", beta_hat_PC_SII[1, 3, j])),
              beta_hat_PC_SII[1, 4, j]))
  cat("\n")
  
}
