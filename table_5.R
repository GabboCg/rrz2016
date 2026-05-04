#!/usr/bin/env Rscript
# ======================================================== #
#
#                  Replication of Table 5 
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

# Replication of Table 5: Out-of-sample test results
# RRZ (2016), "Short interest and aggregate stock returns", JFE
#
# Out-of-sample R^2 statistics (Campbell & Thompson 2008)
# Clark-West (2007) MSFE-adjusted statistics
# Harvey-Leybourne-Newbold (1998) encompassing tests
# OOS period: 1990:01-2014:12

# Load auxiliary scripts
source("load.R")

# ==========================================
# Out-of-sample forecasts (expanding window)
# ------------------------------------------

cat("Computing out-of-sample forecasts...\n")

FC_PM <- rep(NA, P_oos)                                        # prevailing mean
FC_PR <- array(NA, dim = c(P_oos, ncol(GW) + 1, length(h_vec)))  # predictor-based

for (p in 1:P_oos) {
  
  if (p %% 50 == 0) cat("  OOS period", p, "/", P_oos, "\n")

  # Prevailing mean benchmark (uses log excess returns)
  FC_PM[p] <- mean(r[1:(R_oos + (p - 1))])

  for (j in seq_along(h_vec)) {
    
    hh <- h_vec[j]
    n_est <- R_oos + (p - 1) - hh  # estimation sample size

    # --- GW predictor forecasts (use raw GW, not standardized) ---
    for (i in 1:ncol(GW)) {
      
      X_p <- cbind(1, GW[1:n_est, i])
      y_p <- r_h[2:(R_oos + p - hh), j]
      res_p <- ols_fit(y_p, X_p)
      FC_PR[p, i, j] <- c(1, GW[R_oos + (p - 1), i]) %*% res_p$beta
      
    }

    # --- SII forecast (linear detrending, recursively estimated) ---
    trend_p <- 1:(R_oos + (p - 1))
    X_lin_p <- cbind(1, trend_p)
    res_lin <- ols_fit(log_EWSI[1:(R_oos + (p - 1))], X_lin_p)
    SII_p <- scale(res_lin$resid)[, 1]

    X_sii_p <- cbind(1, SII_p[1:n_est])
    y_p <- r_h[2:(R_oos + p - hh), j]
    res_p <- ols_fit(y_p, X_sii_p)
    FC_PR[p, ncol(GW) + 1, j] <- c(1, SII_p[R_oos + (p - 1)]) %*% res_p$beta
    
  }
  
}

# ==========================================
# Evaluate forecasts: R2_OS and Clark-West test
# ------------------------------------------

n_total_pred <- ncol(GW) + 1
R2OS <- matrix(NA, n_total_pred, length(h_vec))
CW_tstat <- matrix(NA, n_total_pred, length(h_vec))

for (j in seq_along(h_vec)) {
  
  hh <- h_vec[j]
  n_oos <- P_oos - (hh - 1)

  actual_j <- r_h[(R_oos + 1):(T_full - (hh - 1)), j]
  u_PM_j <- actual_j - FC_PM[1:n_oos]

  for (i in 1:n_total_pred) {
    
    u_PR_ij <- actual_j - FC_PR[1:n_oos, i, j]

    MSFE_PM <- mean(u_PM_j ^ 2)
    MSFE_PR <- mean(u_PR_ij ^ 2)
    R2OS[i, j] <- 100 * (1 - MSFE_PR / MSFE_PM)

    # Clark-West adjusted statistic
    f_CW <- u_PM_j ^ 2 - u_PR_ij ^ 2 + (FC_PM[1:n_oos] - FC_PR[1:n_oos, i, j]) ^ 2
    res_CW <- nwest_ols(f_CW, matrix(1, n_oos, 1), hh)
    CW_tstat[i, j] <- res_CW$tstat[1]
    
  }
  
}

# ==========================================
# Encompassing tests (Harvey-Leybourne-Newbold 1998)
# ------------------------------------------

lambda_hat <- array(NA, dim = c(ncol(GW), 4, length(h_vec)))

for (j in seq_along(h_vec)) {
  
  hh <- h_vec[j]
  n_oos <- P_oos - (hh - 1)
  L_j <- hh

  actual_j <- r_h[(R_oos + 1):(T_full - (hh - 1)), j]
  u_PM_j <- actual_j - FC_PM[1:n_oos]
  u_PR_SII <- actual_j - FC_PR[1:n_oos, n_total_pred, j]

  for (i in 1:ncol(GW)) {
    
    u_PR_i <- actual_j - FC_PR[1:n_oos, i, j]

    # Lambda on SII-based forecast (does SII encompass predictor i?)
    res_enc <- ols_fit(u_PR_i, cbind(u_PR_i - u_PR_SII))
    lam_hat <- res_enc$beta[1]

    d_ij <- (u_PR_i - u_PR_SII) * u_PR_i
    d_bar <- mean(d_ij)
    Q_ij <- sum(d_ij ^ 2) / n_oos
    
    if (L_j > 0) {
      
      for (l in 1:L_j) {
        
        Q_ij <- Q_ij + (1 / n_oos) * (1 - l / (L_j + 1)) * sum(d_ij[(l + 1):n_oos] * d_ij[1:(n_oos - l)])
        
      }
      
    }
    
    HLN <- sqrt(n_oos) * (Q_ij ^ (-0.5)) * d_bar
    lambda_hat[i, 1:2, j] <- c(lam_hat, HLN)

    # Reverse: does predictor i encompass SII?
    res_enc2 <- ols_fit(u_PR_SII, cbind(u_PR_SII - u_PR_i))
    lam_hat2 <- res_enc2$beta[1]

    d_ij2 <- (u_PR_SII - u_PR_i) * u_PR_SII
    d_bar2 <- mean(d_ij2)
    Q_ij2 <- sum(d_ij2 ^ 2) / n_oos
    
    if (L_j > 0) {
      
      for (l in 1:L_j) {
        
        Q_ij2 <- Q_ij2 + (1 / n_oos) * (1 - l / (L_j + 1)) * sum(d_ij2[(l + 1):n_oos] * d_ij2[1:(n_oos - l)])
        
      }
      
    }
    
    HLN2 <- sqrt(n_oos) * (Q_ij2^(-0.5)) * d_bar2
    lambda_hat[i, 3:4, j] <- c(lam_hat2, HLN2)
    
  }
  
}

# ---------------------------------------------------------------------------
# Display results
# ---------------------------------------------------------------------------

pred_names_gw <- c("DP", "DY", "EP", "DE", "RVOL", "BM", "NTIS", "TBL", "LTY", "LTR", "TMS", "DFY", "DFR", "INFL", "SII")

cat("\n=== TABLE 5: Out-of-sample R^2 statistics ===\n")
cat("OOS period: 1990:01-2014:12\n\n")

# R2_OS and Clark-West
cat(sprintf("%-6s", "Pred"))

for (j in seq_along(h_vec)) cat(sprintf(" %10s", paste0("h=", h_vec[j])))

cat("\n")
cat(strrep("-", 6 + 10 * length(h_vec)), "\n")

for (i in 1:n_total_pred) {
  
  cat(sprintf("%-6s", pred_names_gw[i]))
  
  for (j in seq_along(h_vec)) {
    
    stars <- ""
    pval_cw <- 1 - pnorm(CW_tstat[i, j])
    if (pval_cw <= 0.01) stars <- "***"
    else if (pval_cw <= 0.05) stars <- "**"
    else if (pval_cw <= 0.10) stars <- "*"
    cat(sprintf(" %7.2f%s", R2OS[i, j], stars))
    
  }
  
  cat("\n")
  
}

# Encompassing tests
cat("\n--- Encompassing tests (lambda estimates) ---\n")

for (j in seq_along(h_vec)) {
  
  hh <- h_vec[j]
  cat(sprintf("\nh = %d\n", hh))
  cat(sprintf("%-6s %8s %8s | %8s %8s\n", "Pred", "lam_SII", "HLN", "lam_i", "HLN"))
  cat(strrep("-", 48), "\n")
  
  for (i in 1:ncol(GW)) {
    
    stars1 <- ""
    p1 <- 1 - pnorm(lambda_hat[i, 2, j])
    if (p1 <= 0.01) stars1 <- "***"
    else if (p1 <= 0.05) stars1 <- "**"
    else if (p1 <= 0.10) stars1 <- "*"

    stars2 <- ""
    p2 <- 1 - pnorm(lambda_hat[i, 4, j])
    if (p2 <= 0.01) stars2 <- "***"
    else if (p2 <= 0.05) stars2 <- "**"
    else if (p2 <= 0.10) stars2 <- "*"

    cat(sprintf("%-6s %8.2f %6.2f%s | %8.2f %6.2f%s\n",
                pred_names_gw[i],
                lambda_hat[i, 1, j], lambda_hat[i, 2, j], stars1,
                lambda_hat[i, 3, j], lambda_hat[i, 4, j], stars2))
    
  }
  
}
