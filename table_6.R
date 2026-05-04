#!/usr/bin/env Rscript
# ======================================================== #
#
#                  Replication of Table 6 
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

# Replication of Table 6: Out-of-sample CER gains
# RRZ (2016), "Short interest and aggregate stock returns", JFE
#
# Mean-variance investor with relative risk aversion = 3
# Weight constrained to [-0.5, 1.5]
# Volatility forecast: 10-year rolling window of excess returns
# OOS period: 1990:01-2014:12
# Subperiods: pre-crisis (1990:01-2006:12), crisis (2007:01-2014:12)

# Load auxiliary scripts
source("load.R")

# ==========================================
# Cumulative simple excess returns and risk-free rates
# ------------------------------------------

ER <- R_SP500 - Rfree_lag  # simple excess returns

ER_h <- matrix(NA, T_full, length(h_vec))
R_f_h <- matrix(NA, T_full, length(h_vec))
for (j in seq_along(h_vec)) {
  
  hh <- h_vec[j]
  
  for (t in 1:(T_full - (hh - 1))) {
    
    ER_h[t, j]  <- prod(1 + R_SP500[t:(t + hh - 1)]) - prod(1 + Rfree_lag[t:(t + hh - 1)])
    R_f_h[t, j] <- prod(1 + Rfree_lag[t:(t + hh - 1)]) - 1
    
  }
  
}

# ==========================================
#                 Parameters
# ------------------------------------------

RRA <- 3
w_LB <- -0.5
w_UB <- 1.5
vol_window <- 120  # 10 years

# ==========================================
# Out-of-sample forecasts (expanding window, using simple excess returns)
# ------------------------------------------

cat("Computing OOS forecasts for asset allocation...\n")

n_pred_aa <- ncol(GW) + 1  # 14 GW + SII
FC_PM_aa  <- matrix(NA, P_oos, length(h_vec))
FC_PR_aa  <- array(NA, dim = c(P_oos, n_pred_aa, length(h_vec)))
FC_vol    <- matrix(NA, P_oos, length(h_vec))

for (p in 1:P_oos) {
  
  if (p %% 50 == 0) cat("  OOS period", p, "/", P_oos, "\n")

  for (j in seq_along(h_vec)) {
    
    hh <- h_vec[j]
    n_est <- R_oos + (p - 1) - hh

    # Volatility forecast (rolling or expanding)
    # MATLAB: std(ER_h(1:R+p-h(j),j)) or rolling window
    n_vol_end <- R_oos + p - hh  # = R+p-h in MATLAB (1-indexed)
    if (n_vol_end <= vol_window - 1) {
      
      FC_vol[p, j] <- sd(ER_h[1:n_vol_end, j])
      
    } else {
      
      FC_vol[p, j] <- sd(ER_h[(n_vol_end - (vol_window - 1)):n_vol_end, j])
      
    }

    # Prevailing mean forecast
    FC_PM_aa[p, j] <- mean(ER_h[1:n_vol_end, j])

    # GW predictor forecasts
    for (i in 1:ncol(GW)) {
      
      X_p <- cbind(1, GW[1:n_est, i])
      y_p <- ER_h[2:(R_oos + p - hh), j]
      res_p <- ols_fit(y_p, X_p)
      FC_PR_aa[p, i, j] <- c(1, GW[R_oos + (p - 1), i]) %*% res_p$beta
      
    }

    # SII forecast (recursively detrended)
    trend_p <- 1:(R_oos + (p - 1))
    res_lin <- ols_fit(log_EWSI[1:(R_oos + (p - 1))], cbind(1, trend_p))
    SII_p <- scale(res_lin$resid)[, 1]

    X_sii_p <- cbind(1, SII_p[1:n_est])
    y_p <- ER_h[2:(R_oos + p - hh), j]
    res_p <- ols_fit(y_p, X_sii_p)
    FC_PR_aa[p, n_pred_aa, j] <- c(1, SII_p[R_oos + (p - 1)]) %*% res_p$beta
    
  }
  
}

# ==========================================
# Portfolio weights and returns (non-overlapping rebalancing)
# ------------------------------------------

# Read auxiliary functions
source("R/constrain-weight.R")

# Allocate storage
w_PM <- matrix(NA, P_oos, length(h_vec))
R_PM <- matrix(NA, P_oos, length(h_vec))
ER_PM <- matrix(NA, P_oos, length(h_vec))

w_PR <- array(NA, dim = c(P_oos, n_pred_aa, length(h_vec)))
R_PR <- array(NA, dim = c(P_oos, n_pred_aa, length(h_vec)))
ER_PR <- array(NA, dim = c(P_oos, n_pred_aa, length(h_vec)))

R_BH <- matrix(NA, P_oos, length(h_vec))
ER_BH <- matrix(NA, P_oos, length(h_vec))

for (j in seq_along(h_vec)) {
  
  hh <- h_vec[j]
  n_periods <- P_oos %/% hh  # number of non-overlapping periods

  for (tt in 1:n_periods) {
    
    p_idx <- (tt - 1) * hh + 1  # index into OOS arrays
    oos_idx <- R_oos + (tt - 1) * hh + 1  # index into full sample

    vol_t <- FC_vol[p_idx, j]

    # Prevailing mean portfolio
    w_pm_t <- constrain_weight((1 / RRA) * FC_PM_aa[p_idx, j] / vol_t ^ 2)
    w_PM[p_idx, j] <- w_pm_t
    R_PM[p_idx, j] <- R_f_h[oos_idx, j] + w_pm_t * ER_h[oos_idx, j]
    ER_PM[p_idx, j] <- R_PM[p_idx, j] - R_f_h[oos_idx, j]

    # Predictor-based portfolios
    for (i in 1:n_pred_aa) {
      
      w_pr_t <- constrain_weight((1 / RRA) * FC_PR_aa[p_idx, i, j] / vol_t ^ 2)
      w_PR[p_idx, i, j] <- w_pr_t
      R_PR[p_idx, i, j] <- R_f_h[oos_idx, j] + w_pr_t * ER_h[oos_idx, j]
      ER_PR[p_idx, i, j] <- R_PR[p_idx, i, j] - R_f_h[oos_idx, j]
      
    }

    # Buy-and-hold
    R_BH[p_idx, j] <- R_f_h[oos_idx, j] + ER_h[oos_idx, j]
    ER_BH[p_idx, j] <- ER_h[oos_idx, j]
    
  }
  
}

# ==========================================
#        CER gains and Sharpe ratios
# ------------------------------------------

# Read auxiliary functions
source("R/compute-cer-sharpe.R")

# Full OOS period
CER_gain_full <- matrix(NA, n_pred_aa + 1, length(h_vec))  # +1 for buy-and-hold
Sharpe_full <- matrix(NA, n_pred_aa + 2, length(h_vec))  # +2 for PM and B&H

# Pre-crisis: 1990:01-2006:12
GFC_start <- (2006 - 1989) * 12 + 1  # month 205 in OOS
CER_gain_pre <- matrix(NA, n_pred_aa + 1, length(h_vec))
Sharpe_pre <- matrix(NA, n_pred_aa + 2, length(h_vec))

# Crisis: 2007:01-2014:12
CER_gain_gfc <- matrix(NA, n_pred_aa + 1, length(h_vec))
Sharpe_gfc <- matrix(NA, n_pred_aa + 2, length(h_vec))

for (j in seq_along(h_vec)) {
  
  hh <- h_vec[j]

  # Full period
  res_full <- compute_cer_sharpe(R_PR[, , j], ER_PR[, , j], R_PM[, j], ER_PM[, j], hh)
  CER_gain_full[1:n_pred_aa, j] <- res_full$cer_gain
  Sharpe_full[1, j] <- res_full$sharpe_PM
  Sharpe_full[2:(n_pred_aa + 1), j] <- res_full$sharpe

  # Buy-and-hold
  R_bh <- R_BH[is.finite(R_BH[, j]), j]
  ER_bh <- ER_BH[is.finite(ER_BH[, j]), j]
  R_pm <- R_PM[is.finite(R_PM[, j]), j]
  CER_BH <- (12 / hh) * (mean(R_bh) - 0.5 * RRA * sd(R_bh) ^ 2)
  CER_PM_val <- (12 / hh) * (mean(R_pm) - 0.5 * RRA * sd(R_pm) ^ 2)
  CER_gain_full[n_pred_aa + 1, j] <- 100 * (CER_BH - CER_PM_val)
  Sharpe_full[n_pred_aa + 2, j] <- sqrt(12 / hh) * mean(ER_bh) / sd(ER_bh)

  # Pre-crisis period
  pre_idx <- 1:(GFC_start - 1)
  res_pre <- compute_cer_sharpe(R_PR[pre_idx, , j], ER_PR[pre_idx, , j], R_PM[pre_idx, j], ER_PM[pre_idx, j], hh)
  CER_gain_pre[1:n_pred_aa, j] <- res_pre$cer_gain
  Sharpe_pre[1, j] <- res_pre$sharpe_PM
  Sharpe_pre[2:(n_pred_aa + 1), j] <- res_pre$sharpe

  R_bh_pre <- R_BH[pre_idx, j]; R_bh_pre <- R_bh_pre[is.finite(R_bh_pre)]
  ER_bh_pre <- ER_BH[pre_idx, j]; ER_bh_pre <- ER_bh_pre[is.finite(ER_bh_pre)]
  R_pm_pre <- R_PM[pre_idx, j]; R_pm_pre <- R_pm_pre[is.finite(R_pm_pre)]
  CER_BH_pre <- (12 / hh) * (mean(R_bh_pre) - 0.5 * RRA * sd(R_bh_pre) ^ 2)
  CER_PM_pre <- (12 / hh) * (mean(R_pm_pre) - 0.5 * RRA * sd(R_pm_pre) ^ 2)
  CER_gain_pre[n_pred_aa + 1, j] <- 100 * (CER_BH_pre - CER_PM_pre)
  Sharpe_pre[n_pred_aa + 2, j] <- sqrt(12 / hh) * mean(ER_bh_pre) / sd(ER_bh_pre)

  # Crisis period
  gfc_idx <- GFC_start:P_oos
  res_gfc <- compute_cer_sharpe(R_PR[gfc_idx, , j], ER_PR[gfc_idx, , j], R_PM[gfc_idx, j], ER_PM[gfc_idx, j], hh)
  CER_gain_gfc[1:n_pred_aa, j] <- res_gfc$cer_gain
  Sharpe_gfc[1, j] <- res_gfc$sharpe_PM
  Sharpe_gfc[2:(n_pred_aa + 1), j] <- res_gfc$sharpe

  R_bh_gfc <- R_BH[gfc_idx, j]
  R_bh_gfc <- R_bh_gfc[is.finite(R_bh_gfc)]
  ER_bh_gfc <- ER_BH[gfc_idx, j]
  ER_bh_gfc <- ER_bh_gfc[is.finite(ER_bh_gfc)]
  R_pm_gfc <- R_PM[gfc_idx, j]
  R_pm_gfc <- R_pm_gfc[is.finite(R_pm_gfc)]
  CER_BH_gfc <- (12 / hh) * (mean(R_bh_gfc) - 0.5 * RRA * sd(R_bh_gfc) ^ 2)
  CER_PM_gfc <- (12 / hh) * (mean(R_pm_gfc) - 0.5 * RRA * sd(R_pm_gfc) ^ 2)
  CER_gain_gfc[n_pred_aa + 1, j] <- 100 * (CER_BH_gfc - CER_PM_gfc)
  Sharpe_gfc[n_pred_aa + 2, j] <- sqrt(12 / hh) * mean(ER_bh_gfc) / sd(ER_bh_gfc)
  
}

# ==========================================
#         Display Table 6: CER gains
# ------------------------------------------

pred_names_t6 <- c("DP", "DY", "EP", "DE", "RVOL", "BM", "NTIS", "TBL", "LTY", "LTR", "TMS", "DFY", "DFR", "INFL", "SII", "Buy and hold")

cat("\n=== TABLE 6: Out-of-sample CER gains (basis points, annualized) ===\n\n")

print_cer_table <- function(mat, title, pred_nms) {
  
  cat(title, "\n")
  cat(sprintf("%-14s", "Predictor"))
  
  for (j in seq_along(h_vec)) cat(sprintf(" %8s", paste0("h=", h_vec[j])))
  
  cat("\n")
  cat(strrep("-", 14 + 8 * length(h_vec)), "\n")
  
  for (i in 1:nrow(mat)) {
    
    cat(sprintf("%-14s", pred_nms[i]))
    
    for (j in seq_along(h_vec)) cat(sprintf(" %8.2f", mat[i, j]))
    
    cat("\n")
    
  }
  
  cat("\n")
  
}

print_cer_table(CER_gain_full, "Full OOS: 1990:01-2014:12", pred_names_t6)
print_cer_table(CER_gain_pre,  "Pre-crisis: 1990:01-2006:12", pred_names_t6)
print_cer_table(CER_gain_gfc,  "Crisis: 2007:01-2014:12", pred_names_t6)
