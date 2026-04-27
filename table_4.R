#!/usr/bin/env Rscript
# ======================================================== #
#
#                  Replication of Table 4 
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

# Replication of Table 4: Alternative detrending methods
# RRZ (2016), "Short interest and aggregate stock returns", JFE
#
# Tests SII predictive power with linear, quadratic, cubic, and stochastic
# detrending of log(EWSI). Reports beta, t-stat (NW-HAC), R^2 for h=1,3,6,12

source("load.R")

# ==========================================
# Construct SII under different detrending methods (full sample)
# ------------------------------------------

trend_t <- 1:T_full
trend_t2 <- trend_t ^ 2
trend_t3 <- trend_t ^ 3

# Linear (already in SII from 00_data.R)
SII_linear <- SII

# Quadratic
fit_quad <- lm(log_EWSI ~ trend_t + trend_t2)
SII_quad <- scale(residuals(fit_quad))[, 1]

# Cubic
fit_cub <- lm(log_EWSI ~ trend_t + trend_t2 + trend_t3)
SII_cub <- scale(residuals(fit_cub))[, 1]

# Stochastic (60-month backward-looking moving average)
MA_size <- 60
SII_stoch <- rep(NA, T_full)
for (t in MA_size:T_full) {
  
  SII_stoch[t] <- log_EWSI[t] - mean(log_EWSI[(t - MA_size + 1):t])
  
}

SII_stoch[MA_size:T_full] <- scale(SII_stoch[MA_size:T_full])[, 1]

# ==========================================
# In-sample regressions for each detrending method
# ------------------------------------------

methods <- list(
  "Linear"     = list(sii = SII_linear, start = 1),
  "Quadratic"  = list(sii = SII_quad,   start = 1),
  "Cubic"      = list(sii = SII_cub,    start = 1),
  "Stochastic" = list(sii = SII_stoch,  start = MA_size)
)

cat("\n=== TABLE 4: Alternative detrending methods ===\n")
cat("Sample: 1973:01-2014:12 (stochastic starts 1977:12)\n\n")

results_table <- matrix(NA, nrow = length(methods), ncol = 3 * length(h_vec))
rownames(results_table) <- names(methods)

for (m in seq_along(methods)) {
  
  mname <- names(methods)[m]
  sii_m <- methods[[m]]$sii
  s0 <- methods[[m]]$start

  for (j in seq_along(h_vec)) {
    
    hh <- h_vec[j]

    # Align: y = r_h at t+1:t+h, predictor at t
    # For stochastic, first valid obs is MA_size
    t_start <- s0
    t_end   <- T_full - hh

    y_j  <- 100 * r_h[(t_start + 1):(t_end + 1), j]
    x_j  <- -sii_m[t_start:t_end]  # negate SII
    X_ij <- cbind(1, x_j)

    res <- nwest_ols(y_j, X_ij, hh)

    col_base <- (j - 1) * 3
    results_table[m, col_base + 1] <- res$beta[2]
    results_table[m, col_base + 2] <- res$tstat[2]
    results_table[m, col_base + 3] <- 100 * res$rsqr
    
  }
  
}

# Display
for (j in seq_along(h_vec)) {
  
  hh <- h_vec[j]
  cat(sprintf("--- h = %d ---\n", hh))
  cat(sprintf("%-12s %8s %8s %8s\n", "Method", "beta(-)", "t-stat", "R2(%)"))
  cat(strrep("-", 38), "\n")
  
  for (m in seq_along(methods)) {
    
    col_base <- (j - 1) * 3
    cat(sprintf("%-12s %8.2f  [%5.2f] %8.2f\n",
                names(methods)[m],
                results_table[m, col_base + 1],
                results_table[m, col_base + 2],
                results_table[m, col_base + 3]))
  }
  
  cat("\n")
  
}

# Also display in the paper's format (all horizons side by side)
cat("\n--- Combined Table 4 format ---\n")
cat(sprintf("%-12s", "Method"))

for (j in seq_along(h_vec)) {
  
  cat(sprintf("  %6s %8s %7s", "beta", "[t]", "R2(%)"))
  
}

cat("\n")
cat(strrep("-", 12 + 4 * 23), "\n")

for (m in seq_along(methods)) {
  
  cat(sprintf("%-12s", names(methods)[m]))
  
  for (j in seq_along(h_vec)) {
    
    col_base <- (j - 1) * 3
    cat(sprintf("  %6.2f  [%5.2f] %7.2f",
                results_table[m, col_base + 1],
                results_table[m, col_base + 2],
                results_table[m, col_base + 3]))
    
  }
  
  cat("\n")
  
}
