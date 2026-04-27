###############################################################################
# 00_data.R
# Shared data loading and variable construction for RRZ (2016) replication
# Source this file from each table script: source("00_data.R")
###############################################################################

# Read packages
library(readxl)

# --- File path ---
DATA_FILE <- "data-raw/Returns_short_interest_data.xlsx"

# --- Load GW variables (full history) ---
gw_raw <- read_excel(DATA_FILE, sheet = "GW variables")
# Force all columns to numeric (some are read as character due to scientific notation)
for (col in names(gw_raw)) {
  
  gw_raw[[col]] <- as.numeric(gw_raw[[col]])
  
}

# Sample: 1973:01-2014:12 (rows 1225:1728 in R's 1-indexed data after header)
idx_start <- 1225
idx_end <- 1728
T_full <- idx_end - idx_start + 1  # 504

# --- Equity risk premium ---
# Rfree_lag: lagged risk-free rate (Dec 1972 - Nov 2014)
# R_SP500:   S&P 500 return (Jan 1973 - Dec 2014)
Rfree_lag <- gw_raw$Rfree[(idx_start - 1):(idx_end - 1)]
R_SP500 <- gw_raw$CRSP_SPvw[idx_start:idx_end]
r <- log(1 + R_SP500) - log(1 + Rfree_lag)

# --- Construct 14 Goyal-Welch predictors ---
SP <- gw_raw$Index[idx_start:idx_end]
SP_lag <- gw_raw$Index[(idx_start - 1):(idx_end - 1)]
D12 <- gw_raw$D12[idx_start:idx_end]
E12 <- gw_raw$E12[idx_start:idx_end]

log_DP <- log(D12 / SP)     # 1. Log dividend-price ratio
log_DY <- log(D12 / SP_lag) # 2. Log dividend yield
log_EP <- log(E12 / SP)     # 3. Log earnings-price ratio
log_DE <- log(D12 / E12)    # 4. Log dividend-payout ratio

# 5. Stock return volatility (Mele 2007 estimator)
# MATLAB reads q1215:q1729 (Excel rows) = data rows 1214:1728 = 1972:02-2014:12
# and k1214:k1728 = data rows 1213:1727 = 1972:01-2014:11
SP_R_start <- idx_start - 11   # data row for 1972:02
SP_R <- gw_raw$CRSP_SPvw[SP_R_start:idx_end]
RF_lag_vol <- gw_raw$Rfree[(SP_R_start - 1):(idx_end - 1)]
n_vol <- length(SP_R)
RVOL <- rep(NA, n_vol - 11)

for (t in seq_along(RVOL)) {
  
  RVOL[t] <- mean(abs(SP_R[t:(t + 11)] - RF_lag_vol[t:(t + 11)]))
  
}

RVOL <- sqrt(pi / 2) * sqrt(12) * RVOL  # annualize

BM   <- gw_raw$`b/m`[idx_start:idx_end] # 6. Book-to-market
NTIS <- gw_raw$ntis[idx_start:idx_end]  # 7. Net equity issuance
TBL  <- gw_raw$tbl[idx_start:idx_end]   # 8. Treasury bill rate
LTY  <- gw_raw$lty[idx_start:idx_end]   # 9. Long-term yield
LTR  <- gw_raw$ltr[idx_start:idx_end]   # 10. Long-term return

BAA  <- gw_raw$BAA[idx_start:idx_end]
AAA  <- gw_raw$AAA[idx_start:idx_end]
TMS  <- LTY - TBL # 11. Term spread
DFY  <- BAA - AAA # 12. Default yield spread

CORPR <- gw_raw$corpr[idx_start:idx_end]
DFR   <- CORPR - LTR # 13. Default return spread

INFL_lag <- gw_raw$infl[(idx_start - 1):(idx_end - 1)] # 14. Inflation (lagged)

# Collect into matrix (504 x 14)
GW <- cbind(log_DP, log_DY, log_EP, log_DE, RVOL, BM, NTIS, TBL, LTY, LTR, TMS, DFY, DFR, INFL_lag)
colnames(GW) <- c("DP", "DY", "EP", "DE", "RVOL", "BM", "NTIS", "TBL", "LTY", "LTR", "TMS", "DFY", "DFR", "INFL")

# --- Adjust signs and standardize ---
# Negate NTIS(7), TBL(8), LTY(9), INFL(14) so higher values predict higher returns
GW_adjust <- GW
GW_adjust[, c(7, 8, 9, 14)] <- -GW[, c(7, 8, 9, 14)]
GW_standardize <- scale(GW_adjust)

# --- Principal components (exclude DE=4 and TMS=11) ---
X_pc <- GW_standardize[, -c(4, 11)]
pc_fit <- prcomp(X_pc, center = FALSE, scale. = FALSE)  # already standardized
PC_GW <- scale(pc_fit$x[, 1:3])

# --- Short interest index (SII) ---
si_raw <- read_excel(DATA_FILE, sheet = "Short interest")
EWSI <- si_raw[[2]][1:T_full]
log_EWSI <- log(EWSI)

# Linear detrending: log(EWSI) = a + b*t + u
trend <- 1:T_full
fit_linear <- lm(log_EWSI ~ trend)
SII <- scale(residuals(fit_linear))[, 1]

# --- Cumulative (average) log excess returns ---
h_vec <- c(1, 3, 6, 12)
r_h <- matrix(NA, nrow = T_full, ncol = length(h_vec))
for (j in seq_along(h_vec)) {
  
  hh <- h_vec[j]
  
  for (t in 1:(T_full - (hh - 1))) {
    
    r_h[t, j] <- mean(r[t:(t + hh - 1)])
    
  }
  
}

colnames(r_h) <- paste0("h", h_vec)

# --- Predictor labels ---
pred_names <- c("DP", "DY", "EP", "DE", "RVOL", "BM", "NTIS", "TBL", "LTY", "LTR", "TMS", "DFY", "DFR", "INFL", "SII")

# --- Source helper functions from R/ ---
source("R/nwest-ols.R")
source("R/ols-fit.R")

# --- Key constants ---
IN_SAMPLE_END <- 1989
R_oos <- (IN_SAMPLE_END - 1972) * 12 # 204
P_oos <- T_full - R_oos # 300

cat("Data loaded: T =", T_full, "| R (in-sample) =", R_oos, "| P (out-of-sample) =", P_oos, "\n")
