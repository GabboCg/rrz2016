compute_cer_sharpe <- function(R_mat, ER_mat, R_PM_vec, ER_PM_vec, h_val, RRA = 3) {

  # Returns: CER gain (in %, annualized) and Sharpe ratio
  R_pm <- R_PM_vec[is.finite(R_PM_vec)]
  ER_pm <- ER_PM_vec[is.finite(ER_PM_vec)]
  CER_PM <- (12 / h_val) * (mean(R_pm) - 0.5 * RRA * sd(R_pm) ^ 2)
  Sharpe_PM <- sqrt(12 / h_val) * mean(ER_pm) / sd(ER_pm)

  n_preds <- if (is.matrix(R_mat)) ncol(R_mat) else 1
  cer_gain <- rep(NA, n_preds)
  sharpe <- rep(NA, n_preds)

  for (i in 1:n_preds) {

    R_i <- if (is.matrix(R_mat)) R_mat[, i] else R_mat
    ER_i <- if (is.matrix(ER_mat)) ER_mat[, i] else ER_mat
    R_i <- R_i[is.finite(R_i)]
    ER_i <- ER_i[is.finite(ER_i)]
    CER_i <- (12 / h_val) * (mean(R_i) - 0.5 * RRA * sd(R_i)^2)
    cer_gain[i] <- 100 * (CER_i - CER_PM)
    sharpe[i] <- sqrt(12 / h_val) * mean(ER_i) / sd(ER_i)

  }

  list(cer_gain = cer_gain, sharpe = sharpe, sharpe_PM = Sharpe_PM)

}
