nwest_ols <- function(y, X, nlag) {

  # OLS regression with Newey-West HAC standard errors (Bartlett kernel)
  # y: T x 1 vector, X: T x k matrix (should include intercept column)
  n <- length(y)
  k <- ncol(X)
  XX_inv <- solve(crossprod(X))
  beta <- XX_inv %*% crossprod(X, y)
  resid <- y - X %*% beta
  SSE <- sum(resid^2)
  SST <- sum((y - mean(y))^2)
  rsqr <- 1 - SSE / SST

  # Newey-West HAC covariance
  S <- matrix(0, k, k)
  e_X <- X * as.vector(resid)
  S <- crossprod(e_X) / n

  if (nlag > 0) {

    for (l in 1:nlag) {

      w <- 1 - l / (nlag + 1)
      Gamma_l <- crossprod(e_X[(l + 1):n, , drop = FALSE], e_X[1:(n - l), , drop = FALSE]) / n
      S <- S + w * (Gamma_l + t(Gamma_l))

    }

  }

  V <- n * XX_inv %*% S %*% XX_inv
  se <- sqrt(diag(V))
  tstat <- beta / se

  list(beta = as.vector(beta), tstat = as.vector(tstat),
       se = as.vector(se), resid = as.vector(resid),
       rsqr = rsqr, SSE = SSE)

}
