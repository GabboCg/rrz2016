ols_fit <- function(y, X) {

  XX_inv <- solve(crossprod(X))
  beta <- XX_inv %*% crossprod(X, y)
  resid <- y - X %*% beta

  list(beta = as.vector(beta), resid = as.vector(resid))

}
