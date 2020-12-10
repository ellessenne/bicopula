#' @export
simulate_bisurv <- function(dist1, dist2, formula1, formula2, data, beta1, beta2, maxt1 = NULL, maxt2 = NULL, copula = copula::frankCopula(param = 0)) {
  # Create covariates' matrices
  X1 <- stats::model.matrix(formula1, data = data)
  X2 <- stats::model.matrix(formula2, data = data)

  # Check if 'beta*' coefficients match with formulae
  if (!all(names(beta1) %in% colnames(X1))) {
    stop("'beta1' must include coefficients for all variables in 'formula1'", call. = FALSE)
  }
  if (!all(names(beta2) %in% colnames(X2))) {
    stop("'beta2' must include coefficients for all variables in 'formula2'", call. = FALSE)
  }
  beta1 <- stats::na.omit(beta1[colnames(X1)])
  beta2 <- stats::na.omit(beta2[colnames(X2)])

  # Number of rows/observations to simulate
  nsim <- nrow(data)

  # Draw correlated uniforms from a given copula
  UV <- copula::rCopula(n = nsim, copula = copula)

  # Simulate the margins
  Y1 <- dist1$simulate(u = UV[, 1], X = X1, beta = beta1, maxt = maxt1)
  Y1 <- cbind.data.frame(Y1)
  names(Y1) <- paste0(names(Y1), "1")
  # ---
  Y2 <- dist2$simulate(u = UV[, 2], X = X2, beta = beta2, maxt = maxt2)
  Y2 <- cbind.data.frame(Y2)
  names(Y2) <- paste0(names(Y2), "2")

  # Assemble data.frame to return
  out <- cbind.data.frame(Y1, Y2)

  # Return
  return(out)
}
