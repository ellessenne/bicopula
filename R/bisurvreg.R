#' @export
bisurvreg <- function(formula1, formula2, data, distribution1 = "exponential", distribution2 = "exponential", copula = "frank", init = NULL) {
  # Match marginal distributions and copula
  distribution1 <- .match_distribution(x = distribution1)
  distribution2 <- .match_distribution(x = distribution2)
  copula <- .match_copula(x = copula)
  # Process Surv components
  S1 <- .process_surv_formula(formula = formula1, data = data, which = "y")
  start1 <- S1[, which(grepl("^time|^start", colnames(S1))), drop = FALSE]
  status1 <- S1[, which(grepl("^status", colnames(S1))), drop = FALSE]
  S2 <- .process_surv_formula(formula = formula2, data = data, which = "y")
  start2 <- S2[, which(grepl("^time|^start", colnames(S2))), drop = FALSE]
  status2 <- S2[, which(grepl("^status", colnames(S2))), drop = FALSE]
  # Create matrices of covariates
  X1 <- .process_surv_formula(formula = formula1, data = data, which = "x")
  X2 <- .process_surv_formula(formula = formula2, data = data, which = "x")
  #
  list(start1, status1, start2, status2, X1, X2)
  if (is.null(init)) {
    # Obtain starting values by fitting the null model first...
    empty_formula1 <- stats::update(formula1, ~ -.)
    empty_formula2 <- stats::update(formula2, ~ -.)
    empty_X1 <- stats::model.matrix(empty_formula1[-2], data = data)
    empty_X2 <- stats::model.matrix(empty_formula2[-2], data = data)
    # Starting values for ancillary
    empty_ancillary <- .make_ancillary(start1 = start1, start2 = start2, copula = copula, distribution1 = distribution1, distribution2 = distribution2)
    f_init <- TMB::MakeADFun(
      data = list(
        model = "bisurvreg",
        t1 = start1,
        status1 = status1,
        t2 = start2,
        status2 = status2,
        X1 = empty_X1,
        X2 = empty_X2,
        distribution1 = distribution1,
        distribution2 = distribution2,
        copula = copula
      ),
      parameters = list(
        beta1 = rep(0, ncol(empty_X1)),
        beta2 = rep(0, ncol(empty_X2)),
        ancillary = empty_ancillary
      ),
      silent = TRUE,
      DLL = "bicopula_TMBExports"
    )
    empty_fit <- stats::nlminb(start = f_init$par, objective = f_init$fn, gradient = f_init$gr, hessian = f_init$he)
    init <- list(
      beta1 = c(empty_fit$par[grepl("beta1", names(empty_fit$par))], rep(0, ncol(X1) - 1)),
      beta2 = c(empty_fit$par[grepl("beta2", names(empty_fit$par))], rep(0, ncol(X2) - 1)),
      ancillary = empty_fit$par[grepl("ancillary", names(empty_fit$par))]
    )
  }
  # Fit model
  f <- TMB::MakeADFun(
    data = list(
      model = "bisurvreg",
      t1 = start1,
      status1 = status1,
      t2 = start2,
      status2 = status2,
      X1 = X1,
      X2 = X2,
      distribution1 = distribution1,
      distribution2 = distribution2,
      copula = copula
    ),
    parameters = list(
      beta1 = init$beta1,
      beta2 = init$beta2,
      ancillary = init$ancillary
    ),
    silent = TRUE,
    DLL = "bicopula_TMBExports"
  )
  fit <- stats::nlminb(start = f$par, objective = f$fn, gradient = f$gr, hessian = f$he)
  # Repair names
  names(fit$par)[1:ncol(X1)] <- colnames(X1)
  names(fit$par)[(ncol(X1) + 1):(ncol(X1) + ncol(X2))] <- colnames(X2)
  names(fit$par)[ncol(X1) + ncol(X2) + 1] <- "theta"
  return(fit)
}
#       # Fit the empty model here to improve starting values, if not provided with starting values
#       init <- rep(1, ncol(empty_X) + as.numeric(distribution != "exponential"))
#       f <- switch(distribution,
#                   "exponential" = TMB::MakeADFun(data = list(model = "ll_exp", data = empty_X, time = start, status = status), parameters = list(beta = init), silent = TRUE, DLL = "streg_TMBExports"),
#                   "weibull" = TMB::MakeADFun(data = list(model = "ll_wei", data = empty_X, time = start, status = status), parameters = list(beta = init[-length(init)], logp = init[length(init)]), silent = TRUE, DLL = "streg_TMBExports"),
#                   "gompertz" = TMB::MakeADFun(data = list(model = "ll_gom", data = empty_X, time = start, status = status), parameters = list(beta = init[-length(init)], gamma = init[length(init)]), silent = TRUE, DLL = "streg_TMBExports")
#       )
#       # Fit
#       init <- rep(0, ncol(X) + as.numeric(distribution != "exponential"))
#       init[1] <- empty_fit$par[1]
#       if (length(empty_fit$par) > 1) init[length(init)] <- empty_fit$par[2]
#     }
#
#
#
#     names(init) <- switch(distribution,
#                           "exponential" = colnames(X),
#                           "weibull" = c(colnames(X), "ln_p"),
#                           "gompertz" = c(colnames(X), "gamma"),
#     )
#     # Pick correct likelihood function
#     f <- switch(distribution,
#                 "exponential" = TMB::MakeADFun(data = list(model = "ll_exp", data = X, time = start, status = status), parameters = list(beta = init), silent = TRUE, DLL = "streg_TMBExports"),
#                 "weibull" = TMB::MakeADFun(data = list(model = "ll_wei", data = X, time = start, status = status), parameters = list(beta = init[-length(init)], logp = init[length(init)]), silent = TRUE, DLL = "streg_TMBExports"),
#                 "gompertz" = TMB::MakeADFun(data = list(model = "ll_gom", data = X, time = start, status = status), parameters = list(beta = init[-length(init)], gamma = init[length(init)]), silent = TRUE, DLL = "streg_TMBExports")
#     )
#     # Fit
#     model.fit <- stats::nlminb(start = f$par, objective = f$fn, gradient = f$gr, hessian = f$he)
#     # Hessian
#     model.fit$hessian <- f$he(x = model.fit$par)
#     # Fix names of Hessian matrix
#     names(model.fit$par) <- names(init)
#     colnames(model.fit$hessian) <- names(init)
#     rownames(model.fit$hessian) <- names(init)
#
#     # Create object to output
#     out <- list()
#     out$coefficients <- model.fit$par
#     out$vcov <- solve(model.fit$hessian)
#     # Add sum(log(t)) of uncensored observations to log-likelihood
#     out$loglik <- -model.fit$objective
#     out$n <- nrow(X)
#     out$n.events <- sum(status)
#     out$time.at.risk <- sum(start)
#     out$distribution <- distribution
#     out$convergence <- model.fit$convergence
#     out$formula <- formula
#     out$time.range <- range(start)
#     if (x) out$x <- data
#     if (y) out$y <- S
#     out$call <- cl
#     # Return an object of class streg
#     structure(out, class = "streg")
#   }
