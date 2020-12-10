#' @export
Exponential <- R6::R6Class("Exponential",
  cloneable = FALSE,
  public = list(
    lambda = NULL,
    initialize = function(lambda) {
      self$lambda <- lambda
    },
    simulate = function(u, X, beta, maxt) {
      beta <- c(log(self$lambda), beta)
      t <- -log(u) / exp(X %*% beta)
      .censor(t = t, maxt = maxt)
    }
  )
)

#' @export
Weibull <- R6::R6Class("Weibull",
  cloneable = FALSE,
  public = list(
    lambda = NULL,
    gamma = NULL,
    initialize = function(lambda, gamma) {
      self$lambda <- lambda
      self$gamma <- gamma
    },
    simulate = function(u, X, beta, maxt) {
      beta <- c(log(self$lambda), beta)
      t <- (-log(u) / exp(X %*% beta))^(1 / self$gamma)
      .censor(t = t, maxt = maxt)
    }
  )
)

#' @export
Gompertz <- R6::R6Class("Gompertz",
  cloneable = FALSE,
  public = list(
    lambda = NULL,
    gamma = NULL,
    initialize = function(lambda, gamma) {
      self$lambda <- lambda
      self$gamma <- gamma
    },
    simulate = function(u, X, beta, maxt) {
      beta <- c(log(self$lambda), beta)
      t <- 1 / self$gamma * log(1 - (self$gamma * log(u)) / exp(X %*% beta))
      .censor(t = t, maxt = maxt)
    }
  )
)
