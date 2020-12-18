#' @keywords internal
.censor <- function(t, maxt) {
  if (is.null(maxt)) {
    status <- rep(1, length(t))
  } else {
    status <- as.numeric(t < maxt)
    t <- pmin(t, maxt)
  }
  list(eventtime = t, status = status)
}

#' @keywords internal
.process_surv_formula <- function(formula, data, which) {
  if (which == "y") {
    out <- eval(expr = formula[[2]], envir = data)
  } else if (which == "x") {
    out <- stats::model.matrix(formula[-2], data = data)
  }
  return(out)
}

#' @keywords internal
.match_distribution <- function(x) {
  x <- match.arg(x, choices = c("exponential", "weibull", "gompertz"))
  out <- switch(x,
    "exponential" = 1,
    "weibull" = 2,
    "gompertz" = 3
  )
  return(out)
}

#' @keywords internal
.match_copula <- function(x) {
  x <- match.arg(x, choices = c("frank"))
  out <- switch(x,
    "frank" = 1
  )
  return(out)
}

#' @keywords internal
.make_ancillary <- function(start1, start2, copula, distribution1, distribution2) {
  theta <- switch(copula,
    copula::iTau(copula = copula::frankCopula(dim = 2), tau = stats::cor(start1, start2, method = "kendall"))
  )
  if (distribution1 == 1) {
    anc1 <- NA
  } else {
    anc1 <- 1
  }
  if (distribution2 == 1) {
    anc2 <- NA
  } else {
    anc2 <- 1
  }
  out <- c(theta, anc1, anc2)
  out <- out[is.finite(out)]
  return(out)
}
