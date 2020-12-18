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
