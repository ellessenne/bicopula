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
