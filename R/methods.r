##' @importFrom stats logLik
##' @export
logLik.ctrw <- function(m, ...) {
  if (!is.null(m$rep)) {
    val <- if (m$rep$pdHess) {
      -1 * m$opt$objective
    } else {
      NA
    }
  } else {
    val <- -1 * m$opt$objective
  }

  nobs <- nrow(m$fitted)
  structure(
    val = val,
    nobs = nobs,
    nall = nobs,
    df = length(m$tmb$par),
    class = "logLik"
  )
}
