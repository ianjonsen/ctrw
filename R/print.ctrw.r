##' Print \code{ctrw}
##'
##' @method print ctrw
##'
##' @param x a \code{ctrw} fit object
##' @param digits number of digits to use in display
##' @param ... unused. For compatibility with the generic method.
##'
##' @export

print.ctrw <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  mmod <- class(x$data)[2]
  nbrStates <- nrow(x$predicted)
  nbStates <- nrow(x$fitted)
  parm <- x$par

  print(x$call)

  cat("\n")
  cat("negative log-likelihood:", logLik(x)$val)

}
