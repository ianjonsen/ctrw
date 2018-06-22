##' @importFrom stats pnorm
##' @importFrom lme4 nobars findbars
##' @importFrom lme4 nobars
##' @importFrom dplyr %>%
##' @method summary ctrwSSM
##' @export
summary.ctrwSSM <- function(x, digits = 3, ...) {
  if (length(list(...)) > 0) {
    warning("additional arguments ignored")
  }

  mmod <- class(x$data)[2]
  nbrStates <- nrow(x$predicted)
  nbStates <- nrow(x$fitted)
  parm <- x$par
  resid <- list(x = x$fitted$x - subset(x$data, keep)$x,
                y = x$fitted$y - subset(x$data, keep)$y)

  cat("negative log-likelihood:", x$opt$objective, "\n")
  cat("convergence:", x$opt$message, "\n\n")
  cat("Argos measurement error model:", mmod, "\n")
  cat("number of observations:", nbStates, "\n")
  cat("number of regularised state estimates:", nbrStates, "\n\n")
  cat("parameter estimates\n")
  cat("-------------------\n")
  print(parm, digits = digits, justify = "right");cat("\n")

  cat("quantiles of residuals: x\n")
  cat("----------------------------\n")
  print(quantile(resid$x), digits = digits); cat("\n")

  cat("quantiles of residuals: y\n")
  cat("----------------------------\n")
  print(quantile(resid$y), digits = digits); cat("\n")


}
