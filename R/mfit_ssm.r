##' mfit_ssm
##'
##' Fit ssm to multiple individual tracks
##'
##' @param d input data.frame
##' @param span degree of loess smoothing (range: 0 - 1) to identify potential outliers in prefilter
##' @param min.dist minimum distance from track to define potential outlier locations
##' @param ... arguments to be passed to fit_ssm
##' @importFrom dplyr group_by do rowwise %>%
##'
##' @export
mfit_ssm <- function(d, span = 0.01, min.dist = 100, ...)
{

  fit <- d %>%
    group_by(id) %>%
    try(do(pf = prefilter(., span = span, min.dist = min.dist)), silent = TRUE) %>%
    rowwise() %>%
    do(ssm = fit_ssm(.$pf, ...))

  browser()
}
