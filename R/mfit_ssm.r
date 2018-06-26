##' mfit_ssm
##'
##' Fit ssm to multiple individual tracks
##'
##' @param d input data.frame
##' @param ... arguments to be passed to prefilter and fit_ssm
##' @importFrom dplyr group_by do rowwise %>%
##'
##' @export
mfit_ssm <- function(d, ...)
{

  fit <- d %>%
    group_by(id) %>%
    try(do(pf = prefilter(., ...)), silent = TRUE) %>%
    rowwise() %>%
    do(ssm = fit_ssm(.$pf, ...))

  browser()
}
