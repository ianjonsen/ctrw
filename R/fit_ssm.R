##' @title Fit a Continuous-time Random Walk state-space filter to Argos data
##'
##' @description fits a simple random walk in continuous time to filter Argos KF or LS data and predict
##' locations on a regular time step
##'
##' @param d a data frame of observations including Argos KF error ellipse info
##' @param span degree of loess smoothing (range: 0 - 1) to identify potential outliers in prefilter
##' @param min.dist minimum distance from track to define potential outlier locations in prefilter
##' @param ... arguments passed to sfilter, described below:
##' @param ts the time step, in hours, to predict to
##' @param fit.to.subset a logical vector (default is TRUE) indicating whether the SSM is to be fit to the data subset determined by prefilter
##' @param psi a logical scalar (default is FALSE) indicating whether the KF measurement error model is to be fit with a scale parameter for the error ellipses
##' @param optim numerical optimizer
##' @param verbose report progress during minimization
##' @param f the span parameter for the loess fits used to estimate initial location states
##'
##' @return a list with components
##' \item{\code{call}}{the matched call}
##' \item{\code{predicted}}{a data.frame of predicted location states}
##' \item{\code{fitted}}{a data.frame of fitted locations}
##' \item{\code{par}}{model parameter summmary}
##' \item{\code{data}}{the input data.frame}
##' \item{\code{subset}}{the inpu subset vector}
##' \item{\code{mem}}{the measurement error model used}
##' \item{\code{ts}}{the prediction time interval}
##' \item{\code{opt}}{the object returned by the optimizer}
##' \item{\code{tmb}}{the TMB object}
##' \item{\code{rep}}{TMB sdreport}
##' \item{\code{aic}}{the calculated Akaike Information Criterion}
##'
##' @examples
##' \dontrun{
##' require(dplyr)
##' data(ellie)
##' ## fit KF measurement model
##' fit <- ellie %>%
##'     select(1:8) %>%
##'     prefilter(., min.dist = 100) %>%
##'     fit_ssm(ts = 6)
##'
##' ## fit LS measurement model
##' fit.ls <- ellie %>%
##'     select(1:5) %>%
##'     prefilter(., min.dist = 100) %>%
##'     fit_ssm(ts = 6)
##' }
##' @importFrom dplyr group_by do rowwise %>%
##'
##' @export
fit_ssm <- function(d,
                    span = 0.01,
                    min.dist = 100,
                    ...
                    )
{

  d %>%
    group_by(id) %>%
    do(pf = prefilter(., span = span, min.dist = min.dist)) %>%
    rowwise() %>%
    do(ssm = try(sfilter(.$pf, ...), silent = TRUE))

}
