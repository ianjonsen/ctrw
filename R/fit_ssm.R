##' @title Fit a Continuous-time Random Walk state-space filter to Argos data
##'
##' @description fits a simple random walk in continuous time to filter Argos KF or LS data and predict
##' locations on a regular time step
##'
##' @param d a data frame of observations including Argos KF error ellipse info
##' @param span degree of loess smoothing (range: 0 - 1) to identify potential outliers in prefilter
##' @param min.dt minimum allowable time difference between observations; dt <= min.dt will be ignored by the SSM
##' @param min.dist minimum distance from track to define potential outlier locations in prefilter
##' ##' @param ptime the regular time interval, in hours, to predict to. Alternatively, a vector of prediction times, possibly not regular, must be specified as a data.frame with id and POSIXt dates.
##' @param ... arguments passed to sfilter, described below:
##' @param fit.to.subset fit the SSM to the data subset determined by prefilter (default is TRUE)
##' @param psi estimate scaling parameter for the KF measurement error model error ellipses (0 = no psi, default; 1 = single psi for semi-minor axis)
##' @param pf just pre-filter the data, do not fit the ctrw (default is FALSE)
##' @param optim numerical optimizer (nlminb or optim)
##' @param verbose report progress during minimization
##'
##' @return a list with components
##' \item{\code{call}}{the matched call}
##' \item{\code{predicted}}{a data.frame of predicted location states}
##' \item{\code{fitted}}{a data.frame of fitted locations}
##' \item{\code{par}}{model parameter summmary}
##' \item{\code{data}}{the input data.frame}
##' \item{\code{subset}}{the inpu subset vector}
##' \item{\code{mem}}{the measurement error model used}
##' \item{\code{opt}}{the object returned by the optimizer}
##' \item{\code{tmb}}{the TMB object}
##' \item{\code{rep}}{TMB sdreport}
##' \item{\code{aic}}{the calculated Akaike Information Criterion}
##' \item{\code{time}}{the processing time for sfilter}
##'
##' @examples
##' \dontrun{
##' require(dplyr)
##' data(ellie)
##' ## fit KF measurement error model
##' fkf <- fit_ssm(ellie, min.dist = 150, ptime = 12, psi = 0)
##'
##' ## fit KFp measurement error model
##' fkfp <- fit_ssm(ellie, min.dist = 150, ptime = 12, psi = 1)
##'
##' ## fit LS measurement error model
##' fls <- fit_ssm(ellie[, 1:5], min.dist = 150, ptime = 12)
##' }
##' @importFrom dplyr group_by do rowwise %>% ungroup select mutate slice
##' @importFrom tibble as_tibble
##'
##' @export
fit_ssm <- function(d,
                    span = 0.1,
                    min.dt = 60,
                    min.dist = 100,
                    pf = FALSE,
                    ptime,
                    ...
                    )
{
  if(is.null(ptime)) print("\nNo ptime specified, using 6 h as a default time step")
  else if(length(ptime) > 1 & !is.data.frame(ptime))
    stop("\nptime must be a data.frame with id's when specifying multiple prediction times")
  else if(length(ptime) > 1 & is.data.frame(ptime)) {
    if(sum(!names(ptime) %in% c("id","date")) > 0) stop("\n ptime names must be `id` and `date`")
  }

  fit <- d %>%
    group_by(id) %>%
    do(pf = prefilter(., span = span, min.dt = min.dt, min.dist = min.dist))

  if(pf){
    fit <- do.call(rbind, fit$pf) %>%
      as_tibble()
  } else {
    fit <- fit %>%
      rowwise() %>%
      do(ssm = try(sfilter(.$pf, ptime = ptime, ...), silent = TRUE))

    fail <- which(sapply(fit$ssm, length) == 6)
    if (length(fail) > 0) {
      cat(sprintf("\n%d optimisation failures\n", length(fail)))
    }

    fit <- fit %>%
      ungroup() %>%
      mutate(id = sapply(.$ssm, function(x)
        x$data$id[1])) %>%
      select(id, ssm)
  }

  fit
}
