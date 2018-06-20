##' Continuous-time Random Walk Filter
##'
##' Fit a simple random walk in continuous time to filter Argos KF of LS data and predict
##' locations on a regular time step.
##'
##'
##' @title Correlated Random Walk Filter
##' @param d a data frame of observations including Argos KF error ellipse info
##' @param ts the time step, in hours, to predict to
##' @param fit.to.subset a logical vector (default is TRUE) indicating whether the SSM is to be fit to the data subset determined by prefilter
##' @param optim numerical optimizer
##' @param verbose report progress during minimization
##' @param span the span parameter for the loess fits used to estimate
##'   initial locations
##' @return a list with components
##' \item{\code{predicted}}{a data.frame of predicted location states}
##' \item{\code{fitted}}{a data.frame of fitted locations}
##' \item{\code{par}}{model parameter summmary}
##' \item{\code{data}}{the input data.frame}
##' \item{\code{subset}}{the inpu subset vector}
##' \item{\code{ts}}{the prediction time step}
##' \item{\code{opt}}{the object returned by the optimizer}
##' \item{\code{tmb}}{the TMB object}
##' \item{\code{aic}}{the calculated Akaike Information Criterion}
##'
##' @examples
##' \dontrun{
##' require(dplyr)
##' data(ellie)
##' fit <- ellie %>%
##'     prefilter(., min.dist = 100) %>%
##'     fit_ssm(ts = 6)
##'
##' ## fit LS measurement model
##' fit.ls <- ellie %>%
##'     prefilter(., min.dist = 100) %>%
##'     fit_ssm(ts = 6)
##' }
##'
##' @useDynLib ctrw
##' @importFrom TMB MakeADFun sdreport newtonOption
##' @importFrom stats loess loess.control cov sd predict nlminb
##' @importFrom dplyr mutate filter select full_join arrange lag %>%
##' @importFrom geosphere mercator
##'
##' @export

fit_ssm <-
  function(d,
           ts = 1,
           fit.to.subset = TRUE,
           parameters = NULL,
           optim = c("nlminb", "optim"),
           verbose = FALSE,
           span = 0.1) {

    optim <- match.arg(optim)

    ## drop any records flagged to be ignored
    ## add is.data flag (distinquish obs from reg states)
    d <- d %>%
      if(fit.to.subset) filter(.$keep) %>%
      mutate(isd = TRUE)

    ## Interpolation times - assume on ts-multiple of the hour
    tsp <- ts * 3600
    tms <- (as.numeric(d$date) - as.numeric(d$date[1])) / tsp
    index <- floor(tms)
    ts <- data.frame(date = seq(trunc(d$date[1], "hour"), by = tsp, length.out = max(index) + 2))

    ## merge data and interpolation times
    d.all <- full_join(d, ts, by = "date") %>%
      arrange(date) %>%
      mutate(isd = ifelse(is.na(isd), FALSE, isd)) %>%
      mutate(id = ifelse(is.na(id), na.omit(unique(id))[1], id))

    class(d.all$x) <- NA_real_
    class(d.all$y) <- NA_real_

    ## calc delta times in hours for observations & interpolation points (states)
    dt <- difftime(d.all$date, lag(d.all$date), units = "hours") %>%
      as.numeric() / 24
    dt[1] <- 0

    ## Predict track from loess smooths (state initial values)
    ## TP -data here should only be at obs times.
    fit.x <-
      loess(
        x ~ as.numeric(date),
        data = d,
        span = span,
        na.action = "na.exclude",
        control = loess.control(surface = "direct")
      )
    fit.y <-
      loess(
        y ~ as.numeric(date),
        data = d,
        span = span,
        na.action = "na.exclude",
        control = loess.control(surface = "direct")
      )

    ## Predict track, increments and stochastic innovations
    ## for interp this predict call needs a newdata = all dates ( interp and obs combined times )
    xs <-
      cbind(predict(fit.x, newdata = data.frame(date = as.numeric(d.all$date))),
            predict(fit.y, newdata = data.frame(date = as.numeric(d.all$date)))
            )

    if (is.null(parameters)) {
      ## Estimate stochastic innovations
      es <- xs[-1,] - xs[-nrow(xs),]

      ## Estimate components of variance
      V <- cov(es)
      sigma <- sqrt(diag(V))
      rho <- V[1, 2] / prod(sqrt(diag(V)))

      parameters <-
        list(
          l_sigma = log(pmax(1e-08, sigma)),
          l_rho_p = log((1 + rho) / (1 - rho)),
          X = xs,
          l_tau = c(0,0),
          l_rho_o = 0
        )
    }

    ## TMB - data list
    data <-
      list(
        Y = cbind(d.all$x, d.all$y),
        dt = dt,
        isd = as.integer(d.all$isd),
        obs_mod = ifelse(class(d)[2] == "KF", 1, 0),
        m = d.all$smin,
        M = d.all$smaj,
        c = d.all$eor,
        K = cbind(d.all$amf_x, d.all$amf_y)
      )

    if(class(d)[2] == "KF") map = list(l_tau = factor(c(NA,NA)), l_rho_o = factor(NA))
    else map = list()


    ## TMB - create objective function
    obj <-
      MakeADFun(
        data,
        parameters,
        map = map,
        random = "X",
        DLL = "ctrw",
        hessian = TRUE,
        silent = !verbose
      )
    obj$env$inner.control$trace <- verbose
    obj$env$tracemgc <- verbose

    ## Minimize objective function
    opt <-
      suppressWarnings(switch(
        optim,
        nlminb = nlminb(obj$par, obj$fn, obj$gr),
        optim = do.call("optim", obj)
      ))

    ## Parameters, states and the fitted values
    rep <- sdreport(obj)
    fxd <- summary(rep, "report")

    rdm <-
      matrix(summary(rep, "random"),
             nrow(d.all),
             4,
             dimnames = list(NULL, c("x", "y", "x.se", "y.se")))

    ## Fitted values (estimated locations at observation times)
    fd <- as.data.frame(rdm) %>%
      mutate(id = unique(d.all$id), date = d.all$date, isd = d.all$isd) %>%
      select(id, date, x, y, x.se, y.se, isd) %>%
      filter(isd) %>%
      select(-isd)

    ## Predicted values (estimated locations at regular time intervals, defined by `ts`)
    pd <- as.data.frame(rdm) %>%
      mutate(id = unique(d.all$id), date = d.all$date, isd = d.all$isd) %>%
      select(id, date, x, y, x.se, y.se, isd) %>%
      filter(!isd) %>%
      select(-isd)

    if (optim == "nlminb")
      aic <- 2 * length(opt[["par"]]) + 2 * opt[["objective"]]

    out <- list(
      predicted = pd,
      fitted = fd,
      par = fxd,
      data = d,
      subset = subset,
      opt = opt,
      tmb = obj,
      aic = aic
    )
    class(out) <- append("ctrw", class(out))
  }
