##' Filter Argos LS or KF track data
##'
##' Internal function
##'
##' @title fit ctrw SSM to individual track data
##'
##' @useDynLib ctrw
##' @importFrom TMB MakeADFun sdreport newtonOption
##' @importFrom stats loess loess.control cov sd predict nlminb
##' @importFrom dplyr mutate filter select full_join arrange lag %>%
##'
##' @export

sfilter <-
  function(d,
           ts = 1,
           v = 2,
           fit.to.subset = TRUE,
           parameters = NULL,
           optim = c("nlminb","optim"),
           verbose = FALSE,
           f = 0.1,
           inner.control = NULL
           ) {
    st <- proc.time()
    call <- match.call()
    optim <- match.arg(optim)
    data.class <- class(d)[2]
 #   cat("\nfitting", data.class, "measurement error model\n")
 #   if(data.class == "LS" & psi != 0) cat("psi will be ignored\n")

    if(!v %in% 1:2) stop("v argument must be 1 or 2 - see ?fit_ssm")

    ## drop any records flagged to be ignored, if fit.to.subset is TRUE
    ## add is.data flag (distinquish obs from reg states)
    if(fit.to.subset) {
      dnew <- d %>%
        filter(.$keep) %>%
        mutate(isd = TRUE)
    } else {
      dnew <- d %>%
        mutate(isd = TRUE)
    }

    ## Interpolation times - assume on ts-multiple of the hour
    tsp <- ts * 3600
    tms <- (as.numeric(d$date) - as.numeric(d$date[1])) / tsp
    index <- floor(tms)
    ts <- data.frame(date = seq(trunc(d$date[1], "hour"), by = tsp, length.out = max(index) + 2))

    ## merge data and interpolation times
    d.all <- full_join(dnew, ts, by = "date") %>%
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
        data = dnew,
        span = f,
        na.action = "na.exclude",
        control = loess.control(surface = "direct")
      )
    fit.y <-
      loess(
        y ~ as.numeric(date),
        data = dnew,
        span = f,
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
          l_psi = 0,
          l_tau = c(0,0),
          l_rho_o = 0
        )
    }

    ## TMB - data list
    fill <- rep(1, nrow(d.all))
    if(data.class == "KF" & v == 1) {
      data <-
        list(
          Y = cbind(d.all$x, d.all$y),
          dt = dt,
          isd = as.integer(d.all$isd),
          obs_mod = 1,
          v = 1,
          K = cbind(fill,fill),
          m = d.all$smin,
          M = d.all$smaj,
          c = d.all$eor
        )
    } else if(data.class == "KF" & v == 2) {
      data <-
        list(
          Y = cbind(d.all$x, d.all$y),
          dt = dt,
          isd = as.integer(d.all$isd),
          obs_mod = 1,
          v = 2,
          K = cbind(fill,fill),
          m = d.all$smin,
          M = d.all$smaj,
          c = d.all$eor
        )
    } else if(data.class == "LS") {
      data <-
        list(
          Y = cbind(d.all$x, d.all$y),
          dt = dt,
          isd = as.integer(d.all$isd),
          obs_mod = 0,
          v = 0,
          m = fill,
          M = fill,
          c = fill,
          K = cbind(d.all$amf_x,d.all$amf_y)
        )
    } else stop("Data class not recognised")

    if(data.class == "KF" & v == 1) {
      map <-
        list(
          l_psi = factor(NA),
          l_tau = factor(c(NA, NA)),
          l_rho_o = factor(NA)
        )
    }
    else if (data.class == "KF" & v == 2) {
      map <-
        list(
          l_tau = factor(c(NA, NA)),
          l_rho_o = factor(NA)
        )
    }
    else if (data.class == "LS") {
      map <- list(l_psi = factor(NA))
    }

## TMB - create objective function
    if(is.null(inner.control)) {
      inner.control <- list(maxit = 1, smartsearch = FALSE)
    }
    obj <-
      MakeADFun(
        data,
        parameters,
        map = map,
        random = "X",
        DLL = "ctrw",
        hessian = TRUE,
        silent = !verbose,
        inner.control = inner.control
      )
#    newtonOption(obj, trace = verbose, smartsearch = TRUE)
#    obj$env$inner.control$trace <- verbose
#    obj$env$inner.control$smartsearch <- FALSE
#    obj$env$inner.control$maxit <- 1
    obj$env$tracemgc <- verbose

## add par values to trace if verbose = TRUE
    myfn <- function(x){print("pars:"); print(x); obj$fn(x)}

    ## Minimize objective function
    opt <-
      suppressWarnings(switch(
        optim,
        nlminb = try(nlminb(obj$par, obj$fn, obj$gr)), #myfn
        optim = try(do.call(optim, args = list(par = obj$par, fn = obj$fn, gr = obj$gr, method = "L-BFGS-B"))) #myfn
      ))

  ## if error then exit with limited output to aid debugging
  if (class(opt) != "try-error") {
    ## Parameters, states and the fitted values
    rep <- sdreport(obj)
    fxd <- summary(rep, "report")

    if (data.class == "KF" & v == 1) {
      fxd <- fxd[1:3, ]
    } else if (data.class == "KF" & v == 2) {
      fxd <- fxd[c(1:3,7),]
    } else if (data.class == "LS") {
      fxd <- fxd[1:6,]
    }

    rdm <-
      matrix(summary(rep, "random"),
             nrow(d.all),
             4,
             dimnames = list(NULL, c("x", "y", "x.se", "y.se")))

    ## Fitted values (estimated locations at observation times)
    fd <- as.data.frame(rdm) %>%
      mutate(id = unique(d.all$id),
             date = d.all$date,
             isd = d.all$isd) %>%
      select(id, date, x, y, x.se, y.se, isd) %>%
      filter(isd) %>%
      select(-isd)

    ## Predicted values (estimated locations at regular time intervals, defined by `ts`)
    pd <- as.data.frame(rdm) %>%
      mutate(id = unique(d.all$id),
             date = d.all$date,
             isd = d.all$isd) %>%
      select(id, date, x, y, x.se, y.se, isd) %>%
      filter(!isd) %>%
      select(-isd)

    if (optim == "nlminb") {
      aic <- 2 * length(opt[["par"]]) + 2 * opt[["objective"]]
    } else {
      aic <- aic <- 2 * length(opt[["par"]]) + 2 * opt[["value"]]
    }
    out <- list(
      call = call,
      predicted = pd,
      fitted = fd,
      par = fxd,
      data = d,
      inits = parameters,
      mmod = data.class,
      ts = tsp / 3600,
      opt = opt,
      tmb = obj,
      rep = rep,
      aic = aic,
      time = proc.time() - st
    )
  } else {
    out <- list(
      call = call,
      data = d,
      inits = parameters,
      mmod = data.class,
      ts = tsp / 3600,
      tmb = obj
    )
  }
    class(out) <- append("ctrwSSM", class(out))

    out
  }
