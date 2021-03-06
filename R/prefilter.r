##' @title Prepare Argos data for fitting a ctrw model
##'
##' @description \code{prefilter()} (1) determines Argos data type (LS or KF); (2) converts dates to POSIXt;
##' identifies observations with duplicate dates; (3) orders observations in time; (4) converts
##' longitudes from 0,360 to -180,180; (5) projects lonlat coords to mercator x,y coords (in km);
##' (6) adds location error multiplication factors based on Argos location class (for type LS);
##' and (7) uses a loess smooth to identify potential outlier locations (by distance only) to be ignored when fitting
##' the \code{ctrw} model
##'
##' @details Internal function
##'
##' @param d input data - must have 5 (LS), or 8 (KF) columns (see details)
##' @param span degree of loess smoothing (range: 0 - 1) to identify potential outliers
##' @param min.dt minimum allowable time difference between observations; dt < min.dt will be ignored by the SSM
##' @param min.dist minimum distance from track to define potential outlier locations
##' @param time.gap not currently implemented
##' @importFrom lubridate ymd_hms
##' @importFrom stats loess
##' @importFrom dplyr mutate distinct arrange filter select %>% left_join
##' @importFrom rgdal project
##'
##' @export

prefilter <- function(d, span = 0.01, min.dt = 0, min.dist = 100, time.gap = NULL) {

  # check input data
  if(!ncol(d) %in% c(5,8)) stop("Data can only have 5 (for LS data) or 8 (for KF data) columns")

  if((ncol(d) == 5 &
      !isTRUE(all.equal(names(d), c("id", "date", "lc", "lon", "lat")))) ||
     (ncol(d) == 8 &
      !isTRUE(all.equal(
        names(d),
        c("id", "date", "lc", "lon", "lat", "smaj", "smin", "eor")
      )))) stop("Unexpected column names in Data, type `?prefilter` for details on data format")

  if(length(unique(d$id)) > 1) stop("Multiple individual tracks in Data, use fit_ssm")

  if(!is.null(d$id)) d <- d %>% mutate(id = as.character(id))

  if(ncol(d) == 5 | (ncol(d) == 8 & sum(is.na(d$smaj))/nrow(d) == 1)) data.type <- "LS"
  else data.type <- "KF"


  ##  convert dates to POSIXt
  ##  flag any duplicate date records,
  ##  order records by time,
  ##  set lc to ordered factor
  ##  convert lon from 0,360 to -180,180
  d <- d %>%
    mutate(date = ymd_hms(date, tz = "GMT")) %>%
    mutate(keep = difftime(date, lag(date), units = "secs") > min.dt) %>%
    mutate(keep = ifelse(is.na(keep), TRUE, keep)) %>%
    arrange(order(date)) %>%
    mutate(lc = factor(lc, levels = c(3,2,1,0,"A","B","Z"), ordered = TRUE)) %>%
    mutate(lon = ifelse(lon > 180, 180 - lon, lon))

  ## reproject from longlat to mercator x,y (km)
  prj <- "+proj=merc +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"
  d[, c("x", "y")] <- as_tibble(project(as.matrix(d[, c("lon", "lat")]), proj = prj))

  if(data.type == "KF") {
    ##  flag any records with smaj/smin = 0
    d <- d %>%
      mutate(keep = ifelse(smaj == 0 | smin == 0, FALSE, keep))
  }

  f <- sum(!d$keep)
#  cat(sprintf("%d observations with duplicate dates will be ignored \n", f))

  ## Prepare operations specific to data type
  switch(data.type,
         LS = {
           ## Argos error multiplication factors
           amf <- data.frame(
             lc = factor(
               c("3", "2", "1", "0", "A", "B", "Z"),
               levels = c("3", "2", "1", "0", "A", "B", "Z"),
               ordered = TRUE
             ),
             amf_x = c(1, 1.54, 3.72, 23.9, 13.51, 44.22, 44.22),
             amf_y = c(1, 1.29, 2.55, 103.7, 14.99, 32.53, 32.53)
           )
           d <- d %>%
             left_join(., amf, by = "lc")
           if(sum(is.na(d$lc)) > 0) stop("\n NA's found in location class values,\n
                                         perhaps your input lc's != c(3,2,1,0,`A`,`B`,`Z`)?")
         },
         KF = {
           ##  convert error ellipse axes from m to km
           ##  convert ellipse orientation from deg to rad
           d <- d %>%
             mutate(smaj = smaj / 1000, smin = smin / 1000) %>%
             mutate(eor = eor / 180 * pi)
         })


  ## use loess smooth to identify outliers by distance (faster than speed filters)
  lo.x <-
    loess(
      x ~ as.numeric(date),
      data = d,
      span = span,
      na.action = "na.exclude",
      control = loess.control(surface = "direct")
    )
  res.x <- lo.x$residuals

  lo.y <-
    loess(
      y ~ as.numeric(date),
      data = d,
      span = span,
      na.action = "na.exclude",
      control = loess.control(surface = "direct")
    )
  res.y <- lo.y$residuals

 ## flag extreme outlier locations, if residuals > min.dist, to be ignored by SSM filter
  d <- d %>%
    mutate(keep = ifelse((
      (res.x <= quantile(res.x, 0.01) |
        res.x >= quantile(res.x, 0.99)) & abs(res.x) > min.dist
    ) |
      (
        (res.y <= quantile(res.y, 0.01) |
          res.y >= quantile(res.y, 0.99)) &
          abs(res.y) > min.dist
      ),
    FALSE,
    keep
    ))

  f1 <- sum(!d$keep) - f
#  cat(sprintf("%d potential outlier locations with residuals > %d km will be ignored \n", f1, min.dist))

  switch(data.type,
         LS = {
           d <- d %>% select(id, date, lc, lon, lat, x, y, amf_x, amf_y, keep)
         },
         KF = {
           d <- d %>% select(id, date, lc, lon, lat, smaj, smin, eor, x, y, keep)
         })
  class(d) <- append("ctrwData", class(d))
#  cat("Data is of class: ", class(d)[1], "  ", class(d)[2], sep = "")

  return(d)
}
