##' Prefilter
##'
##' Prepare Argos data for fitting a ctrw model
##'
##' @importFrom lubridate ymd_hms
##' @importFrom dplyr mutate distinct arrange filter select %>%
##'
##' @export

prefilter <- function(d, min.dist = 100, time.gap = NULL) {

  ## Argos error multiplication factors
  amf <- data.frame(
    lc = factor(
      c("3", "2", "1", "0", "A", "B"),
      levels = c("3", "2", "1", "0", "A", "B"),
      ordered = TRUE
    ),
    amf_x = c(1, 1.54, 3.72, 23.9, 13.51, 44.22),
    amf_y = c(1, 1.29, 2.55, 103.7, 14.99, 32.53)
  )

  ##  convert dates to POSIXt
  ##  flag any duplicate date records,
  ##  order records by time,
  ##  remove Z-class locations
  ##  set lc to ordered factor
  ##  convert lon from 0,360 to -180,180
  ##  reproject lon,lat to x,y in km
  ##  convert error ellipse axes from m to km
  ##  convert ellipse orientation from deg to rad
  ##  add Argos error multiplication factors
  d <- d %>%
    mutate(date = ymd_hms(date, tz = "GMT")) %>%
    mutate(keep = difftime(date, lag(date), units = "secs") > 10) %>%
    arrange(order(date)) %>%
    mutate(keep = ifelse(lc == "Z", FALSE, keep)) %>%
    mutate(lc = factor(lc, levels = c(3,2,1,0,"A","B"), ordered = TRUE)) %>%
    mutate(lon = ifelse(lon > 180, 180 - lon, lon)) %>%
    mutate(x = geosphere::mercator(cbind(.$lon,.$lat), r = 6378.137)[,1]) %>%
    mutate(y = geosphere::mercator(cbind(.$lon,.$lat), r = 6378.137)[,2]) %>%
    mutate(smaj = smaj / 1000, smin = smin / 1000) %>%
    mutate(eor = eor / 180 * pi) %>%
    left_join(., amf, by = "lc")

  ## use loess smooth to identify outliers by distance (faster than speed filters)
  lo.x <-
    loess(
      x ~ as.numeric(date),
      data = d,
      span = 0.01,
      na.action = "na.exclude",
      control = loess.control(surface = "direct")
    )
  res.x <- lo.x$residuals

  lo.y <-
    loess(
      y ~ as.numeric(date),
      data = d,
      span = 0.01,
      na.action = "na.exclude",
      control = loess.control(surface = "direct")
    )
  res.y <- lo.y$residuals

 ## flag extreme outlier locations to be ignored by SSM filter
  d <- d %>%
    mutate(keep = ifelse((
      res.x <= quantile(res.x, 0.001) &
        res.x >= quantile(res.x, 0.999) & abs(res.x) >= min.dist
    ) &
      (
        res.y <= quantile(res.y, 0.001) &
          res.y >= quantile(res.y, 0.999) &
          abs(res.y) >= min.dist
      ),
    FALSE,
    keep
    ))

 d %>% select(id, date, lc, lon, lat, smaj, smin, eor, x, y, amf_x, amf_y, keep)
}
