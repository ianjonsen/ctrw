prefilter <- function(d, time.gap = 21) {
  
  ##  convert dates to POSIXt
  ##  remove any duplicate time records,
  ##  order records by time,
  ##  remove Z-class locations
  ##  remove locations in N Hemisphere
  d <- d %>%
    mutate(date = lubridate::ymd_hms(date, tz = "GMT")) %>%
    distinct(date, .keep_all = TRUE) %>%
    arrange(order(date)) %>%
    filter(lc != "Z") %>%
    filter(lat < 0)

  ## use loess smooth to identify outliers by distance (faster than speed filters)
  res.lon <-
    loess(
      lon ~ as.numeric(date),
      data = d,
      span = 0.1,
      na.action = "na.exclude",
      control = loess.control(surface = "direct")
    ) %>%
    .$residuals 
 
  res.lat <-
    loess(
      lat ~ as.numeric(date),
      data = d,
      span = 0.1,
      na.action = "na.exclude",
      control = loess.control(surface = "direct")
    ) %>%
    .$residuals 
  
 ## flag extreme outlier locations to be ignored by SSM filter
 d <- d %>% 
    mutate(dist.keep = ifelse((res.lon > quantile(res.lon, 0.001) & res.lon < quantile(res.lon, 0.999)) & 
           (res.lat > quantile(res.lat, 0.001) & res.lat < quantile(res.lat, 0.999)), TRUE, FALSE))
 
 ## convert ellipse orientation from deg to rad
 ## convert lon from 0,360 to -180,180
 ## convert coords from lonlat to mercator (km)
 ## convert ellipse axes from m to km
 d <- d %>% 
   mutate(eor = eor / 180 * pi) %>%
   mutate(lon = ifelse(lon > 180, 180 - lon, lon)) 
  
 ## identify time gaps (in days) between observations
 d <- d %>%
   mutate(diff = (date - lag(date)) / 86400) %>%
   mutate(time.keep = TRUE)
 
 
 ## flag segments after long time gaps
 q <- which(d$diff > time.gap)
 if(length(q) > 0) {
   q <- min(q)
   d$time.keep[q:nrow(d)] <- FALSE
 }

 d %>% select(id, date, lc, lon, lat, smaj, smin, eor, dist.keep, time.keep)
}