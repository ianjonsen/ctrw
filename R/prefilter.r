##' @title Prepare Argos data for fitting a ctrw model
##'
##' @description \code{prefilter()} (1) determines Argos data type (LS or KF); (2) converts dates to POSIXt;
##' identifies observations with duplicate dates; (3) orders observations in time;
##' (4) removes duplicate observations; (5) removes observations occurring within 60 s of one another (keeps first);
##' (6) shifts longitudes that straddle -180,180 to 0,360 and vice-versa; (7) projects lonlat coords to mercator x,y
##' coords (in km); (8) adds location error multiplication factors based on Argos location class (for type LS);
##' and (9) uses a loess smooth to identify potential outlier locations (by distance only) to be ignored when fitting
##' the \code{ctrw} model
##'
##' @details Internal function
##'
##' @param d input data - must have 5 (LS), or 8 (KF) columns (see details)
##' @param vmax max travel rate (m/s) passed to argosfilter::sdafilter to define outlier locations
##' @param min.dt minimum allowable time difference between observations; dt < min.dt will be ignored by the SSM
##' @param time.gap not currently implemented
##' @importFrom lubridate ymd_hms
##' @importFrom stats loess
##' @importFrom dplyr mutate distinct arrange filter select %>% left_join lag
##' @importFrom rgdal project
##' @importFrom argosfilter sdafilter
##'
##' @export

prefilter <- function(d, vmax = 10, min.dt = 60, time.gap = NULL) {

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
  ##  if lon spans -180,180 then shift to
  ##    0,360; else if lon spans 360,0 then shift to
  ##    -180,180

  d <- d %>%
    mutate(date = ymd_hms(date, tz = "GMT")) %>%
    mutate(keep = difftime(date, lag(date), units = "secs") > min.dt) %>%
    mutate(keep = ifelse(is.na(keep), TRUE, keep)) %>%
    arrange(order(date)) %>%
    mutate(lc = factor(lc, levels = c(3,2,1,0,"A","B","Z"), ordered = TRUE))

  ## Use argosfilter::sdafilter to identify outlier locations
  filt <- sdafilter(d$lat, d$lon, d$date, d$lc, ang=-1, vmax=vmax)
  d <- d %>%
    mutate(keep = ifelse(filt == "removed", FALSE, keep))

  if(min(d$lon, na.rm = TRUE) < 0 & diff(range(d$lon, na.rm = TRUE)) > 350) {
    d <- d %>%
      mutate(lon = wrap_lon(lon, 0)) %>%
      mutate(cntr = 0)
  } else if (min(d$lon) < 0 & max(d$lon) > 0){
    d <- d %>%
      mutate(lon = wrap_lon(lon, -180)) %>%
      mutate(cntr = 180)
  } else {
    d <- d %>%
      mutate(cntr = 90)
  }

  ## reproject from longlat to mercator x,y (km)
  if(d$cntr[1] == 0 | d$cntr[1] == 90){
    prj <- "+proj=merc +lat_0=0 +lon_0=180 +datum=WGS84 +units=km +no_defs"
  } else if(d$cntr[1] == 180) {
    prj <- "+proj=merc +lat_0=0 +lon_0=0 +datum=WGS84 +units=km +no_defs"
  }
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

  f1 <- sum(!d$keep) - f
#  cat(sprintf("%d potential outlier locations with residuals > %d km will be ignored \n", f1, min.dist))

  switch(data.type,
         LS = {
           d <- d %>% select(id, date, lc, lon, lat, x, y, amf_x, amf_y, keep, cntr)
         },
         KF = {
           d <- d %>% select(id, date, lc, lon, lat, smaj, smin, eor, x, y, keep, cntr)
         })
  class(d) <- append("ctrwData", class(d))
#  cat("Data is of class: ", class(d)[1], "  ", class(d)[2], sep = "")

  return(d)
}
