##' @importFrom geosphere mercator
##' @importFrom dplyr %>% tbl_df arrange mutate select
##' @export
extract <- function(x, what = "fitted", ...) {
  if (length(list(...)) > 0) {
    warning("additional arguments ignored")
  }

  if(!what %in% c("fitted","predicted","data")) stop("Only `fitted`, `predicted` or `data` objects can be extracted")

  switch(what,
         fitted = {
           lapply(x$ssm, function(.) .$fitted) %>%
             do.call(rbind, .) %>%
             tbl_df() %>%
             arrange(id) %>%
             mutate(lon = mercator(cbind(x,y), r = 6378.137, inverse = TRUE)[,1],
                    lat = mercator(cbind(x,y), r = 6378.137, inverse = TRUE)[,2]) %>%
             select(id, date, lon, lat, x, y, x.se, y.se)
         },
         predicted = {
           lapply(x$ssm, function(.) .$predicted) %>%
             do.call(rbind, .) %>%
             tbl_df() %>%
             arrange(id) %>%
             mutate(lon = mercator(cbind(x,y), r = 6378.137, inverse = TRUE)[,1],
                    lat = mercator(cbind(x,y), r = 6378.137, inverse = TRUE)[,2]) %>%
             select(id, date, lon, lat, x, y, x.se, y.se)
         },
         data = {
           lapply(x$ssm, function(.) .$data) %>%
             do.call(rbind, .) %>%
             tbl_df() %>%
             arrange(id)
         })

}
