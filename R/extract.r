##' @importFrom geosphere mercator
##' @importFrom dplyr %>% tbl_df arrange mutate select
##' @importFrom tibble as_tibble
##' @export
extract <- function(x, what = "fitted", ...) {
  if (length(list(...)) > 0) {
    warning("additional arguments ignored")
  }

  if(!what %in% c("fitted","predicted","data"))
    stop("Only `fitted`, `predicted` or `data` objects can be extracted")

  switch(what,
         fitted = {
           lapply(x$ssm, function(.) .$fitted) %>%
             do.call(rbind, .) %>%
             as_tibble() %>%
             arrange(id) %>%
             select(id, date, lon, lat, x, y, x.se, y.se)
         },
         predicted = {
           lapply(x$ssm, function(.) .$predicted) %>%
             do.call(rbind, .) %>%
             as_tibble() %>%
             arrange(id) %>%
             select(id, date, lon, lat, x, y, x.se, y.se)
         },
         data = {
           lapply(x$ssm, function(.) .$data) %>%
             do.call(rbind, .) %>%
             as_tibble() %>%
             arrange(id)
         })

}
