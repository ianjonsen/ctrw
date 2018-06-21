##' Visualise ctrw fit to track data
##'
##' @title plot
##' @param m a fitted object of class ctrw
##' @importFrom ggplot2 ggplot geom_point geom_qq aes ggtitle theme_bw
##' @importFrom gridExtra grid.arrange
##' @method plot ctrw
##' @export

plot.ctrw <- function(m)
{
  dat <- subset(m$data, keep)
  ndat <- subset(m$data, !keep)

  p1 <- ggplot() +
    geom_point(data = dat, aes(x, y), shape = 19, col = grey(0.85)) +
    geom_point(data = ndat, aes(x, y), shape = 20, size = 0.75, col = "firebrick") +
    geom_point(data = fit$fitted, aes(x, y), size = 0.2, shape = 20, col = "dodgerblue") +
    theme_bw() +
    ggtitle("Fitted values")

  p2 <- ggplot() +
    geom_point(data = dat, aes(x, y), shape = 19, col = grey(0.85)) +
    geom_point(data = ndat, aes(x, y), shape = 20, size = 0.75, col = "firebrick") +
    geom_point(data = fit$predicted, aes(x, y), size = 0.2, shape = 20, col = "dodgerblue") +
    theme_bw() +
    ggtitle("Predicted values")

  p3 <- ggplot() +
    geom_point(data = dat, aes(date, x), shape = 19, col = grey(0.85)) +
    geom_point(data = ndat, aes(date, x), shape = 20, size = 0.75, col = "firebrick") +
    geom_point(data = fit$fitted, aes(date, x), size = 0.2, shape = 20, col = "dodgerblue") +
    theme_bw()

  p4 <- ggplot() +
    geom_point(data = dat, aes(date, y), shape = 19, col = grey(0.85)) +
    geom_point(data = ndat, aes(date, y), shape = 20, size = 0.75, col = "firebrick") +
    geom_point(data = fit$fitted, aes(date, y), size = 0.2, shape = 20, col = "dodgerblue") +
    theme_bw()

  p5 <- ggplot() +
    geom_qq(data = fit$fitted, aes(sample = x - dat$x), size = 0.5) +
#    geom_qq_line(data = fit$fitted, aes(sample = x - dat$x)) +
    ggtitle("x_hat - x")

  p6 <- ggplot() +
    geom_qq(data = fit$fitted, aes(sample = y - dat$y), size = 0.5) +
#    geom_qq_line(data = fit$fitted, aes(sample = y - dat$y)) +
    ggtitle("y_hat - y")


  grid.arrange(p1, p2, p3, p4, p5, p6, layout_matrix = matrix(
    c(1, 2, 1, 2, 3, 3, 4, 4, 5, 6),
    nrow = 5,
    ncol = 2,
    byrow = T
    ))
}
