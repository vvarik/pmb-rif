#' Geometric standard deviation
#'
#' Find geometric standard deviation as a multiplication factor. Missing values are removed by default.
#' @param x a vector
#' @examples
#' vec = c(1:10, 100)
#' gmu(vec); gsd(vec)
#' @export


gsd = function(x) {
  exp(sd(log(x), na.rm = T))
}
