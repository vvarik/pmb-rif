#' Geometric confidence interval
#'
#' Find geometric confidence interval as a multiplication factor. Missing values
#' are removed by default. Currently, defaults to 95% confidence level.
#' @param x a vector
#' @examples
#' vec = c(1:10, 100)
#' gmu(vec); gci(vec)
#' @export

gci = function(x) {
  exp(ci(log(x)))
}
