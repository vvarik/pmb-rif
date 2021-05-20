#' Geometric standard error
#'
#' Find geometric standard error as a multiplication factor. Missing values are removed by default.
#' @param x a vector
#' @examples
#' vec = c(1:10, 100)
#' gmu(vec); gse(vec)
#' @export


gse = function(x) {
  exp(se(log(x)))
}

