#' Geometric mean
#'
#' Find geometric mean. Missing values are removed by default.
#' @param x a vector
#' @param trim the fraction (0 to 0.5) of observations to be trimmed from each
#' end of ‘x’ before the mean is computed. Values of trim outside that range are
#' taken as the nearest endpoint.
#' @examples
#' vec = c(1:10, 100)
#' mean(vec); gmu(vec)
#' @export

gmu = function(x, trim=0) {
  exp(mean(log(x), na.rm=T, trim=trim))
}

