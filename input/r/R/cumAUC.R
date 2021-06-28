#' Cumulative Area Under the Curve
#'
#' Find cumulative area under the curve 
#' @param dat 
#' @examples
#' cumAUC(subset(beaver1, day == 346))
#' @export
cumAUC = function (x, y) {
  out = vector(length = length(x))
  for(i in seq_along(x)) {
    out[i] = DescTools::AUC(x[1:i], y[1:i])
  }
  out
}
