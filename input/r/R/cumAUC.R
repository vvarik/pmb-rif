#' Cumulative Area Under the Curve
#'
#' Find cumulative area under the curve 
#' @param dat 
#' @examples
#' cumAUC(subset(beaver1, day == 346))
#' @export
cumAUC = function (dat, x, y) {
  pars = as.list(match.call()[-1]) %>% lapply(., as.character)

  out = vector(length = nrow(dat))
  for(i in seq_along(dat[, pars$x])) {
    out[i] = DescTools::AUC(dat[, pars$x][1:i], dat[, pars$y][1:i])
  }
  out
}
