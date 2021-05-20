#' Geometric summaries
#'
#' Derive some geometric summaries
#' @param x A vector
#' @param dec Decimal place to round to
#' @return  'gsum' returns geometric versions of mean ('gmu'), standard
#' deviation ('gsd'), standard error ('gse'), and 95\% confidence interval
#' ('gci').
#' @details Geometric estimates of spread are factors of
#' multiplication/division.
#' @examples
#' gsum(c(1:10, 100))
#' @export

gsum = function (x, dec=3) {
  mu = gmu(x)
  sd = gsd(x)
  se = gse(x)
  ci = gci(x)
  out = list(gmu = mu, gsd = sd, gse = se, gci = ci)
  lapply(out, round, dec)
}
