#' Confidence interval
#'
#' Find the confidence interval. Missing values are removed by default.
#' @param x A vector
#' @param myInt A confidence interval, defaults to 0.975 i.e. 95% confidence
#' interval
#' @examples
#' ci(1:10)
#' @export

ci = function(x, myInt=0.975) {
  mean  = mean(x, na.rm=T)
  sd    = sd(x, na.rm=T)
  n     = sum(!is.na(x))
  ci    = stats::qt(myInt, df=n-1)*sd/sqrt(n)
  ci
}

