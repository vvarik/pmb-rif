#' Standard error
#'
#' Find the standard error of the mean. Missing values are removed by default.
#' @param x a vector
#' @examples
#' se(1:10)#' 
#' @export
 
se = function(x) {
  stats::sd(x, na.rm=T)/sqrt(sum(!is.na(x)))
}
