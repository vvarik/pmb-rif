#' Customized sum 
#'
#' Set na.rm to default and sum NAs to 0
#' @param x A vector
#' @export

mySum = function(x) {
  if(all(is.na(x))) 0
  else sum(x, na.rm=TRUE)
}
