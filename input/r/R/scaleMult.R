#' Scale multiplicatively
#'
#' Find geometric mean. Missing values are removed by default.
#' @param x a vector
#' @param target.values a vector, can be of length one
#' @export
scaleMult = function(x, target.values){
  
  this.median = median(as.double(x), na.rm=T)
  target.median = median(as.double(target.values), na.rm=T)
  
  if(this.median==0)
    stop("Division by zero")
  
  #this assumes non-zero-centered values
  multiplicative.factor = target.median / this.median
  
  return(as.double(x) * multiplicative.factor)
  
}
