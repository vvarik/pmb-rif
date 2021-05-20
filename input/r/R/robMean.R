#' Robust mean
#'
#' Find the robust mean
#' @param x A vector
#' @examples
#' robMean(1:10)
#' @export

robMean = function(x){
  # smhuber does not like 0, probably has to do with 'tol'
  x = x[x!=0]  
  # if zeros are all there is, return zero
  if(length(x)==0) out = 0  
  # obtain robust mean
  else if (length(x) > 1) out = smoothmest::smhuber(x)$mu 
  # if what's left is just one nonzero value
  else out = x  
  return(out)
}
