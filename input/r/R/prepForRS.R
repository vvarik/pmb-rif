#' Prepare for response surface analysis
#' 
#' Take the variables needed and format otherwise
#' @export
prepForRS = function (x, swapnames = F) {
  x = x[, .(mut, effect, cond, d1, d2, ab1 = abb1, ab2 = abb2)]
  
  if(swapnames) {
    old = c('d1', 'd2', 'ab1', 'ab2')
    new = c('d2', 'd1', 'ab2', 'ab1')
    setnames(x, old, new)
  } 

  # some point down the road handles only df's without complaint
  x = as.data.frame(x) 

  return(x)
}
