#' Insert
#'
#' Insert elements into a vector using index
#' @param vec a vector
#' @param element the element to insert 
#' @param idx the index at which to insert 
#' @examples
#' vec = rep(TRUE, 5) 
#' ins(vec, FALSE, 2)
#' ins(vec, FALSE, c(2, 4))
#' @export
ins = function (vec, element, idx) {
  for(i in idx) vec = append(vec, element, i-1)
  vec
}

