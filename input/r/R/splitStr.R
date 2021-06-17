#' Split strings
#'
#' Split strings (or vectors of strings) into vectors or extract a single string
#' based on index 
#' @export
splitStr <- function(string.vector, delimiter, index=NA){
  if(is.na(index)){
    return(unlist(strsplit(as.character(string.vector), split=delimiter,fixed=TRUE)))
  } else {
    return(sapply(strsplit(as.character(string.vector), split=delimiter,fixed=TRUE), "[[", index))
  }
}
