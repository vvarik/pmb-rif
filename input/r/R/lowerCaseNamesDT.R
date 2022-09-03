#' Lower case names in data.table
#'
#' Set names to lowercase in data.table
#' @param x a data.table
#' @export

lowerCaseNamesDT = function(x) setnames(x, tolower(names(x)))
