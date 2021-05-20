#' Turn strings into CamelCase
#'
#' Often, the variables in data come in different forms. If needed, camel case 
#' to uniform them.
#' @param s a vector of strings
#         strict a boolean, if TRUE will capitalize also
#         capitalized abbreviations
#' @examples
#' camelCase('this_is_fun')
#' @export

camelCase = function(s, strict = FALSE) {
  capWords = function(x, ...) {
    cap = function(x) {
      paste(
        toupper(substring(x, 1, 1)), 
        {x = substring(x, 2); if(strict) tolower(x) else x}, 
        sep = "", collapse = " " 
      )
    }

    sapply(strsplit(x, split = " "), cap, USE.NAMES = !is.null(names(x)))
  }

  out = gsub('_', ' ', s)
  out = sapply(out, capWords)
  gsub(' ', '', out)
}
