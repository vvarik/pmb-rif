#' Geometric sequence
#'
#' Create geometric sequence
#' @param start the first term in sequence
#' @param end the last term in sequence
#' @param ratio the common ratio (i.e. the multiple that defines the sequence)
#' @examples
#' gseq(1, 16, 2)
#' @export
gseq = function(start, end, ratio){
  power = 0
  result = c()
  more = TRUE
  while (more) {
    result = c(result, start*ratio**power)
    power = power + 1
    more = max(result) < end
  }
  result[result < end]
}

