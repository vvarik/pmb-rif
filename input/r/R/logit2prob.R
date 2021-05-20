#' Logit to probability
#'
#' Turn logits into probability
#' @param x A vector
#' @export

logit2prob = function(x, y=0) exp(x + y) / (1 + exp(x + y))
