#' Fit surface model
#'
#' Convenience function to fit with most common parameters
#' Skips to next step in case of an error
#'
#' @export
fitSurf = function(x, dr, model='loewe2', B=B.CP) {
  try(fitSurface(x, dr, null_model=model, statistic = "both",
      B.CP = B, parallel = FALSE))
}
