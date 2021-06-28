#' Graphical parameters for base plot
#'
#' Find cumulative area under the curve 
#' @export
setPar = function(...) {
# my favourite parameters to use
  par(cex.axis=1.5, cex.lab=1.5, cex.main=1.5,
    tcl=-0.3, lwd=2, las=1, mfrow=c(2,3),
    #around the plot, below, left, top, right
    mar=c(4.1, 4.1, 1.5, 0.5),  
    # closeness of axis title, label, line
    mgp=c(2.5, 0.75, 0), ...
  )
}
