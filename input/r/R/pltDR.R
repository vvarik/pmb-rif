#' Plot dose-response for monotherapies
#'
#' @param formula A formula for fit
#' @param dat Data in data.table
#' @param col The color of lines, points
#' @param by A grouping variable for fit
#' @param type A type of drc plot, one of 'average', 'all', 'none',
#' 'confidence', 'bars', 'obs'.
#' @export
pltDR = function(formula, dat, col='red3', by=cond, type=average, ci95=T, 
  zero = 10, ...) {
  
  # get the arguments ready for fit
  type = deparse(substitute(type))
  
  # dummy variable for drm
  dat$CURVE = factor(dat[, eval(substitute(by))], levels = c('intra', 'pH5.5', 'pH7.4'))

  # fit the 4PL model 
  fit = drm(formula, curveid = CURVE, data = subset(dat),
     fct=LL2.4(names=c("hill", "Emax", "Emin", "logEC50"),
               fixed=c(1, NA, NA, NA)),
     pmodels = data.frame(CURVE, CURVE, CURVE)
  )

  par(lwd=1.5, cex=1.25, tcl=-0.25, #pty='s', # asp = 1
            mar=c(3.1, 4.1, 0.5, 0.5), # around the plot, below, left, top, right
            mgp=c(1.7, 0.5, 0)  # closeness of axis title, label, line
  )
  
  plot(fit, type = type, ...,
    pch=c(19, 21, 0), col=col, lty=1, broken = T,
    ylab = expression(log[10]*'('*over('CFU'['24h'], 'CFU'['0h'])*')'),
    ylim = c(-6, 5)
  )


  # confidence intervals; the original ones from drm do not always work
  addCI = function(fit) {
    for (i in unique(fit$origData$cond)) {

      vec = setNames(
        c(with(fit$dataList, range(dose[dose != 0])), 1.05),
        c('start', 'end', 'ratio')
      )
      vec[1] = vec[1]/zero

      tmp = do.call(gseq, as.list(vec))

      preds = predict(fit, newdata = data.frame(tmp, CURVE=i), 
                       interval = 'confidence')
      polygon(c(rev(tmp), tmp), c(rev(preds[ ,3]), preds[ ,2]), 
              col = adjustcolor(col, 0.2), border = NA)
    }
  }

  addCI(fit)
}
