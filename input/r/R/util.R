#' Utility functions
#'

#' @export
scale_colour_discrete = function(...) {
  scale_colour_manual(..., values = cbPalette)
}

#' @export
scale_fill_discrete = function(...) {
  scale_fill_manual(..., values = cbPalette)
}


#' @export
getReady = function () {
  
  theme_set(theme_cowplot())
  
  options(datatable.print.nrows=20)
  
  theme_set(theme_pdf())
  update_geom_defaults("point", list(size = 3.5))
  update_geom_defaults("line", list(size = 1))
  theme_update(legend.position = "right")
  .text_size = 6
  .point_size = 4
}


#' @export
contour.ResponseSurface = function (x, ...) 
{
    if (!exists("maxR", x)) 
        stop("maxR statistics were not found.")
    cpdNames = if (!is.null(x$names)) 
        x$names
    else c("Compound 1", "Compound 2")
    args = list(...)
    if (!exists("xlab", args)) 
        args$xlab = paste0("Dose (", cpdNames[[1]], ")")
    if (!exists("ylab", args)) 
        args$ylab = paste0("Dose (", cpdNames[[2]], ")")
    if (!exists("colorPalette", args)) {
        args$colorPalette = c("red", "white", "blue")
        if (x$fitResult$coef["b"] >= x$fitResult$coef["m1"] && 
            x$fitResult$coef["b"] >= x$fitResult$coef["m2"]) {
            args$colorPalette = rev(args$colorPalette)
        }
    }
    args$x = x$maxR
    do.call(pltMaxR, args)
}


#' @export
pltMaxR = function (x, main = "Contour plot for maxR", 
  xlab = "Dose (Compound 1)", ylab = "Dose (Compound 2)", 
  show.magnitude = F,
  rev.axes = F,
  extend.axes = F,
  colorPalette = c("blue", "white", "red"), 
  logScale = TRUE,  
  zTransform = function(z) {z}, 
  #plevels = c(0.7, 0.8, 0.9, 0.95, 0.99, 0.999),
  plevels = c(0.95, 0.99, 0.999, 0.9999),
  cutoff = max(plevels), maxshow = NULL, ...) 
{
    uniqueDoses <- with(x$Ymean, list(d1 = sort(unique(d1)), 
        d2 = sort(unique(d2))))
    doseGrid <- expand.grid(uniqueDoses)
    log10T <- function(z) log10(z + 0.5 * min(z[z != 0]))
    transformF <- if (logScale) log10T else function(z) z
    maxRvalues <- x$Ymean$R
    maxRvalues <- maxRvalues[match(with(doseGrid, paste(d1, d2)), 
        with(x$Ymean, paste(d1, d2)))]
    maxRvalues[is.na(maxRvalues)] <- 0

    if (is.null(maxshow)) {
        maxshow <- quantile(attr(x$Ymean, "distr"), cutoff)
        if (is.null(maxshow)) {
            maxshow <- 3.5
            warning("No `maxshow` parameter specified, so 3.5 is used.")
        }
    }
    origMaxRvalues <- maxRvalues
    maxRvalues[maxRvalues > maxshow] <- maxshow
    maxRvalues[maxRvalues < -maxshow] <- -maxshow
    zlevels <- sort(
      quantile(attr(x$Ymean, "distr"), 
        c(plevels[plevels < cutoff], cutoff)
      ) %o% c(-1, 1))

    pointscale = 1.5 

    if(show.magnitude) {
      pointsize = matrix(
        abs(origMaxRvalues)/pointscale, sapply(uniqueDoses, length)
      )
      pointsize[pointsize < 2] = 1
    } else {
      pointsize = 1.5
    }

    xlim = transformF(uniqueDoses$d1)
    ylim = transformF(uniqueDoses$d2)
    if(extend.axes) {
      xlim = extendrange(xlim)
      ylim = extendrange(ylim)
    } else {
      xlim = range(xlim)
      ylim = range(ylim)
    }

    if(rev.axes) {
      xlim = rev(xlim)
      ylim = rev(ylim)
    }

    filled.contour(
        x = transformF(uniqueDoses$d1),
        y = transformF(uniqueDoses$d2), 
        z = matrix(maxRvalues, sapply(uniqueDoses, length)), 
        levels = zlevels,
        color.palette = colorRampPalette(colorPalette), 
        #frame.plot = F,
        plot.axes = {
            axis(1, at = transformF(uniqueDoses$d1), tck=-0.02,
                mgp = c(2.5, 0.75, 0), 
                labels = format(uniqueDoses$d1, ...), cex.axis = 1.25
            )
            axis(2, at = transformF(uniqueDoses$d2), tck=-0.02,
                mgp = c(2.5, 0.75, 0),
                labels = format(uniqueDoses$d2, ...), cex.axis = 1.25
            )
            points(
                expand.grid(x = transformF(uniqueDoses$d1), 
                y = transformF(uniqueDoses$d2)
                ), 
                col = rgb(0, 0, 0, 0.3), 
                cex = pointsize, pch = 20
               #col = rgb(0, 0, 0, pointalpha), cex = 2
            )
        }, 
        key.axes = {
            axis(4, at = c(0, zlevels), tck = 0, #-0.02, 
              mgp = c(3, 0.3, 0), cex.axis = 1, 
              labels = c(paste0("↑", 
                if (colorPalette[1] == "blue") "Syn" else "Ant", 
                "\n\n \n\n", "↓", 
                if (colorPalette[1] == "blue") "Ant" else "Syn"), 
                paste0(" ", " "),
                1 - rev(plevels[plevels < cutoff]), 
                1 - plevels[plevels < cutoff], 
                paste0(" ", " ")
              )
            )

        }, 
        cex.main = 1.5,
        cex.lab = 1.25,
        key.title = title(main = "p-value", line = 1, cex.main = 1.25), 
        # xlim = rev(range(transformF(uniqueDoses$d1))), 
        # ylim = rev(range(transformF(uniqueDoses$d2))), 
        xlim = xlim,
        ylim = ylim,
        zlim = maxshow * c(-1, 1), main = main, xlab = xlab, 
        ylab = ylab)

      if(show.magnitude) {
        legend('topleft', inset = 0.01,
        legend=c("0-2", "4", "6"), #bty='n',
        box.lwd = 0, bg = 'white',
        title='Z-score', 
        pch=20, 
        pt.cex = c(1, c(4, 6)/pointscale), 
        col = rgb(0, 0, 0, 0.3))
      }

}


#' @export
pltDDI = function(x, ...) {
  pal = RColorBrewer::brewer.pal(9, 'RdBu')
  pal[5] = "#FFFFFF"
  contour(x, xlab = 'PMB, xMIC', ylab = 'RIF, xMIC', 
    colorPalette = pal,
    xlim = c(0.01, 30), ylim = c(0.01, 30), ...)
}

#' @export
getCFUt0 = function(dat) {
  cfu_t0 = dat[dat$time_h==0, "cfu"]*10^6 # cfu's are recorded cfu/10^6
  exp(mean( log(cfu_t0) ) )
}


#' @export
getDetLim = function(dat, vol=5) {
  # finds detection limit, returns for effect and for cfu
  cfu_t0 = findCFUt0(dat)
  det_lim_effect = log10(1*1000/vol) - log10(cfu_t0)
  det_lim_cfu = log10(cfu_t0) + det_lim_effect
  list(effect=det_lim_effect, cfu=det_lim_cfu)
}
