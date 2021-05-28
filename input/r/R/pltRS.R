#' Plot response surface
#'
#' domesticated BIGL::plot.ResponseSurface
#' note that it calls a subfunction pltResponseSurface (code below)
#'
#' @export
pltRS <- function (x, color = c("z-score", "maxR", "occupancy", "confint"), ...) 
{
  color <- match.arg(color)
  inputs <- as.list(substitute(list(...)))[-1L]
  if (!exists("colorPalette", inputs)) {
    inputs$colorPalette <- c("#2166AC", rep("grey70", 2), "#B2182B")
    #if (x$fitResult$coef["b"] > x$fitResult$coef["m1"]) {
    if (x$fitResult$coef["b"] >= x$fitResult$coef["m1"] && 
        x$fitResult$coef["b"] >= x$fitResult$coef["m2"]) {
      inputs$colorPalette <- rev(inputs$colorPalette)
    }
  }
  if (color == "z-score") {
    boundary <- sd(x$offAxisTable[["z.score"]])
    inputs$colorBy <- x$offAxisTable[, c("d1", "d2", "z.score")]
    if (!exists("breaks", inputs)) 
      inputs$breaks <- c(-Inf, -boundary, 0, boundary, 
                         Inf)
    if (!exists("main", inputs)) 
      inputs$main <- "Z-scores"
  }
  else if (color == "maxR") {
    inputs$colorBy <- x$maxR$Ymean[, c("d1", "d2", "R")]
    q <- attr(x$maxR$Ymean, "q")
    if (!exists("breaks", inputs)) 
      inputs$breaks <- c(-Inf, -q, 0, q, Inf)
    if (!exists("main", inputs)) 
      inputs$main <- "maxR"
  }
  else if (color == "occupancy") {
    inputs$colorBy <- x$occupancy
    inputs$colorPalette <- c("#EFF3FF", "#BDD7E7", "#6BAED6", 
                             "#2171B5")
    if (!exists("breaks", inputs)) 
      inputs$breaks <- c(0, 0.25, 0.5, 0.75, 1)
    if (!exists("main", inputs)) 
      inputs$main <- "Occupancy rate"
  }
  else if (color == "confInt") {
      if (is.null(x$confInt)) 
          stop("No confidence intervals were calculated")
      x$confInt$offAxis = cbind(x$confInt$offAxis,
      t(sapply(rownames(x$confInt$offAxis), 
          function(y) strsplit(y, split = "_")[[1]])))
      colnames(x$confInt$offAxis)[5:6] = c("d1", "d2")
      inputs$colorBy = x$confInt$offAxis[, c("d1", "d2", "call")]
      if (!exists("breaks", inputs)) 
          inputs$breaks <- seq_len(4)
      if (!exists("main", inputs)) 
          inputs$main <- "Calls from confidence intervals"
  }
  inputs$data <- x$data
  inputs$fitResult <- x$fitResult
  inputs$transforms <- x$transforms
  inputs$null_model <- x$null_model
  do.call(pltResponseSurface, inputs)
}


#' @export
pltResponseSurface = function (data, fitResult = NULL, 
	transforms = fitResult$transforms, 
    predSurface = NULL, null_model = c("loewe", "hsa", "bliss", "loewe2"), colorPalette = c("blue", 
        "grey70", "red"), colorBy = "none", colorPoints = c("black", 
        "sandybrown", "brown", "white"), breaks = c(-Inf, 0, 
        Inf), radius = NULL, logScale = TRUE, colorfun = median, 
    zTransform = function(x) x, add = FALSE, main = "", legend = TRUE, 
    lit=T, zlab="Response", xlab="Cpd1", ylab="Cpd2", zlim = NULL,
    xat = "pretty", yat = "pretty", plotfun = NULL, ...) 
{
    null_model <- match.arg(null_model)
    if (missing(fitResult) & missing(predSurface)) 
        stop("Marginals fit result or predicted surface need to be supplied.")
    if (is.character(colorBy) & all(colorBy %in% colors())) {
        colorPalette <- colorBy
        colorBy <- "colors"
    }
    uniqueDoses <- with(data, list(d1 = sort(unique(d1)), d2 = sort(unique(d2))))
    doseGrid <- expand.grid(uniqueDoses)
    logT <- function(z) log(z + 0.5 * min(z[z != 0]))
    log10T <- function(z) log10(z + 0.5 * min(z[z != 0]))
    transformF <- if (logScale) 
        log10T
    else function(z) z
    zGrid <- predSurface
   if (is.null(predSurface)) {
        respSurface <- BIGL:::predictResponseSurface(doseGrid, fitResult, 
            null_model = null_model, transforms = transforms)
        if (!is.null(transforms)) {
            predSurface <- with(transforms, InvPowerT(respSurface, 
                compositeArgs))
        }
        else {
            predSurface <- respSurface
        }
        zGrid <- predSurface
    }
    if (missing(radius)) {
        avgEffect <- abs(mean(zTransform(zGrid)))
        radius <- max(10^(log10(avgEffect) - 1.5), 0.05)
    }
    if (inherits(colorBy, c("matrix", "data.frame"))) {
        stopifnot(c("d1", "d2") %in% colnames(colorBy))
        colorVec <- colorBy
        colorBy <- "asis"
        coloredBy <- colorVec
        cols <- colnames(coloredBy)
        pCols <- which(!(cols %in% c("d1", "d2")))
        if (length(pCols) > 1) 
            pCols <- min[pCols]
        if (any(duplicated(coloredBy[, c("d1", "d2")]))) {
            coloredBy <- aggregate(coloredBy[, ..pCols], by = coloredBy[, 
                c("d1", "d2")], FUN = colorfun)
        }
        # cols <- colnames(coloredBy)
        # pCols <- which(!(cols %in% c("d1", "d2")))
        if (length(pCols) > 1) 
            pCols <- min[pCols]
        colorVec <- rep(NA, nrow(doseGrid))
        for (i in 1L:nrow(doseGrid)) {
            ind <- match(paste(doseGrid[i, ], collapse = ";"), 
                apply(coloredBy[, c("d1", "d2")], 1, function(x) paste(x, 
                  collapse = ";")))
            if (!is.na(ind)) 
                colorVec[i] <- coloredBy[[pCols]][ind]
        }
    }
    dataOffAxis <- with(data, data[d1 & d2, , drop = FALSE])
    predOffAxis <- predSurface[cbind(match(dataOffAxis$d1, uniqueDoses$d1), 
        match(dataOffAxis$d2, uniqueDoses$d2))]
    if (nrow(dataOffAxis) == 0) {
        warning("No off-axis observations were found. Surface won't be custom colored..")
        colorBy <- "none"
    }
    surfaceColors <- colorRampPalette(colorPalette)(length(breaks) - 
        1)
    getFF = function(response) {
        if (is.numeric(response)) {
            cut(response, breaks = breaks, include.lowest = TRUE)
        }
        else if (is.factor(response)) {
            response
        }
        else if (is.character(response)) {
            factor(response, levels = c("Syn", "None", "Ant"), 
                labels = c("Syn", "None", "Ant"), ordered = TRUE)
        }
    }
    surfaceColor <- function(response) {
        ff <- getFF(response)
        zcol <- surfaceColors[ff]
        return(zcol)
    }
    getLabels <- function(response) {
        ff <- getFF(response)
        labels <- gsub(",", ", ", levels(ff))
        return(labels)
    }
    if (colorBy == "asis") {
        if (is.numeric(colorVec)) 
            colorVec[is.na(colorVec)] <- 0
        zcol <- surfaceColor(colorVec)
        labels <- getLabels(colorVec)
    }
    else if (colorBy == "colors") {
        zcol <- rep(colorPalette, length(zGrid))
    }
    else {
        zGridFloor <- floor(100 * zGrid)
        col <- terrain.colors(diff(range(zGridFloor, na.rm = TRUE)))
        zcol <- col[zGridFloor - min(zGridFloor, na.rm = TRUE) + 
            1]
    }
    labnames <- c(zlab, xlab, ylab)
    if (!is.null(attr(data, "orig.colnames"))) 
        labnames <- attr(data, "orig.colnames")
    if (!add) {
        if (!is.null(plotfun)) 
            data <- aggregate(effect ~ d1 + d2, data, FUN = plotfun)[, 
                names(data)]
        plot3d(transformF(data$d1), transformF(data$d2), zTransform(data$effect), 
          zlim = zlim,
            xlab = "", ylab = "", zlab = "", box = FALSE, axes = FALSE)
        if (!is.numeric(xat)) {
            xat <- match.arg(xat, c("pretty", "actual"))
            if (xat == "pretty") {
                xlab <- axisTicks(range(transformF(uniqueDoses$d1)), 
                  log = logScale, nint = 3)
                if (logScale && length(xlab) > 4) 
                  xlab <- xlab[!(log10(xlab)%%1)]
            }
            else xlab <- uniqueDoses$d1
            xat <- transformF(xlab)
        }
        else {
            xlab <- xat
            xat <- transformF(xat)
        }
        if (!is.numeric(yat)) {
            yat <- match.arg(yat, c("pretty", "actual"))
            if (yat == "pretty") {
                ylab <- axisTicks(range(transformF(uniqueDoses$d2)), 
                  log = logScale, nint = 3)
                if (logScale && length(ylab) > 4) 
                  ylab <- ylab[!(log10(ylab)%%1)]
            }
            else ylab <- uniqueDoses$d2
            yat <- transformF(ylab)
        }
        else {
            ylab <- yat
            yat <- transformF(yat)
        }
        with(data, spheres3d(transformF(d1), transformF(d2), 
            zTransform(effect), radius = radius, add = TRUE, 
            lit=lit,
            col = colorPoints[1 + 1 * (d2 == 0) + 2 * (d1 == 
                0)]))
    }
#    planes3d(0, 0, 1, zTransform(0), col = "grey", lit = FALSE, 
#        alpha = 0.3)
    persp3d(transformF(uniqueDoses$d1), transformF(uniqueDoses$d2), 
        zTransform(zGrid), add = TRUE, col = zcol, alpha = 0.6, 
        lit = FALSE, aspect = FALSE, lwd = 2)
    bgplot3d({
        plot.new()
        title(main = main)
        if (legend && is.character(labels)) 
            legend("topleft", legend = labels, pch = 16, cex = 1.5, 
                col = surfaceColors)
    })
    if (!add) {
        rgl.bbox(xlen = 0, ylen = 0, zlen = 0, expand = 1.03, 
            color = "#000000", front = "lines", back = "cull")
        rgl.bbox(color=c("black", "black"),          # grey60 surface and black text
                 emission="grey50",       # emission color is grey50
                 marklen = 30,
                 #xunit = 'pretty', yunit = 'pretty', zunit='pretty',
                 front = 'lines',
                 back = 'lines',
                 xlen = 0, ylen = 0, zlen = 0,
                 # xat = xat, xlab = xlab, 
                 # yat = yat, ylab = ylab
                 expand = 1.03
                 )
        axis3d(edge = "x+-", at = xat, labels = xlab)
        axis3d(edge = "y+-", at = yat, labels = ylab)
        axis3d(edge = "z+-")
        mtext3d(labnames[2], edge = "x+-", line = 2)
        mtext3d(labnames[3], edge = "y+-", line = 2)
        mtext3d(labnames[1], edge = "z+-", line = 2)
    }
    persp3d(transformF(uniqueDoses$d1), transformF(uniqueDoses$d2), 
        zTransform(zGrid), add = TRUE, col = zcol, alpha = 0.6, 
        lit = FALSE, aspect = FALSE, lwd = 2)
}
