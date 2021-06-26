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
getRS = function (dat, name) {
  vec = as.character(unique(dat$strain))
  out = setNames(vector(mode = "list", length = length(vec)), vec)
  
  for(i in vec) {
    cat('*** Took on ', i, '\n')
    out[[i]] = analyzeRS(dat[strain == i], 'cond')
  }

  for(i in seq_along(out))
    for(j in seq_along(out[[i]]))
      names(out[[i]])[j] = unique(out[[i]][[j]]$rsl$data$cond)

  saveRDS(out, paste0('input/dat/tmp/', name, '.rds'))

  out
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


#' @export
remNonCharcoal = function (dat) {
  ## remove non-char results for combination and high rif
  ## duplicates
  idx_dup1 = duplicated(subset(dat, select=c(date, time_h, d_1, d_2, rep)))
  idx_dup2 = duplicated(subset(dat, select=c(date, time_h, d_1, d_2, rep)), fromLast = T)
  idx = idx_dup1 | idx_dup2 
  idx = idx & dat$time_h!=0  # remove time 0, where duplicates were necessary
  idx = idx & (dat$d_1==48 | (dat$d_1==4.8 & dat$d_2==1) )  # take only high rif and combo
  idx = idx & dat$char==0  # take only no-char cases
  ## apply the filter
  dat[!idx, ]
}


#' @export
remBadPlates = function (dat) {
  bad.plates = 
    fread("input/dat/raw/iris_files/qc.csv", skip = 3) %>% 
    filter(is.na(row)) %>% 
    #filter(notes=="remove the whole plate")
    filter(!str_detect(notes, "normalize"))
  
  dat[!filename %in% bad.plates$file]
}


#' @export
flagBadColonies = function (dat) {

  bad.colonies = 
    fread("input/dat/raw/iris_files/qc.csv", skip = 3) %>% 
    filter(!is.na(row)) %>% 
    dplyr::select(filename=file, row, column) %>% 
    mutate(filtered.out=T)

   merge(dat, bad.colonies, by=c("filename", "row", "column"), 
     all.x=T)
  
}


#' @export
turnBadColoniesToNAs = function (dat) {
  dat[filtered.out==T, opacity:=NA]
  dat[filtered.out==T, size:=NA]
}


#' @export
setMisdetectedColoniesToNAs = function (dat) {
  
  # Find the mutants that are absent only a few times (and thus probably a
  # misdetection than an empty spot). Absent mutants have size zero.
  arbitrary.threshold = 0.4
  colonies.sometimes.missing = 
    dat[, .(frequency.absent = nrow(.SD[size==0])/.N), colony] %>% 
    .[frequency.absent > 0 & frequency.absent < arbitrary.threshold] %>% 
    .[, colony]
  
  #set to NA these possible misdetections
  dat[size==0 & colony %in% colonies.sometimes.missing, c('opacity', 'size') := NA ]

}


#' @export
getPairwiseReplicateCorrelations = function(dat) {
  #get all replicate correlations per condition, and per plate number

  correlations = list()

  for(cond in unique(dat$condition)){
    for(plate in unique(dat$plate.number)){
      #get all pairwise combinations of the conditions
      sset = dat[condition==cond][plate.number==plate]
      length(unique(sset$biorep))
  
      if(length(unique(sset$biorep))<2){
        #there's no replicate reproducibility with one replicate
        cat(paste0("Warning: only one replicate for ", cond, " in library plate", plate))
        cat("\n")
        next
      }
      
      rep.pairs = combn(unique(sset$biorep), 2)
      
      for(k in 1:dim(rep.pairs)[2]){
        this.correlation = cor(
          sset[biorep==rep.pairs[1,k]]$opacity,
          sset[biorep==rep.pairs[2,k]]$opacity,
          use="na.or.complete"
        )
        
        correlations = c(correlations, 
          list(data.table(correlation=this.correlation,
            replicate.A = rep.pairs[1,k],
            replicate.B = rep.pairs[2,k],
            filename.A = as.character(unique(sset[biorep==rep.pairs[1,k]]$filename)),
            filename.B = as.character(unique(sset[biorep==rep.pairs[2,k]]$filename)),
            condition = cond,
            plate.number = plate)))
      }
    }
  }
  return(rbindlist(correlations))
}


#' @export
getCandidateFilesToKickOut = function(dat) {
  #Count number of times a specific filename has correlations under a threshold
  #vs over the threshold (all times). Kick out filenames that have 2
  #correlations under 0.5 (disagree with both other replicates)
  bad.files = data.table(table(c(dat[correlation<0.5]$filename.A, 
                                 dat[correlation<0.5]$filename.B)))
  all.files = data.table(table(c(dat$filename.A, 
                                 dat$filename.B)))
  setnames(bad.files, c("filename", "times.discordant"))
  setnames(all.files, c("filename", "times.present"))
  out = merge(bad.files, all.files, by="filename")
  out[,ratio.discordant.vs.all:=times.discordant/times.present]
  return(out)
}


#' @export
correctOpacityForPlates = function (dat) {
  dat[, opacity.raw := opacity]  #store original
  median.to.set = round(median(dat$opacity, na.rm=T), -2)
  dat[,opacity:=scaleMult(x=opacity, target.values=median.to.set), by=filename]
}


#' @export
long2wide = function(dat) {
  #There's 2 LB "conditions", each with 5 replicates as control for each plate
  #position. Take the median of both LB conditions across five replicates.
  control.medians = dat[grepl('noAb', condition), .(
    control.median = as.numeric(median(opacity, na.rm=T)),
    count = sum(!is.na(opacity))
    ), .(media, colony)]
  
  condition.iris.files = dat[!grepl('noAb', condition)]
  
  condition.vs.control = merge(condition.iris.files, control.medians, by=c("media", "colony"))

  return(condition.vs.control)
}


#' @export
getEpsilons = function (dat) {
  out = setkey(dat, NULL)

  out = out[grepl("Rifampicin16Polymyxin2", condition), 
    .(colony, biorep, media, combination.opacity=opacity,
      combination.f=f.condition)] %>% 

    merge(out[antibiotic=="Rifampicin16", 
      .(colony, biorep, media, Rifampicin16.opacity=opacity,
        Rifampicin16.f=f.condition)]) %>% 

    merge(out[antibiotic=="Polymyxin2", 
      .(colony, biorep, media, Polymyxin2.opacity=opacity,
        Polymyxin2.f=f.condition)]) 
  
  
  #this table will bring f.values of drug12, drug1, and drug2 in different
  #columns, so it will be easy to then perform the t-test
  out[,expected.f:=Rifampicin16.f*Polymyxin2.f]
  out[,epsilon:=combination.f-expected.f]
  
  out[,media.mean.epsilon:=mean(epsilon, na.rm=T), by=media]
  
  return(out)
  
}


#' @export
addAnnotation = function (dat) {
  PA14.map = fread('input/dat/raw/PA14_computational_v2_updated_annotations.csv') %>% 
    .[,.(
      plate = `1536 plate`, 
      row = `1536 row`, 
      column = `1536 column`, 
      locus = `""Active"" Gene Locus`, 
      gene.id = `""Active"" GeneID`, 
      gene.name = `""Active"" Gene Name`,
      PAO1.ortholog = `PAO1 ortholog of ""Active"" gene`,
      Tn.pos.bp = `Tn insertion position within ""Active"" Gene (bp)` 
    )] %>% 
    .[,colony := plate*10000 + row*100 + column]
  
  # add information about n of unique Tn mutants per gene
  PA14.map[!is.na(gene.id), Tn.per.gene:=uniqueN(Tn.pos.bp), locus]
  PA14.map[, Tn.mutant.id:=paste(locus, Tn.pos.bp, sep='.')]
  PA14.map[!is.na(gene.id), copy.n:=uniqueN(colony), Tn.mutant.id]
  PA14.map[!is.na(gene.id), copy.id:=seq_len(.N), Tn.mutant.id]

  PA14.map[,gene.name.to.show := ifelse(gene.name=='', gsub('PA14_', '', locus),
    gene.name)]
  
  merge(dat, PA14.map, by='colony', all.x=T)
}

#' @export
tTest = function (dat) {
#compare each gene's epsilon to medium's mean epsilon -- a one-sample t-test

  .tTest = function(x, ...){
    try({
      result = t.test(x, ...)
      return(list(as.numeric(result$statistic), result$p.value))
    },silent=T)
    return(list(NA_real_, NA_real_))
  }

  dat[,c("t.test.statistic", "t.test.pvalue") := .tTest(epsilon,
    mu=unique(media.mean.epsilon)), .(Tn.mutant.id, media)]
  
  dat[,c("t.test.statistic.gene", "t.test.pvalue.gene") := .tTest(epsilon,
    mu=unique(media.mean.epsilon)), .(locus, media)]
  
  #also get how many non-NA values each test is calculated on
  dat[,t.test.sample.size:=length(which(!is.na(epsilon))), 
    .(Tn.mutant.id, media)]
  
  dat[,t.test.sample.size.gene:=length(which(!is.na(epsilon))), 
    .(locus, media)]
  
  #also calculate a median per gene in condition
  dat[,epsilon.median.mutant.condition:=median(epsilon, na.rm=T), 
    .(Tn.mutant.id, media)]
  
  dat[,epsilon.median.gene.condition:=median(epsilon, na.rm=T), 
    .(locus, media)]
  
  dat[,epsilon.median.condition:=median(epsilon, na.rm=T), .(media)]
}


#' @export
getTTestResults = function (dat) {
  out = dat[,
    .(Tn.mutant.id, gene.name, gene.name.to.show, colony, locus, media,
      Tn.per.gene, copy.n, copy.id, t.test.statistic, t.test.pvalue,
      t.test.sample.size, t.test.statistic.gene, t.test.pvalue.gene,
      t.test.sample.size.gene, epsilon.median.mutant.condition,
      epsilon.median.gene.condition)]
  
   out[!duplicated(out)]
  
}


#' @export
correctForMultipleTesting = function (dat) {
  #do a multiple testing correction per media

  dat[,t.test.q.value:=p.adjust(t.test.pvalue, method="BH"), media]
  dat[,t.test.q.value.gene:=p.adjust(t.test.pvalue.gene, method="BH"), media]
  
  #do also a multiple testing correction using the new IHW method
  dat[,t.test.IHW.adj.p.value := IHW::adj_pvalues(IHW::ihw(t.test.pvalue,
      t.test.statistic, alpha = 0.05)), media]
  dat[,t.test.IHW.adj.p.value.gene :=
    IHW::adj_pvalues(IHW::ihw(t.test.pvalue.gene, t.test.statistic.gene, 
        alpha = 0.05)), media]
  
  unique(dat, by = c('Tn.mutant.id', 'media'))
  
}

# Reverse Genetics Validation --------------------

#' @export
addDrugRatio = function (dat) {
  lut = dat[step == 1 & cond == 3, .(rat = unique(c2/c1)), date]
  dat[lut, on = 'date']
}


#' @export
addAUC = function (dat) {
  dat[, AUC:=DescTools::AUC(time_h, OD620), .(date, well)]
}

