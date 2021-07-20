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
getPA14Map = function () {
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
}


#' @export
addAnnotation = function (dat) {
  PA14.map = getPA14Map()
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
  dat[, AUC:=cumAUC(time_h, OD620), .(date, well)]
}


#' @export
addFitness = function (dat) {
  # full growth of strain, OD_bug
  dat[, AUC_bug  := gmu(AUC[step == 8]), .(date, mut, time_h)]
  # fitness, fit_AUC
  lut = dat[, .(fit_AUC=gmu(AUC/AUC_bug), AUC_bug), 
    .(date, mut, cond, step, time_h)]
  lut[, fit_AUC := ifelse(AUC_bug < 1, NA, fit_AUC)]
  
  merge(dat, lut, all.x=T)
}


#' @export
remDissimilarExperimentalConditions = function (dat) {
  dat[!(mut!=28 & date == '2019-12-01')] %>%  # o/n cells there, different inhibition
  .[!(rat == 8 & cond_f %in% c('combo', 'RIF'))] %>% 
  # for slowly growing mutant, I grew the culture o/n prior to performing the
  # experiment on the following dates
  .[!(mut == 28 & !date %in% as.Date(c('2019-12-01', '2019-12-08', '2019-12-20', '2019-12-21')))]
}


#' @export
fitPA14MutDR = function(mydat, x=dose, y=fit_AUC,
                 curveid=grp, ...) {
  # Fit dose-responses of fitness
  y = deparse(substitute(y))
  x = deparse(substitute(x))
  curveid = substitute(curveid)
  eq = reformulate(x, response=y)

  try(
    drm(eq, curveid=eval(curveid), data=mydat,
    fct=LL2.4(names=c("hill", "Emax", "Emin", "logEC50"),
                    fixed=c(NA, NA, NA, NA)),
    pmodels = data.frame(eval(curveid), 1, 1, eval(curveid)),
    lowerl=c(1, 0, -Inf, -Inf),
    upperl=c(10, Inf, 1.2, Inf),
    ...
    ) 
  )
}


#' @export
pltPA14MutDR = function(mydat, cex=1.5, ylim=c(0,1.1), lty=1, 
                 showname = TRUE, ...
) {
  plottype = substitute(plottype)
  mut_name = unique(mydat$origData$mut)
  label_x = 'PMB, RIF/48, µg/mL'
  try({
    plot(mydat, cex.axis=1.5, lty=lty, ylim=ylim, cex=cex,
         xlab=label_x, ylab='Fitness', broken=TRUE, ...
    ) 
    if(showname) mtext(mut_name)
  })
}

#' @export
addPA14MutAnnotation = function (dat) {
  candidates = fread('input/dat/raw/pa14_mut_to_validate.csv') %>% 
    .[, .(locus, mut=no)]
  pa14_map = getPA14Map() %>% 
    .[, gene := ifelse(gene.name=='', locus, gene.name)]
  setkey(pa14_map, locus)
  lut = pa14_map[candidates, on = 'locus'][, .(mut, gene)] %>% 
    unique(., by = 'mut')
  lut[is.na(gene), gene := 'wild-type']
  setkey(lut, mut)
  dat[lut]
}


#' @export
plt45PA14MutantDR = function (organize=T) {

  # order
  if(organize == T) {
    vec = c(
      # no synergy
      'relA', 'PA14_02150', 'PA14_26590',           

      # less than wild-type
      'lldP', 'PA14_56840', 'PA14_60490', 'PA14_66480', 

      # wild-type
      'wild-type',                                 

      # more than wild-type 
      'maiA', 'pmbA', 'PA14_03760', 'PA14_43270',   
      'PA14_43670', 'PA14_51310', 'PA14_62230',   

      # same as wild-type
      'ahpF', 'fdnH', 'mexB', 'mmsR', 'mucB', 'pilF', 'pmrB', 'purE',
      'PA14_11130', 'PA14_12450', 'PA14_15540', 'PA14_17990', 'PA14_18680',
      'PA14_19150', 'PA14_19670', 'PA14_21710', 'PA14_22120', 'PA14_26960',
      'PA14_34820', 'PA14_36150', 'PA14_37120', 'PA14_40730', 'PA14_42900',
      'PA14_47260', 'PA14_49880', 'PA14_53380', 'PA14_59630', 'PA14_65520',
      'PA14_65670', 'PA14_68670'
    )
    foo = sapply(fit, function(x) unique(x$orig$gene))
    idx = c(101:103, match(vec, foo))
  } else {
    idx = seq_along(fit)
  }

  setPar(mfrow = c(8, 6))
  for (i in idx) {

    if(i > 100) {
      plot.new()
      next
    }

    lapply(fit[1], pltPA14MutDR, xlim=c(0.3, 10), lwd=6, 
      pch = 19, type = 'none', main = unique(fit[[i]]$ori$gene),
      col=scales::alpha(drugCols[c('PmbWt','RifWt','CmbWt')], 0.75), 
      showname = FALSE, 
      legend = FALSE, axes = FALSE
    )
    axis(side = 1, at = c(0.3, 1, 3, 10))
    my_ylab = c(0.0,  0.5, 1.0)
    axis(2, las = 1, at = my_ylab, labels = my_ylab)
    # axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0))

    if(i == 1) next

    lapply(fit[i], pltPA14MutDR, showname = FALSE,
      xlim=c(0.3, 10), col=drugCols[c('PmbMono', 'RifMono', 'Cmb')], 
      cex=1.5, pch = 19, add = T, legend = F # legendPos = c(12, 1.15)
    )
  }
}


#' @export
plt45PA14MutantComparisonToLoeweNull = function(fname, organize=T) {
  fname = paste0('input/dat/rds/', fname, '.rds') 
  rsl = readRDS(fname) %>% lapply(., function(x) x$rsl)

  # mut 28 failed surface analysis
  fitNew = fit[-29]

  # order
  if(organize == T) {
    vec = c(
      # no synergy
      'relA', 'PA14_02150', 'PA14_26590',           

      # less than wild-type
      'lldP', 'PA14_56840', 'PA14_60490', 'PA14_66480', 

      # wild-type
      'wild-type',                                 

      # more than wild-type 
      'maiA', 'pmbA', 'PA14_03760', 'PA14_43270',   
      'PA14_43670', 'PA14_51310', 'PA14_62230',   

      # same as wild-type
      'ahpF', 'fdnH', 'mexB', 'mmsR', 'mucB', 'pilF', 'pmrB', 'purE',
      'PA14_11130', 'PA14_12450', 'PA14_15540', 'PA14_17990', 'PA14_18680',
      'PA14_19150', 'PA14_19670', 'PA14_21710', 'PA14_22120', 'PA14_26960',
      'PA14_34820', 'PA14_36150', 'PA14_37120', 'PA14_40730', 'PA14_42900',
      'PA14_47260', 'PA14_49880', 'PA14_53380', 'PA14_59630', 'PA14_65520',
      'PA14_65670', 'PA14_68670'
    )
    foo = sapply(fitNew, function(x) unique(x$orig$gene))
    idx = c(101:103, match(vec, foo))
  } else {
    idx = seq_along(fitNew)
  }

  setPar(mfrow=c(8,6))
  for (i in idx) {

    if(i > 100) {
      plot.new()
      next
    }

    foo = rsl[[1]]$offAxisTable  # wild-type
    bar = list(drm(predicted ~ d1, data = foo, fct = LL2.4()))
    lapply(bar, pltPA14MutDR, xlim=c(0.3, 3), 
      main = unique(fitNew[[i]]$ori$gene), showname=F,
      type = 'none', lwd = 6, col=alpha(drugCols[c('LoeweNullWt')], 0.75), axes
      = FALSE, cex=1.5, pch = 19, legendPos = c(1, 0.5))
    plot(fitNew[[1]], level = 'combo', add = TRUE, col = drugCols[c('CmbWt')], 
      cex = 1.5, pch=19)
    
    axis(side = 1, at = c(0.3, 1, 3, 10))
    my_ylab = c(0.0,  0.5, 1.0)
    axis(2, las = 1, at = my_ylab, labels = my_ylab)
    # axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0))

    if(i == 1) next
  
    foo = rsl[[i]]$offAxisTable
    bar = list(drm(predicted ~ d1, data = foo, fct = LL2.4()))
    lapply(bar, pltPA14MutDR, xlim=c(0.3, 3), type = 'none', lwd = 6,
           col=alpha(drugCols[c('LoeweNullMut')], 0.5), add = TRUE,
           cex=1.5, pch = 19)
    
    lapply(fitNew[i], pltPA14MutDR, showname = F, 
      level = 'combo', xlim=c(0.3, 10),
      col=drugCols[c('Cmb')], cex=1.5, pch = 19, add = T
    )
  
  }
}


#' @export
getPa14RS = function (dat) {
  dat_rs = dat[time_h == 14] %>% 
    .[, mut := mut + 1] %>% 
    dplyr::select(mut, d1 = c1, d2 = c2, ab1 = abb1, ab2 = abb2, effect =
      fit_AUC, cond) %>% 
    as.data.frame()  # some point down the road handles only df's without complaint
  
  # there is a problem with mut24
  #dat_rs = dat_rs[dat_rs$mut != 23, ]

  cores_n = parallel::detectCores()
  doParallel::registerDoParallel(cores_n)
  `%dopar%` = foreach::`%dopar%`

  set.seed(9)
  rs = foreach::foreach (i = unique(dat_rs$mut)) %dopar% {
  # for (i in 1:3) {
    cat('Working on', i, '\n')
    ## Select experiment
    data = subset(dat_rs, mut==i)
    data = data[complete.cases(data), ]
  
    ## Fit joint marginal model
    marginal_fit = BIGL::fitMarginals(data, method='nlslm')
    ## Predict response surface based on generalized Loewe model
    rs_loewe = try(BIGL::fitSurface(data, marginal_fit, null_model='loewe2',
                      statistic = "both", B.CP = 500,
                      parallel = FALSE))
  
    rs_bliss = try(BIGL::fitSurface(data, marginal_fit, null_model='bliss',
                      statistic = "both", B.CP = 500,
                      parallel = FALSE))
    # MeanR, tests the overall fit of the data to zero interaction model
    # MaxR identifies concentrations where synergy or antagonism is present
    maxR_loewe = tryCatch(summary(rs_loewe[['maxR']])[['totals']], error=function(err) NA)
    maxR_bliss = tryCatch(summary(rs_bliss[['maxR']])[['totals']],
    error=function(err) NA)
    
    list(
      mth = marginal_fit,  # monotherapy
      rsl = rs_loewe,
      rsb = rs_bliss,
      maxRl = maxR_loewe,
      maxRb = maxR_bliss
    )
  
  }
  doParallel::stopImplicitCluster()

  return(rs)
}


#' @export
addNewVariablesToTimeKill = function (dat) {
  dat[, new_time_h := ifelse(time_h == 24, 10, time_h)] 
  dat[, group := factor(group, labels=c("no drug", "low RIF", "low PMB",
      "low RIF + low PMB", "10x higher RIF", "10x higher PMB", "low PMB", "10x higher PMB"))]
  dat[, outl := isOutlier(effect), .(broth, group, time_h)]
}


#' @export
pltTimeKill = function (dat) {
  myColors = rep(drugCols[c(7, 1, 5, 3, 1, 5)], 3)
  names(myColors) = NULL
  dat[!(outl) & !(group == 'low RIF' & new_time_h %in% c(0.5, 0.75))] %>% 
    group_by(broth, group, new_time_h) %>% 
    summarise(response=mean(effect, na.rm=T), ci=se(effect)) %>% 
    ggplot(aes(new_time_h, response, group=interaction(broth, group), shape=broth, color=group)) + 
    geom_line(na.rm=T, size = 0.75) +
    geom_pointrange(aes(ymin = response - ci, ymax = response + ci), size=0.75, fatten=5, fill='white') + 
    scale_shape_manual(name = '', values=c(19, 21, 22)) +
    facet_wrap(vars(group), nrow=1) +
    scale_color_manual(name="", values=myColors, guide=F) +
    scale_x_continuous(breaks=c(0, 2, 4, 6, 10), labels=c(0, 2, 4, 6, 24)) +
    scale_y_continuous(breaks=-4:4) +
    labs(x = 'Time, h', y = expression(log[10]*'('*over('CFU'['24h'], 'CFU'['0h '])*')')) +
    theme(legend.position='top', strip.text.x = element_text(size = 16))
}


#' @export
getTopGO = function (dat) {

  out = vector('list', length=2)

  # GO terms need to be in certain format for topGO
  pa14_go = fread('input/dat/raw/ucbpp-pa14_go.csv')
  gene_go = pa14_go[ , c('Locus Tag', 'Accession')] 
  
  # Create list with element for each gene, containing vectors with all terms for each gene
  gene2GO = tapply(gene_go$'Accession', gene_go$'Locus Tag', function(x) x)
  
  # Focus on the genes 
  dat = unique(dat, by = c('media', 'gene.name.to.show'))
  
  fdr_cutoff = 0.05
  tmp = dat[media=='LB' & !is.na(t.test.q.value)]  # KS test below fails with NAs
  gene_list = ifelse(tmp$t.test.q.value<fdr_cutoff, 1, 0)
  names(gene_list) = tmp$locus
  
  # Create topGOData object. Ontology can be one of: BP - biological process, CC
  # - cellular component, BP - molecular function
  GOdata = new("topGOdata",
      ontology = "CC",
      allGenes = gene_list,
      geneSelectionFun = function(x) (x==1),
      annot = annFUN.gene2GO, 
      gene2GO = gene2GO
  )
  
  # Kolmogorov-Smirnov testing. topGO preferentially tests more specific terms,
  # utilizing the topology of the GO graph. Details on algorithms:
  # https://academic.oup.com/bioinformatics/article/22/13/1600/193669
  resultKS = runTest(GOdata, algorithm = 'weight01', statistic = "ks")

  out[['CC']] = GenTable(GOdata, raw.p.value = resultKS, 
    topNodes = length(resultKS@score), numChar = 120
  )
  

  GOdata = new("topGOdata",
      ontology = "BP",
      allGenes = gene_list,
      geneSelectionFun = function(x) (x==1),
      annot = annFUN.gene2GO, 
      gene2GO = gene2GO)
  
  # Kolmogorov-Smirnov testing. topGO preferentially tests more specific terms,
  # utilizing the topology of the GO graph. Details on algorithms:
  # https://academic.oup.com/bioinformatics/article/22/13/1600/193669
  resultKS = runTest(GOdata, algorithm = 'weight01', statistic = "ks")

  out[['BP']] = GenTable(GOdata, raw.p.value = resultKS, 
    topNodes = length(resultKS@score), numChar = 120
  )
  
  out
  
}


#' @export
getValidHitTable = function () {
  tmpFun = function (x) unique(dat[mut %in% x], by = 'gene')$gene
  valid_hits = data.table(
    Gene = c(
      tmpFun(c(23, 31, 33)),  # no synergy
      tmpFun(c(1, 16, 29, 35)),  # less synergy than wt
      tmpFun(c(4, 5, 14, 19, 40, 45, 12))  # more synergy
    ),
    Synergy = c(
      rep('none', 3),
      rep('less', 4),
      rep('more', 7)
      )
  )
  
  valid_hits = getPA14Map() %>% 
    .[, .(Locus=locus, Gene=gene.name.to.show)] %>% 
    unique() %>% 
    .[valid_hits, on = 'Gene'] %>% 
    .[, Locus := ifelse(is.na(Locus), Gene, Locus)]
  
  foo = fread('input/dat/raw/PA14_computational_v2_updated_annotations.csv')
  setnames(foo, c('\"\"Active\"\" Gene Locus', '\"\"Active\"\" Gene Description'),
      c('Locus', 'Description'))
  
  out = foo[, .SD, .SDcols=grepl('Desc|V29|Locus$', names(foo))
    ][valid_hits, on = 'Locus'] %>% 
    unique() %>% 
    .[, Domain := gsub('.*\\|(.*) \\[.*]', '\\1', V29)] %>% 
    dplyr::select(-'V29')  %>% 
    mutate(Description =
      ifelse(grepl('hypothetical|pyrophosphokinase', Description), Domain,
        Description),
      Description = gsub('PpGpp', '(p)ppGpp', Description),
      Description = gsub(' catalytic domain', '', Description),
      Gene = gsub('PA14_', '', Gene)
    )  %>% 
    dplyr::select(-Domain)

  return(out)
  
}


#' @export
addEcOrthologs = function (dat) {
  ort = fread('input/dat/raw/pa14_ecMG1655_orthologs.csv')
  bar = str_detect(names(ort), 'Locus Tag ')
  ort = ort[, ..bar]
  names(ort) = c('locus_ec', 'locus_pa')
  
  ec_map = fread('input/dat/raw/ec_genome_ref_2020-01-24.txt') %>% 
    setnames(., 'Locus tag', 'locus_ec')
  
  out = ort[dat, on = c(locus_pa = 'Locus')]
  
  out = ec_map[, .(locus_ec, Locus)][out, on = 'locus_ec']
  out = dplyr::rename(out, 'EC gene' = Locus, 'PA gene'=Gene) %>% 
    .[, -'locus_ec'] %>% 
    dplyr::mutate('EC gene' = ifelse(is.na(`EC gene`), '', `EC gene`)) %>% 
    dplyr::select('EC gene', 'PA gene', Synergy, Description, locus_pa) 
  
}


#' @export
getGoTermTable = function () {
  # go terms quick and dirty; output -- I will concatenate/edit them manually and
  # reload back in below
  # pa14_go = fread('input/dat/raw/ucbpp-pa14_go.csv') %>% 
  #   dplyr::rename(locus_pa = `Locus Tag`)
  # fwrite(pa14_go[, .(`GO Term`, Namespace, locus_pa)][tab, on = 'locus_pa'][Namespace == 'biological_process'], 'input/dat/raw/tab_go_terms_biol_processes.csv')
  # fwrite(pa14_go[, .(`GO Term`, Namespace, locus_pa)][tab, on = 'locus_pa'][Namespace == 'cellular_component'], 'input/dat/raw/tab_go_terms_cellular_component.csv')
  # fwrite(pa14_go[, .(`GO Term`, Namespace, locus_pa)][tab, on = 'locus_pa'][Namespace == 'molecular_function'], 'input/dat/raw/tab_go_terms_molec_func.csv')
  
  foo = fread('input/dat/raw/hits_go_terms_molec_func.csv') %>% setkey(., 'PA_gene')
  bar = fread('input/dat/raw/hits_go_terms_biol_processes.csv') %>% setkey(., 'PA_gene')
  baz = fread('input/dat/raw/hits_go_terms_cellular_component.csv') %>% setkey(., 'PA_gene')
  
  merge(merge(foo, bar, all = TRUE), baz, all = TRUE) %>% 
    dplyr::select('PA gene' = PA_gene, Location, Process, Function)
}
