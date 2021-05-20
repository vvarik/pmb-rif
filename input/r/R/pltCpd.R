#' Plot compound
#'
#' Plot individual trace against the backdrop of all others. Highlight DMSO
#' control
#' @export
pltCpd = function(cpd, spc=NULL, DT=dat, legPos='none', nRow=3) {

  if(is.null(spc)) {
    spc_id = unique(DT$Species)
  } else if(length(spc)>1){
    spc_id = spc 
  } else {
    spc_id = DT[grepl(spc, Species, ignore.case=T), unique(Species)]
  }


  plt_id = DT[Species %in% spc_id, unique(PltID)]
  exc_id = DT[grepl(cpd, Excipient, ignore.case=T), unique(Excipient)]

  clean_idx = mth[PltID %in% plt_id & Species %in% spc_id, .(PltID, Well)]
  dat = DT[clean_idx, on = c('PltID', 'Well')]
  dat = dat[(Keep)]

  dat[, grp := 'All']
  dat[Excipient == 'DMSO', grp := 'DMSO']
  dat[Excipient == exc_id, grp := Excipient]
  dat[, grp := factor(grp, levels = c('All', 'DMSO', exc_id))]


  p1 = ggplot(dat, aes(TimeH, OD, group = interaction(PltID, Well), col=grp)) +
    scale_y_continuous(trans = 'log2', breaks=myBreaks) +
    labs(y = expression('OD'[578]), x = 'Time, h') +
    facet_wrap(~Species, scales='free_x', nrow=nRow) +
    theme(legend.position=legPos,
      strip.text = element_text(face = "italic", size = .facet_size))
    
  
  colo = c("All" = cbPalette[1], "DMSO" = "grey30", cbPalette[2])
  names(colo)[3] = exc_id

  p1 + geom_line() +
   geom_line(data = . %>% .[Excipient == 'DMSO']) +
   geom_line(data = . %>% .[Excipient == exc_id]) +
   scale_color_manual(name='', values=colo)

}
