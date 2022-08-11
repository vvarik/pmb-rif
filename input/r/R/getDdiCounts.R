#' Get drug-drug interaction counts
#'
#' Function extracts the number of synergies and factor encoding to ready for
#' plotting.
#' @param rs a list of a list of response surface objects
#' @export

getDdiCounts = function(rs) {
  ret = lapply(seq_along(rs), function(i) {
    cond = names(rs[[i]])
    out = lapply(seq_along(rs[[i]]), function(j) {
      summary(rs[[i]][[j]][['rsl']])[['maxR']][['totals']]
    })
    names(out) = cond
    out
  })
  
  names(ret) = names(rs)

  lapply(ret, rbindlist, idcol = 'Cond') %>% 
    rbindlist(idcol = 'Strain') %>% 
    mutate(Strain = fct_relevel(Strain, 'ATCC27853', after = Inf)) %>% 
    mutate(Cond = factor(Cond, levels = c('intra', 'pH7.4', 'pH5.5'))) %>% 
    mutate(Name = recode(Strain, 
      A112      = 'A. baumannii A112',
      A113      = 'A. baumannii A113',
      E15       = 'E. cloacae E15',
      E51       = 'E. cloacae E51',
      K58       = 'K. pneumoniae K58',
      K74       = 'K. pneumonia K74',
      ATCC27853 = 'P. aeruginosa ATCC27853',
      PAO1      = 'P. aeruginosa PAO1',
      PA14      = 'P. aeruginosa PA14'
    ))
}
