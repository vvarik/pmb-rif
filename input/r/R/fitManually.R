#' Fit marginal manually
#' 
#' More here https://github.com/openanalytics/BIGL/issues/4 
#' Somehow fitting monotherapies (marginal curves) with drc and feeding into
#' BIGL::fitMarginals just the problematic coefficient (Emax, BIGL's 'm') works
#'
#' @export
fitManually = function(grp, grp.val, dat = dat_rs) {

  foo = as.data.table(dat)[get(grp) %in% grp.val]
  
  bar = rbind(
    foo[d2 == 0, .(dose = d1, effect, ab = ab1)],
    foo[d1 == 0, .(dose = d2, effect, ab = ab2)]
  )
  
  ## Fit joint marginal model
  getParm = function(x) {
    fit = drc::drm(effect ~ dose, ab, data = x,  
      fct=drc::LL2.4(names=c("h", "m", "b", "e"), fixed=c(1, NA, NA, NA)),
      lowerl=c(min(x$effect), -Inf, -Inf),
      pmodels = data.frame(ab, 1, ab),
      robust = 'median'
    )

    # take Emax, BIGL calls it 'm'
    #coef(fit)[grepl('m', names(coef(fit)))] %>% 
      #setNames(., c('m1', 'm2'))
    setNames(coef(fit), c('m2', 'm1', 'b', 'e2', 'e1' ))

  }

  foo = as.data.frame(foo)
 
  marginal_fit = fitMarginals(as.data.frame(foo), 
    # play around with parms 1:3 (Emax1, Emax2, Emin(common))
    fixed = c('h1' = 1, 'h2' = 1, getParm(bar)[1:3])
  )

  return(marginal_fit)

}
