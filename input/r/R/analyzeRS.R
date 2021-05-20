#' Perform response surface analysis
#' 
#' Old name of the function, for reference purposes: responseSurface 
#' @export
analyzeRS = function(dat = dat_rs, grp, B.CP=500, x1='d1', x2='d2') {

  cores_n = parallel::detectCores()
  doParallel::registerDoParallel(cores_n)
  `%dopar%` = foreach::`%dopar%`
  dat[, d1 := dat[, get(x1)]]
  dat[, d2 := dat[, get(x2)]]

  set.seed(9)
  out = foreach::foreach (i = unique(dat[, get(grp)])) %dopar% {
    cat('=== Working on', i, ' ===\n')
    ## Select experiment
    data = subset(dat, get(grp)==i)
    data = data[complete.cases(data[, .(effect, d1, d2)]), ]

    ## Fit joint marginal model
    marginal_fit = fitMarginals(data, fixed = c('h1' = 1, 'h2' = 1))
    ## Predict response surface based on generalized Loewe model
    rs_loewe = fitSurf(data, marginal_fit, B=B.CP)

    if(length(rs_loewe) == 1){
      marginal_fit = fitManually(grp, i, data)
      rs_loewe = fitSurf(data, marginal_fit, B=B.CP)
    }
  
    rs_bliss = fitSurf(data, marginal_fit, 'bliss', B.CP)
    # MeanR, tests the overall fit of the data to zero interaction model
    # MaxR identifies concentrations where synergy or antagonism is present
    maxR_loewe = tryCatch(
      summary(rs_loewe[['maxR']])[['totals']], error=function(err) NA
    )
    maxR_bliss = tryCatch(
      summary(rs_bliss[['maxR']])[['totals']], error=function(err) NA
    )
    
    list(
      mth = marginal_fit,  # monotherapy
      rsl = rs_loewe,
      rsb = rs_bliss,
      maxRl = maxR_loewe,
      maxRb = maxR_bliss
    )
  
  }
  doParallel::stopImplicitCluster()

  out

}
