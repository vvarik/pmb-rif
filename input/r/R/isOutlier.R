#' Find outliers
#'
#' Find outliers based on MAD, IQR, Dixon test, or Qn.
#' @param x A vector
#' @param method A character, one of 'mad'(default), 'iqr', 'dixon', or 'qn'.
#' @param alpha A numeric, significance level for Dixon test, defaults to 0.05.
#' @param pvalue A logical, whether to return pvalue for Dixon test.
#' @param thr A nonnegative constant that determines the extension of bounds 
#'   (alpha level). Commonly used values are 2, 2.5, 3 (default), and 3.5.
#' @param lod A numeric (defaults to NULL), lower boundary for signal that can
#' reliably distinguished from noise i.e. Limit of Detection. Make sure that if
#' x is log transformed, you do the same to lod.
#' @return Returns a logical vector by default. If test is Dixon and pvalue is
#'   set TRUE, returns a list consisting of logical vector and p-value.
#' @examples
#' x = c(0.5, 1, 2, 10)
#' isOutlier(x)
#' isOutlier(x, 'dixon')
#' isOutlier(x, 'dixon', 0.1)
#' isOutlier(x, 'dixon', pvalue = TRUE)
#' 
#' x = c(0.5, 1, 2, 9)
#' isOutlier(x)            # one outlier
#' isOutlier(x, 'dixon')   # no outliers
#' 
#' # cautionary example
#' set.seed(12)  
#' x = rnorm(4)  # normal distribution
#' isOutlier(x)                # yet one outlier
#' isOutlier(x, 'dixon')       # no outliers
#' isOutlier(x, 'dixon', 0.1)  # one outlier
#' isOutlier(x, 'dixon', pvalue = TRUE)
#' @references Boris Iglewicz and David Hoaglin (1993), Volume 16: How to
#' Detect and Handle Outliers, The ASQC Basic References in Quality Control:
#' Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
#' 
#' Rousseeuw PJ, Croux C. Alternatives to the median absolute deviation. J Amn
#' Statist Assoc 1993, 88:1273– 1283.
#' @export
isOutlier = function(x, method = 'mad', alpha = 0.05, pvalue = FALSE, 
  thr = 3, lod = NULL) {
  if (!method %in% c('iqr', 'dixon', 'mad', 'qn')) stop('Unknown method')
 
  if (!is.null(lod))
      if(diff(range(x)) < abs(lod)) return(rep(FALSE, length(x)))
  
  if (any(duplicated(x))) set.seed(9); x = jitter(x, amount = 0.0001)

  upper = quantile(x, 0.75, na.rm=T) + 1.5 * IQR(x)
  lower = quantile(x, 0.25, na.rm=T) - 1.5 * IQR(x)

  if (method == 'iqr') out = x > upper | x < lower

  else if(method == 'dixon'){
      pval = tryCatch(outliers::dixon.test(x)$p.value, 
                   error=function(cond) return(NA_real_)
      )
      if(is.na(pval)) outl = rep(NA, length(x))
      else if (pval < alpha) outl = outliers::outlier(x, logical = TRUE)
      else outl = rep(FALSE, length(x))
      
      if(pvalue) out = list(is_outl = outl, pval = rep(pval, length(x)))
      else out = outl 

  } else if(method == 'mad') {
      med = median(x)
      dev = abs(x - med)
      mad = median(dev)
      modified_z_score = 0.6745 * dev / mad
      out = modified_z_score > thr

  } else if(method == 'qn') {
      med = median(x)
      dev = abs(x - med)
      qn  = robustbase::Qn(x)
      modified_z_score = dev / qn
      out =  modified_z_score > thr
  }

  return(out)
}
