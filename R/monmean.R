#' Compute the monthly mean
#' 
#' @param x input array (2-d or 3-d, i.e. obs or forecasting)
#' @param time dates of input array (2-d, e.g. lead time x years)
#' 
#' @keywords util
#' @export
monmean <- function(x, time){
  datstr <- format(time, '%Y %m')
  nmon <- length(unique(gsub('.* ', '', datstr)))
  nyears <- ncol(x)
  if (length(dim(x)) == 2){
    xout <- tapply(x, datstr, mean, na.rm=T)
    if (length(xout) != nmon*nyears) stop('Number of months per year not stable')
    xout <- array(xout[unique(datstr)], c(nmon, nyears))
  } else {
    xout <- apply(x, 3, tapply, datstr, mean, na.rm=T)
    if (nrow(xout) != nmon*nyears) stop('Number of months per year not stable')
    xout <- array(xout[unique(datstr), ], c(nmon, nyears, dim(x)[3]))
  }
  rownames(xout) <- unique(gsub('.* ', '', datstr))
  return(xout)
}