#' debias
#' 
#' Applies bias correction derived from forecast and observation data to forecast data set
#' 
#' @param fcst n x m x k array of n lead times, m forecasts, of k ensemble members
#' @param obs n x m matrix of veryfing observations
#' @param method character string with bias correction method name
#' @param crossval logical, should leave-one-out crossvalidation be used 
#' (see details)?
#' @param blocklength block length for moving blocks cross-validation 
#' (defaults to 5)
#' @param forward logical, should only past hindcasts be used for calibration?
#' @param fcst.out array of forecast values to which bias correction
#' should be applied (defaults to \code{fcst})
#' @param fc.time forecast dates of class 'Date' (for monthly correction, 
#' see \code{\link{monthly}})
#' @param fcout.time forecast dates of class 'Date' (for monthly correction, 
#' see \code{\link{monthly}})
#' @param ... additional arguments passed to bias correction methods
#' 
#' @details
#' If \code{crossval} is set to \code{TRUE}, the debiasing for years in block \code{i} are
#' computed based on the forecast and observation data set excluding years in block \code{i}. 
#' If, in addition, there are more years in the output set \code{fcst.out} than in the input 
#' set \code{fcst}, the bias correction for the remaining years in \code{fcst.out} is computed 
#' based on all years in \code{fcst}.
#' 
#' @examples
#' ## initialise forcast observation pairs
#' fcst <- array(rnorm(30*215*51, mean=1, sd=rep(seq(0.5,2, length=30), each=215)), 
#' c(215, 30, 51)) + 0.5*sin(seq(0,4,length=215))
#' obs <- array(rnorm(30*215, mean=2), c(215, 30)) + sin(seq(0,4, length=215))
#' fcst.debias <- debias(fcst, obs, 'unbias')
#' ## should be exactly zero
#' range(rowMeans(obs, dims=1) - rowMeans(fcst.debias, dims=1))
#' 
#' @keywords util
#' @export
debias <- function(fcst, obs, method='unbias', crossval=FALSE, blocklength=5, forward=FALSE, fcst.out=fcst, fc.time=NULL, fcout.time=fc.time, ...){
  ## get name of bias correction function
  dfun <- try(get(method), silent=TRUE)
  if (class(dfun) == 'try-error') stop('Bias correction method has not been implemented yet')
  if (crossval & forward) stop('Chose only one of forward or cross-validation')
  if (forward & !all(fcst == fcst.out)){
    stop('Forward only works with default fcst.out as of yet')
  }
  if (!all(dim(fcout.time) == dim(fcst.out)[1:2])){
    warning('Time for fcst.out not know -- inferred from time for fcst (fc.time)')
  }
  
  ## apply bias correction function
  if (crossval){
    fcst.debias <- array(NA, dim(fcst.out))
    ## figure out number of blocks
    for (i in seq(1, ncol(fcst), blocklength)){
      ii <- seq(i, min(i+blocklength - 1, ncol(fcst)))
      fcst.debias[,ii,] <- dfun(fcst=fcst[,-ii,,drop=FALSE], 
                               obs=obs[,-ii,drop=FALSE], 
                               fcst.out=fcst.out[,ii,,drop=FALSE], 
                               fc.time=if (is.null(fc.time)) NULL else fc.time[,-ii,drop=FALSE], 
                               fcout.time=if (is.null(fcout.time)) NULL else fcout.time[,ii,drop=FALSE],
                               ...)
    }
    ## compute the bias for the remaining years from full set
    ## if there are more years in the output than in the input
    if (ncol(fcst.out) > ncol(fcst)){
      ii <- seq(ncol(fcst)+1, ncol(fcst.out))
      fcst.debias[,ii,] <- dfun(fcst=fcst, 
                                obs=obs, 
                                fcst.out=fcst.out[,ii,,drop=F], 
                                fc.time=if (is.null(fc.time)) NULL else fc.time[,1:ncol(fcst),drop=FALSE], 
                                fcout.time=if (is.null(fcout.time)) NULL else fcout.time[,ii,drop=FALSE], 
                                ...)
    }
  } else if (forward) {
    fcst.debias <- array(fcst.out, dim(fcst.out))
    ## loop through n-years of forecasts to debias all the forecast in year i
    ## with all the forecasts before year i. The first forecast is not debiased at all.
    if (ncol(fcst.out) > 5){
      for (i in 5:ncol(fcst.out)){
        fcst.debias[,i,] <- dfun(fcst=fcst[,seq(1,i-1),, drop=FALSE],
                                 obs=obs[,seq(1,i-1),drop=FALSE],
                                 fcst.out=fcst.out[,i,,drop=FALSE],
                                 fc.time=if(is.null(fc.time)) NULL else fc.time[,1:i,drop=FALSE], 
                                 fcout.time=if (is.null(fcout.time)) NULL else fcout.time[,i,drop=FALSE], 
                                 ...)
      }      
    }
  } else {
    fcst.debias <- dfun(fcst=fcst, 
                        obs=obs, 
                        fcst.out=fcst.out, 
                        fc.time=fc.time, 
                        ...)
  }
  return(fcst.debias)
}