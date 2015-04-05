#' Calibration of daily time series
#' 
#' Applies bias correction derived from forecast and observation data to
#' forecast data set
#' 
#' @param fcst n x m x k array of n lead times, m forecasts, of k ensemble
#'   members
#' @param obs n x m matrix of veryfing observations
#' @param method character string with bias correction method name
#' @param crossval logical, should leave-one-out crossvalidation be used (see
#'   details)?
#' @param blocklength block length for moving blocks crossvalidation (defaults
#'   to 1 for leave-one-out crossvalidation)
#' @param forward logical, should only past hindcasts be used for calibration?
#' @param fcst.out array of forecast values to which bias correction should be
#'   applied (defaults to \code{fcst})
#' @param fc.time forecast dates of class 'Date' (for monthly correction, see
#'   \code{\link{monthly}})
#' @param fcout.time forecast dates of class 'Date' (for monthly correction, see
#'   \code{\link{monthly}})
#' @param nforward number of forecasts to debias normally in forward mode (see
#'   details).
#' @param ... additional arguments passed to bias correction methods
#'   
#' @details If \code{crossval} is set to \code{TRUE}, the debiasing for years in
#' block \code{i} are computed based on the forecast and observation data set
#' excluding years in block \code{i}. If, in addition, there are more years in
#' the output set \code{fcst.out} than in the input set \code{fcst}, the bias
#' correction for the remaining years in \code{fcst.out} is computed based on
#' all years in \code{fcst}.
#' 
#' If \code{forward} is set to \code{TRUE}, the debiasing for forecast \code{i}
#' is computed based on all previous forecast observation pairs. The first
#' \code{nforward} forecasts, however, are debiased normally (as with
#' \code{forward = FALSE}).
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
debias <- function(fcst, obs, method='unbias', crossval=FALSE, blocklength=1, forward=FALSE, fcst.out=fcst, fc.time=NULL, fcout.time=fc.time, nforward=10, ...){
  ## get name of bias correction function
  dfun <- try(get(method), silent=TRUE)
  if (class(dfun) == 'try-error') stop('Bias correction method has not been implemented yet')
  if (crossval & forward) stop('Chose only one of forward or cross-validation')
  if (forward){
    ## check on forecasts    
    flen <- min(ncol(fcst), ncol(fcst.out))
    if (!all(fcst[,1:flen,] == fcst.out[,1:flen,])){
      stop('Forward only works with default fcst.out as of yet')      
    } else if (!all(dim(fcst) == dim(fcst.out))) {
      warning("fcst.out is assumed to start at same time as fcst")
    }
  }
  if (!is.null(fc.time) & !all(dim(fcout.time) == dim(fcst.out)[1:2])){
    warning('Time for fcst.out not know -- inferred from time for fcst (fc.time)')
  }

  ## missing values are not tolerated in forecast
  if (method != 'ccr'){
    stopifnot(!is.na(fcst), !is.na(fcst.out))
    stopifnot(!is.na(obs))
    if (!is.null(fc.time)) stopifnot(!is.na(fc.time), !is.na(fcout.time))
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
    fcst.debias <- fcst.out
    nforward <- min(nforward, ncol(fcst.out))
    nfcst <- ncol(fcst)
    
    ## new approach: first nforward years are calibrated backwards
    for (i in seq(1, min(nforward, ncol(fcst.out)))){
      fcst.debias[,i,] <- dfun(fcst=fcst[,seq(i+1,nfcst),, drop=FALSE],
                               obs=obs[,seq(i+1, nfcst),drop=FALSE],
                               fcst.out=fcst.out[,i,,drop=FALSE],
                               fc.time=if(is.null(fc.time)) NULL else fc.time[,seq(i+1, nfcst),drop=FALSE], 
                               fcout.time=if (is.null(fcout.time)) NULL else fcout.time[,i,drop=FALSE], 
                               ...)
    }      
    
#     ## loop through n-years of forecasts to debias all the forecast in year i
#     ## with all the forecasts before year i. 
#     ## The first nforward forecasts are debiased normally (in sample)
#     fcst.debias[,1:nforward,] <- dfun(fcst=fcst[,1:nforward,,drop=FALSE],
#                                       obs=obs[,1:nforward,drop=FALSE],
#                                       fcst.out=fcst.out[,1:nforward,,drop=FALSE],
#                                       fc.time=if(is.null(fc.time)) NULL else fc.time[,1:nforward,drop=FALSE], 
#                                       fcout.time=if (is.null(fcout.time)) NULL else fcout.time[,1:nforward,drop=FALSE], 
#                                       ...)
    if (ncol(fcst.out) > nforward){
      for (i in seq(nforward + 1, ncol(fcst.out))){
        fcst.debias[,i,] <- dfun(fcst=fcst[,seq(1,min(i-1, nfcst)),, drop=FALSE],
                                 obs=obs[,seq(1,min(i-1, nfcst)),drop=FALSE],
                                 fcst.out=fcst.out[,i,,drop=FALSE],
                                 fc.time=if(is.null(fc.time)) NULL else fc.time[,seq(1,min(i-1, nfcst)),drop=FALSE], 
                                 fcout.time=if (is.null(fcout.time)) NULL else fcout.time[,i,drop=FALSE], 
                                 ...)
      }      
    }
  } else {
    fcst.debias <- dfun(fcst=fcst, 
                        obs=obs, 
                        fcst.out=fcst.out, 
                        fc.time=fc.time, 
                        fcout.time=fcout.time,
                        ...)
  }
  return(fcst.debias)
}