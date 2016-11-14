#' @name debiasApply
#'   
#' @title Apply Bias Correction to Large Ensemble Forecast Data Sets
#'   
#' @description This wrapper allows to apply the bias correction methods 
#'   supplied with \code{\link{debias}} to large arrays of ensemble forecast and
#'   observation data.
#'   
#' @param fcst array of ensemble forecast data with arbitrary dimensions. The 
#'   final dimensions contain lead times, forecast instances (usually years) and
#'   ensemble members in this order.
#' @param obs array of observation used for calibration. Dimensions correspond
#'   to the \code{fcst} array without the last dimension.
#' @param method character string with bias correction method name
#' @param fcst.out array of ensemble forecast data that is to be calibrated. All
#'   but the forecast instance and ensemble dimensions have to be the same as
#'   for \code{fcst}.
#' @param ... Additional parameters passed on to \code{\link{debias}}
#' @param atomic logical, should debiasing be applied to each forecast instance
#'   separately (in contrast to jointly over lead times)?
#' 
#' @keywords util
#' @export
debiasApply <- function(fcst, obs, method='unbias', fcst.out=fcst, ..., atomic=FALSE) {

  ## set number of minimum dimensions
  nn <- ifelse(atomic, 2, 3)
  
  # enquire dimensions of fcst, obs, and fcst.out
  fdims <- dim(fcst)
  ndims <- length(fdims)
  odims <- dim(obs)
  foutdims <- fodims <- dim(fcst.out)
  # add ensemble member of one if none present
  if (ndims == length(odims)){
    fcst <- array(fcst, c(fdims, 1))
    fcst.out <- array(fcst.out, c(fodims, 1))
  }
  fdims <- dim(fcst)
  fodims <- dim(fcst.out)
  ndims <- length(fdims)
  stopifnot(ndims >= nn)
  stopifnot(odims == fdims[seq(along=odims)])
  stopifnot(fodims[seq(1, ndims - nn + 1)] == fdims[seq(1, ndims - nn + 1)])
  
  ## the trivial case
  if (ndims == nn){
    return(debias(fcst=fcst, obs=obs, method=method, fcst.out=fcst.out, ...))
  } 
  
  ## deal with the non-trivial case of multiple instances
  resdims <- fdims[seq(1,ndims - nn)]
  indims <- c(prod(resdims), fdims[ndims - seq(nn - 1, 0)])
  outdims <- c(prod(resdims), fodims[ndims - seq(nn - 1, 0)])
  
  ftmp <- array(fcst, indims)
  otmp <- array(obs, indims[-length(indims)])
  fotmp <- array(fcst.out, outdims)

  ## dimensions for arrays in debias
  ifd <- fdims[ndims - seq(nn - 1, 0)]
  iod <- odims[ndims - seq(nn - 1, 1)]
  ofd <- fodims[ndims - seq(nn - 1, 0)]

  if (nn == 3){
    fout <- sapply(1:nrow(ftmp), function(i) debias(fcst=array(ftmp[i,,,], ifd), 
                                                    obs=array(otmp[i,,], iod), 
                                                    method=method,
                                                    fcst.out=array(fotmp[i,,,], ofd), ...), 
                   simplify='array')
    
    fcst.debias <- array(aperm(fout, c(4,1,2,3)), fodims)
  } else if (nn == 2) {
    fout <- sapply(1:nrow(ftmp), function(i) debias(fcst=array(ftmp[i,,], ifd), 
                                                    obs=array(otmp[i,], iod), 
                                                    method=method,
                                                    fcst.out=array(fotmp[i,,], ofd), ...), 
                   simplify='array')
    
    fcst.debias <- array(aperm(fout, c(3,1,2)), fodims)    
  }
 
  return(array(fcst.debias, foutdims))
}