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
#' @param fcst.out array of ensemble forecast data that is to be calibrated. All
#'   but the forecast instance and ensemble dimensions have to be the same as
#'   for \code{fcst}.
#' @param ... Additional parameters passed on to \code{\link{debias}}
#' 
#' @keywords util
#' @export
debiasApply <- function(fcst, obs, fcst.out=fcst, ...) {
  # enquire dimensions of fcst, obs, and fcst.out
  fdims <- dim(fcst)
  ndims <- length(fdims)
  odims <- dim(obs)
  fodims <- dim(fcst.out)
  stopifnot(odims == fdims[-ndims])
  stopifnot(fodims[seq(1, ndims - 2)] == fdims[seq(1, ndims - 2)])
  stopifnot(ndims >= 3)
  
  ## the trivial case
  if (ndims == 3){
    return(debias(fcst=fcst, obs=obs, fcst.out=fcst.out, ...))
  } 
  
  ## deal with the non-trivial case of multiple instances
  resdims <- fdims[seq(1,ndims - 3)]
  indims <- c(prod(resdims), fdims[ndims - 2:0])
  outdims <- c(prod(resdims), fodims[ndims - 2:0])
  
  ftmp <- array(fcst, indims)
  otmp <- array(obs, indims[-length(indims)])
  fotmp <- array(fcst.out, outdims)

  ## dimensions for arrays in debias
  ifd <- fdims[ndims - 2:0]
  iod <- odims[ndims - 2:1]
  ofd <- fodims[ndims - 2:0]

  fout <- sapply(1:nrow(ftmp), function(i) debias(fcst=array(ftmp[i,,,], ifd), 
                                                  obs=array(otmp[i,,], iod), 
                                                  fcst.out=array(fotmp[i,,,], ofd), ...), 
                 simplify='array')
  
  fcst.debias <- array(aperm(fout, c(4,1,2,3)), fodims)
  
  return(fcst.debias)
}