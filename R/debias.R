#' @name debias
#' 
#' @title
#' Calibration of Daily Time Series
#' 
#' @description
#' Applies bias correction derived from forecast and observation data to
#' forecast data set
#' 
#' @param fcst n x m x k array of n lead times, m forecasts, of k ensemble
#'   members
#' @param obs n x m matrix of veryfing observations
#' @param method character string with bias correction method name
#' @param fcst.out array of forecast values to which bias correction should be
#'   applied (defaults to \code{fcst})
#' @param fc.time forecast dates of class 'Date' (for monthly correction, see
#'   \code{\link{monthly}})
#' @param fcout.time forecast dates of class 'Date' (for monthly correction, see
#'   \code{\link{monthly}})
#' @param strategy keyword for out-of-sample strategy (defaults to "none", i.e. 
#'   all values in obs/fcst are used to debias fcst.out), named list of parameters, or
#'   list of indices of obs/fcst pairs to be used for each instance in fcst.out 
#'   (see \code{\link[easyVerification]{indRef}}).
#' @param ... additional arguments passed to bias correction methods
#'   
#' @details 
#' No missing values are tolerated in either `obs` or `fcst` to ensure consistency
#' of calibration. Missing ensemble members, however, are tolerated in `fcst.out`,
#' thereby allowing calibration of non-homogeneous ensembles.
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
debias <- function(fcst, obs, method='unbias', fcst.out=fcst, 
                   fc.time=NULL, fcout.time=fc.time, strategy='none',
                   ...){
  ## check dimensions of fcst, obs, and fcst.out
  fodims <- dim(fcst.out)
  if (length(dim(fcst)) == 2 & length(obs) == nrow(fcst)){
    fcst <- array(fcst, c(1, dim(fcst)))
    obs <- array(obs, c(1, length(obs)))
    fcst.out <- array(fcst.out, c(1, dim(fcst.out)))
  }
  
  ## check for missing values
  nisna <- sum(!is.na(obs) & apply(!is.na(fcst), 1:2, any))
  if (nisna < 5) return(fcst.out*NA)

  if (!is.null(fc.time)) stopifnot(!is.na(fc.time), !is.na(fcout.time))
  
  
  ## get name of bias correction function
  dfun <- try(get(method), silent=TRUE)
  if (class(dfun) == 'try-error') stop('Bias correction method has not been implemented yet')
  if (!is.null(fc.time) & !all(dim(fcout.time) == dim(fcst.out)[1:2])){
    warning('Time for fcst.out not know -- inferred from time for fcst (fc.time)')
  }
  
  
  ## missing values are not tolerated in output for useqmap
  if (method == 'useqmap' & any(is.na(fcst.out))) stop("Missing values not tolerated in useqmap")
  if (method == 'useqmap' & any(is.na(obs))) stop("Missing values not tolerated in useqmap")
  
  ## deparse out-of-sample options for calibration
  if (length(strategy) == 1 & is.character(strategy)){
    cali.ind <- easyVerification::indRef(ncol(fcst.out), 
                                         type=strategy, 
                                         indices=1:ncol(fcst))
  } else if (is.list(strategy)){
    if (!is.null(names(strategy))){
      if (is.null(strategy$nfcst)) strategy$nfcst <- ncol(fcst.out)
      if (is.null(strategy$indices)) strategy$indices <- 1:ncol(fcst)
      cali.ind <- do.call(easyVerification::indRef, strategy)  
    } else if (is.numeric(unlist(strategy)) & all(sapply(strategy, function(x) x%%1) == 0)){
      cali.ind <- strategy
    } else {
      stop("wrong format of calibration strategy")
    }
  } else {
    stop("wrong format of calibration strategy")
  }
  
  
  ## apply bias correction function (in-sample or out-of-sample)
  fcst.debias <- array(NA, dim(fcst.out))
  cali.ind <- lapply(cali.ind, sort)
  unicali <- unique(cali.ind)
  for (cali in seq(along=unicali)){
    ii.cal <- unicali[[cali]]
    ii.out <- which(sapply(cali.ind, identical, ii.cal))
    fcst.debias[,ii.out,] <- dfun(fcst=fcst[,ii.cal,,drop=FALSE], 
                                  obs=obs[,ii.cal,drop=FALSE], 
                                  fcst.out=fcst.out[,ii.out,,drop=FALSE], 
                                  fc.time=if (is.null(fc.time)) NULL else fc.time[,ii.cal,drop=FALSE], 
                                  fcout.time=if (is.null(fcout.time)) NULL else fcout.time[,ii.out,drop=FALSE],
                                  ...)
  }

  return(array(fcst.debias, fodims))
}