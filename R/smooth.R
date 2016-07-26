#' @name smooth
#' @aliases smooth
#' @aliases smoothobs
#' @aliases unbias
#' 
#' @title
#' Mean De-biasing With Smoothing of Daily Climatology
#' 
#' @description
#' Computes mean de-biasing (with and without loess smoothing)
#' 
#' @param ... arguments passed to \code{\link{linmod}}
#' 
#' @keywords util
#' @rdname smooth
#' 
unbias <- function(...){
  return(linmod(..., formula=obs ~ offset(fcst) - 1, 
                smoothobs=FALSE, smooth=FALSE, recal=FALSE))
}
#' @rdname smooth
#' 
smoothobs <- function(...){
  return(linmod(..., formula=obs ~ offset(fcst) - 1, 
                smoothobs=TRUE, smooth=FALSE, recal=FALSE))
}
#' @rdname smooth
#' 
smooth <- function(...){
  return(linmod(..., formula=obs ~ offset(fcst) - 1, 
                smoothobs=TRUE, smooth=TRUE, recal=FALSE))
}

#' @rdname smooth
#' 
unbiasRecal <- function(...){
  return(linmod(..., formula=obs ~ offset(fcst) - 1, 
                smoothobs=FALSE, smooth=FALSE, recal=TRUE))
}
#' @rdname smooth
#' 
smoothobsRecal <- function(...){
  return(linmod(..., formula=obs ~ offset(fcst) - 1, 
                smoothobs=TRUE, smooth=FALSE, recal=TRUE))
}
#' @rdname smooth
#' 
smoothRecal <- function(...){
  return(linmod(..., formula=obs ~ offset(fcst) - 1, 
                smoothobs=TRUE, smooth=TRUE, recal=TRUE))
}
