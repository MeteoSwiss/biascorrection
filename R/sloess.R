#' Compute loess fit
#' 
#' @description
#' This function computes a loess fit on the given vector if more than 5 non-missing
#' values are available. Else, the raw values are returned.
#' 
#' @param x input vector
#' @param ... additional parameters handed over to loess
#' 
#' @keywords util
sloess <- function(x, ...){
  if (sum(!is.na(x)) > 5){
    xout <- loess(x ~ seq(along=x), ...)$fit
  } else {
    xout <- x
  }
  return(xout)
} 