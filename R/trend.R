#' @name trend
#' @aliases trend
#' @aliases trendRecal
#' 
#' @title
#' Bias with linear time trend
#' 
#' @description
#' Computes mean de-biasing with linear time trend
#'  
#' @param ... arguments passed to \code{\link{linmod}}
#' 
#' @details
#' This bias correction method assumes that the bias can be decomposed 
#' into a stationary seasonal cycle (as in method \code{\link{smoothobs}}) 
#' and a linear time trend estimated from the residuals.
#' 
#' @examples
#' ## initialise forcast observation pairs
#' fcst <- array(rnorm(215*30*51), c(215, 30, 51)) + 
#' 0.5*sin(seq(0,4,length=215)) + 
#' rep(seq(0,1,length=30), each=215)
#' obs <- array(rnorm(215*30, mean=2), c(215, 30)) + 
#' sin(seq(0,4, length=215)) + 
#' rep(seq(0,3,length=30), each=215)
#' fc.time <- outer(1:215, 1981:2010, function(x,y) as.Date(paste0(y, '-11-01')) - 1 + x)
#' fcst.debias <- biascorrection:::trend(fcst[,1:20,], 
#' obs[,1:20], fcst.out=fcst, fc.time=fc.time[,1:20], fcout.time=fc.time, span=0.5)
#' 
#' @seealso linmod
#' 
#' @rdname trend
#' @keywords util
trend <- function(...){
  return(linmod(..., formula=obs ~ offset(fcst) + year + 
                  year:poly(lead,3) + year:exp(-lead/5), recal=FALSE))
}


#' @rdname trend
trendRecal <- function(...){
  return(linmod(..., formula=obs ~ offset(fcst) + year + 
                  year:poly(lead,3) + year:exp(-lead/5), recal=TRUE))
}
