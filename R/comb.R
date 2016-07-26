#' @name comb
#' @aliases combRecal
#' 
#' @title
#' Conditional Bias With Linear Time Trend
#' 
#' @description
#' Compute calibration for biases that are conditional on ensemble mean forecast
#' and exhibit a linear time trend
#' 
#' @param ... arguments passed to \code{\link{linmod}}
#'   
#' @seealso linmod
#'      
#' @examples
#' ## initialise forcast observation pairs
#' seasonal <- sin(seq(0,4,length=215))
#' signal <- outer(seasonal, rnorm(30), '+')
#' fcst <- array(rnorm(215*30*51), c(215, 30, 15)) + 
#'   2*c(signal)
#' obs <- array(rnorm(215*30, mean=2), c(215, 30)) +
#'   signal
#' fc.time <- outer(1:215, 1981:2010, function(x,y) as.Date(paste0(y, '-11-01')) - 1 + x)
#' fcst.debias <- biascorrection:::comb(fcst[,1:20,], 
#' obs[,1:20], fcst.out=fcst, fc.time=fc.time[,1:20], fcout.time=fc.time, span=0.5)
#' 
#' @rdname comb
#' @keywords util
comb <- function(...){
  return(linmod(..., formula=obs ~ fcst + fcst:poly(lead,3) + fcst:exp(-lead/5) + 
                  year + year:poly(lead,3) + year:exp(-lead/5), recal=FALSE))
}


#' @rdname comb
combRecal <- function(...){
  return(linmod(..., formula=obs ~ fcst + fcst:poly(lead,3) + fcst:exp(-lead/5) + 
                  year + year:poly(lead,3) + year:exp(-lead/5), recal=TRUE))
}

