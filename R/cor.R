#' Unbiased correlation
#' 
#' This function computes the unbiased population estimate of the correlation
#' 
#' @param x a numeric vector
#' @param y a numeric vector
#' @param method a character string indicating which correlation coefficient 
#' is to be computed. The unbiased version is only available for 
#' \code{"pearson"} (the default).
#' @param ... additional arguments passed to \code{\link{cor}}
#' 
#' The correlation function introduced here uses the bias correction from
#' Olkin and Pratt (1958) as shown below where \eqn{\hat{\rho}} is the 
#' unbiased correlation estimator and \eqn{r} is the standard sample 
#' correlation. The unbiased correlation is
#' \deqn{\hat{\rho} = r \left[ 1 + \frac{1 - r^2}{2(n - 3)} \right]}
#' 
#' @references Olkin, I. and Pratt, J.W. (1958). Unbiased estimation 
#'     of certain correlation coefficients. Annals of Mathematical 
#'     Statistics, 29, 201-211.
#'     
#' @examples
#' x <- rnorm(20)
#' y <- rnorm(20) + 0.5*x
#' biascorrection:::cor(x,y)
#' stats::cor(x,y)
#' 
cor <- function(x,y, method='pearson', ...){
  if (method == 'pearson'){
    r <- stats::cor(x, y, method=method, ...)
    n <- sum(!is.na(x) & !is.na(y))
    rho <- r * (1 + (1 - r**2) / 2 / (n - 3))
    
  } else {
    rho <- stats::cor(x, y, method=method, ...)
  }
  return(rho)
}