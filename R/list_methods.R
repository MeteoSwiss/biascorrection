#' @name list_methods
#' @aliases print.biascorrection
#' 
#' @title List Calibration Methods
#' 
#' @description 
#' \code{list_methods} lists calibration methods available in \code{biascorrection}
#' and \code{print.biascorrection} is used to nicely format these for display.
#' 
#' @examples
#' list_methods()
#' 
#' @export
list_methods <- function(){
  info <- library(help='biascorrection')$info[[2]]
  methods <- gsub(" .*", "", info)
  desc <- gsub("^[a-zA-Z_]* *", "", info)
  m.i <- which(methods %in% c('monmean', 'biascorrection', 'sloess', 'debias', 'list_methods', 'cor', 'debiasApply'))
  m.i <- c(m.i, which(methods[m.i + 1] == ""), which(desc == "Forecast Data Sets"))
  mout <- data.frame(METHODS=methods, DESCRIPTION=desc)[-m.i,]
  class(mout) <- c('biascorrection')
  return(mout)
}

#' @rdname list_methods
#' @param df data frame of methods to be used with \code{\link{debias}}
#' @export
print.biascorrection <- function(df){
  df$METHODS <- paste0(df$METHODS, ' ')
  class(df) <- 'data.frame'
  print.data.frame(df, right=FALSE, row.names=FALSE)
}