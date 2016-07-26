#' List Calibration Methods
#' 
#' List calibration methods available in \code{biascorrection}
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
  return(mout)
}