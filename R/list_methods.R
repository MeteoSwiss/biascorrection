#' List calibration methods
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
  desc <- gsub("^[a-z_]* *", "", info)
  m.i <- which(! methods %in% c('month', 'sloess', 'debias', 'list_methods', 'cor'))
  mout <- data.frame(METHODS=methods, DESCRIPTION=desc)[m.i,]
  return(mout)
}