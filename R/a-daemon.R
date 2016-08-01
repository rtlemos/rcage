#' rcage: daemon for this package
#'
#' @import rcvirtual
#' @export rcage
#' @exportClass rcage
#'
rcage <- setRefClass(
  Class = 'rcage',
  contains = c('rcvirtual.daemon')
)
