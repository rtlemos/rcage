#' plotter class
#'
#' @import rcvirtual
#' @export rcage.plotter
#' @exportClass rcage.plotter
#'
rcage.plotter <- setRefClass(
  Class = 'rcage.plotter',
  contains = 'rcvirtual.ggplotter'
)
