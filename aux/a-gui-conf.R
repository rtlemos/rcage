#' Edit configuration with a Shiny app
#'
#' @import rcvirtual
#' @import shiny
#' @export rcage.guiconf
#' @exportClass rcage.guiconf
#'
rcage.guiconf <- setRefClass(
  Class = 'rcage.guiconf',
  contains = 'rcvirtual.guiconf'
)