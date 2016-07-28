#' Strategy class
#'
#' @import rcvirtual
#' @export rcage.strategy
#' @exportClass rcage.strategy
#'
rcage.strategy <- setRefClass(
  Class = "rcage.strategy",
  contains = "rcvirtual.strategy"
)