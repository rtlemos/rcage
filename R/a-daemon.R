#' rcage: daemon for this package
#' 
#' @import rcvirtual
#' @export rcage
#' @exportClass rcage
#'
rcage <- setRefClass(
  Class = "rcage",
  contains = c("rcvirtual.daemon"),
  methods = list(
    initialize = function(uconf = NULL, verbose = TRUE, autoconstruct = TRUE){
      "Initialize an rcage daemon object"
      
      p.name <- 'rcage'
      if (is.null(uconf)) {
        uconf <- rcage.conf$new(package.name = p.name, 
                                object.name = paste0(p.name, ".conf"))
        uconf$construct()
      }        
      
      callSuper(package.name = p.name,
       object.name = ".daemon",
       conf = uconf)
    }
  )
)
  