#' Template configuration for this package
#'
#' @import rcvirtual
#' @export rcage.conf
#' @exportClass rcage.conf
#'
rcage.conf <- setRefClass(
  Class = 'rcage.conf',
  contains = "rcvirtual.conf",
  methods = list(
    construct = function(){
      
      callSuper()
      
      .self$package.name <- "rcage"
      t.posix <- as.POSIXlt(c("1964/01/01", "2004/01/01"), tz = "GMT")
      t.bnd <- as.numeric(t.posix)
      
      #-------------------------------------------------------------------------
      # Daemon conf
      #-------------------------------------------------------------------------
      .self$daemon <- list(timeb = t.bnd)
      
      #-------------------------------------------------------------------------
      # Strategy conf
      #-------------------------------------------------------------------------
      .self$strategy <- list(timeb = t.bnd)
      
      #-------------------------------------------------------------------------
      # Plotter conf
      #-------------------------------------------------------------------------
      .self$plotter <- list(timeb = t.bnd)
      
      #-------------------------------------------------------------------------
      # Parameter conf
      #-------------------------------------------------------------------------
      df <- .self$get.conf.template()
      sk <- function(name) {
        as.numeric(mapply(name, FUN = function(name) which(df$name == name)))
      }
      df$type[sk(c("Y", "t"))] <- "fixed"
      df$type[sk("K")] <- "derived"
      df$type[sk("phi")] <- "unknown"
      df$type[sk("theta")] <- "processes"
      df$type[sk("chi")] <- "conjugate normal parameter"
      
      df$initial[sk("phi")] <- 0.1
      
      df$lbound[sk("phi")] <- 0
      
      df$ubound[sk("phi")] <- 1
      
      df$long.name[sk("Y")] <- "catches"
      df$long.name[sk("t")] <- "time"
      df$long.name[sk("K")] <- "carrying capacity"
      df$long.name[sk("phi")] <- "elasticity"
      df$long.name[sk("theta")] <- "log-abundance at age"
      df$long.name[sk("chi")] <- "log-catchability"
      
      df$is.spatial.avg[sk("Y")] <- TRUE
      
      df$store.in.ram[sk("Y")] <- TRUE
      
      df$units[sk("t")] <- "year"
      
      out.dir <- paste0(tempdir(), "/")
      df$output.file[sk("Y")] <- paste0(out.dir, "observations.csv")
      
      .self$parameters <- df[df$type != "-", ]
    }
  )
)