#' Template configuration for this package
#'
#' @import rcvirtual
#' @export rcage.conf
#' @exportClass rcage.conf
#'
rcage.conf <- setRefClass(
  Class = 'rcage.conf',
  contains = 'rcvirtual.conf',
  methods = list(
    construct = function(){

      callSuper()

      t.posix <- as.POSIXlt(c('1964/01/01', '2004/01/01'), tz = 'GMT')
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
      df$in.model[sk(c('Y', 'E', 't', 'c', 'm', 'u', 's',
                       'U', 'G', 'q', 'h', 'o', 'i', 'v', 'r', 'l', 
                       'phi', 'rho', 'xi', 'chi', 'theta'))] <- TRUE
      
      df$type[sk(c('Y', 'E', 't', 'c', 'm', 'u', 's'))] <- 'fixed'
      df$type[sk(c('U', 'G', 'q', 'h', 'o', 'i', 'v', 'r', 'l'))] <- 'derived'
      df$type[sk(c('phi', 'rho', 'xi', 'chi'))] <- 'unknown'
      df$type[sk('theta')] <- 'processes'
      df$type[sk('chi')] <- 'conjugate normal parameter'

      df$initial[sk('phi')] <- 0.5
      df$initial[sk('rho')] <- 5
      df$initial[sk('xi')] <- 0.1
      df$initial[sk('chi')] <- -4

      df$lbound[sk('phi')] <- 0
      df$lbound[sk('rho')] <- 0
      df$lbound[sk('xi')] <- 0
      df$lbound[sk('chi')] <- -10
      df$lbound[sk('t')] <- t.bnd[1]

      df$ubound[sk('phi')] <- 1
      df$ubound[sk('rho')] <- 15
      df$ubound[sk('xi')] <- 1
      df$ubound[sk('chi')] <- -1
      df$ubound[sk('t')] <- t.bnd[2]

      df$long.name[sk('Y')] <- 'log-cpue'
      df$long.name[sk('U')] <- 'system intercept'
      df$long.name[sk('G')] <- 'evolution matrix'
      df$long.name[sk('E')] <- 'fishing effort'
      df$long.name[sk('t')] <- 'time'
      df$long.name[sk('phi')] <- 'elasticity'
      df$long.name[sk('rho')] <- 'egg production ratio'
      df$long.name[sk('xi')] <- 'natural mortality rate'
      df$long.name[sk('theta')] <- 'log-abundance at age'
      df$long.name[sk('chi')] <- 'log-catchability'

      df$units[sk('t')] <- 'year'

      df$input.file[sk('Y')] <- paste0(tempdir(), '/Y.RData')
      df$input.file[sk('E')] <- paste0(tempdir(), '/E.RData')
      
      out.dir <- paste0(tempdir(), '/')
      df$output.file[sk('Y')] <- paste0(out.dir, 'cpue.csv')

      .self$parameters <- df
    }
  )
)
