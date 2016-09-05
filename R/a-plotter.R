#' plotter class
#'
#' @import rcvirtual
#' @export rcage.plotter
#' @exportClass rcage.plotter
#'
rcage.plotter <- setRefClass(
  Class = 'rcage.plotter',
  contains = 'rcvirtual.ggplotter',
  methods = list(
    obs.vs.pred = function(by.age = TRUE){
      obs <- log(C / E[, 1])
      pred <- .self$strategy$get.value('r')$f
      if (by.age) {
        par(mfrow = c(2, ceiling(ncol(obs) / 2)))
        for (i in 1:ncol(obs)) matplot(cbind(obs[, i], pred[, i]), ty = 'l',
                                       ylab = 'log-cpue')
      } else {
        par(mfrow = c(1,1))
        n <- nrow(obs)
        f <- 0.5 * (1:n) / n
        cols <- c(grey(f), rgb(1, 1, f)) 
        matplot(0:(ncol(obs) - 1), t(rbind(obs, pred)), ty = 'l', 
                ylab = 'log-cpue', xlab = 'age', col = cols, lty = 1)
      }
    },
    
    selectivity.plot = function() {
      par(mfrow = c(1,1))
      s <- .self$strategy$get.value('s')
      plot(0:(length(s) - 1), s, ty = 'b', xlab = 'age', ylab = 'selectivity')
    },
    
    posterior.density.plots = function(par.names = NULL) {
      
      par.mean <- .self$strategy$model.res$par
      par.prec <- -.self$strategy$model.res$hessian
      par.sd <- 1 / sqrt(diag(par.prec))
      npar <- if (is.null(par.names)) length(par.mean) else length(par.names)
      par(mfrow = rep(npar, 2))
      id <- if (is.null(par.names)) 1:npar else {
        mapply(par.names, FUN = function(myname) which(names(par.mean) == myname))
      }
      for (i in id) for(j in id) {
        if (i == j) {
          xx <- par.mean[i] + seq(-4 * par.sd[i], 4 * par.sd[i], length = 30)
          if (length(.self$strategy$mcmc.res) > 0) {
            plot(density(.self$strategy$mcmc.res[i, ]), col = 'red', lty = 2,
                 xlab = names(par.mean[i]), ylab = '', main = '')
            lines(xx, dnorm(xx, par.mean[i], par.sd[i]))
          } else {
            plot(xx, dnorm(xx, par.mean[i], par.sd[i]), ty = 'l',
                 xlab = names(par.mean[i]), ylab = '')
          }
        } else {
          x1 <- par.mean[i] + seq(-3 * par.sd[i], 3 * par.sd[i], length = 20)
          x2 <- par.mean[j] + seq(-3 * par.sd[j], 3 * par.sd[j], length = 20)
          z <- mapply(x2, FUN = function(xx2) {
            mapply(x1, FUN = function(xx1) {
              err <- c(xx1 - par.mean[i], xx2 - par.mean[j])
              d <- exp(-0.5 * crossprod(err, par.prec[c(i, j), c(i, j)] %*% err))
              return(as.numeric(d))
            })
          })
          if (length(.self$strategy$mcmc.res) > 0) {
            est <- KernSmooth::bkde2D(
              x = t(gg$strategy$mcmc.res[c(i, j),]), 
              bandwidth = c(5 * KernSmooth::dpih(gg$strategy$mcmc.res[i,]),
                            5 * KernSmooth::dpih(gg$strategy$mcmc.res[j,])))
            contour(est$x1, est$x2, est$fhat, drawlabels = FALSE, 
                    col = 'red', lty = 2, nlevels = 5)
            points(gg$strategy$mcmc.res[i, ], gg$strategy$mcmc.res[j, ], 
                   pch = '.')
            contour(x1, x2, z, nlevels = 5, drawlabels = FALSE, add = TRUE)
          } else {
            contour(x1, x2, z, nlevels = 5, drawlabels = FALSE, add = TRUE)
          }
        }
      }
    },
    
    mcmc.traceplots = function(par.names = NULL) {
      all.pnames <- names(.self$strategy$model.res$par)
      if (is.null(par.names)) par.names <- all.pnames
      npar <- length(par.names)
      par(mfrow = c(ceiling(npar / 3), 3))
      for(param.name in par.names) {
        id <- which(param.name == all.pnames)
        plot(.self$strategy$mcmc.res[id, ], ty = 'l', ylab = param.name,
             xlab = 'MCMC iteration')
      }
    }
  )
)
