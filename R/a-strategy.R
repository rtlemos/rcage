#' strategy class
#'
#' @import rcvirtual
#' @export rcage.strategy
#' @exportClass rcage.strategy
#'
rcage.strategy <- setRefClass(
  Class = 'rcage.strategy',
  contains = 'rcvirtual.strategy',
  methods = list(
  
    # --------------------------------------------------------------------------
    # - Initializer methods ----------------------------------------------------
    # --------------------------------------------------------------------------
    
    # --------------------------------------------------------------------------
    # - Set methods ------------------------------------------------------------
    # --------------------------------------------------------------------------
    
    # --------------------------------------------------------------------------
    # - Get methods ------------------------------------------------------------
    # --------------------------------------------------------------------------
    get.compute.U = function(c, s, phi) {
      # data
      nC <- .self$get.data(long.name = 'nC')
      aS <- .self$get.data(long.name = 'aS')
      phi <- .self$get.data(param.name = 'phi')
      # computing
      U <- c(-aS * (1 - phi), rep(1, nC - 1))
      # storing
      return(U)
    },
    
    get.compute.G = function(c, m, u, s, phi) {
      # data
      nC <- .self$get.data(long.name = 'nC')
      aM <- .self$get.data(long.name = 'aM')
      aU <- .self$get.data(long.name = 'aU')
      nS <- .self$get.data(long.name = 'nS')
      phi <- .self$get.data(param.name = 'phi')
      # computing
      G <- matrix(nrow = nC, ncol = nC, 0)
      G[2:nC, 1:(nC - 1)] <- 1
      G[1, (aM + 1):(aU + 1)] <- (1 - phi) / nS
      return(G)
    },
    
    get.compute.i = function(c) {
      # data
      nC <- .self$get.data(long.name = 'nC')
      # computing
      i <- rep(0, nC)
      return(i)
    },
    
    get.q = function(c, eta, alpha) {
      # data
      nC <- .self$get.data(long.name = 'nC')
      eta <- .self$get.data(param.name = 'eta')
      alpha <- .self$get.data(param.name = 'alpha')
      # computing
      age <- c(1:nC)
      q <- 1 / (1 + exp(-eta * (age - alpha)))
      return(q)
    },
    
    get.compute.v = function(t, h, o, c, G) {
      # data
      tvec <- .self$get.data(param.name = 't')
      h <- .self$get.data(param.name = 'h')
      o <- .self$get.data(long.name = 'o')
      G <- .self$get.data(param.name = 'G')
      nc <- .self$get.data(param.name = 'nC')
      omega <- .self$get.data(param.name = 'omega')
      
      # computing
      nt <- length(tvec)
      V <- Rinv <- Q <- vector('list', length = nt)
      O <- diag(o, nC, nC)
      Oinv <- diag(1 / o, nC, nC)
      eps <- diag(0.0001, nC, nC)
      W <- diag(omega, nC, nC)
      tG <- t(G)
      R <- G %*% diag(h, nC, nC) %*% tG + W
      tt <- 1
      Q[[tt]] <- R + O
      Rinv[[tt]] <- solve(R)
      vm <- solve(Rinv[[tt]] + Oinv)
      V[[tt]] <- 0.5 * vm + 0.5 * t(vm) + eps
      for (tt in (2:nt)) {
        R <- G %*% V[[tt - 1]] %*% tG + W
        Q[[tt]] <- R + O
        Rinv[[tt]] <- solve(R)
        vm <- solve(Rinv[[tt]] + Oinv)
        V[[tt]] <- 0.5 * vm + 0.5 * t(vm) + eps
      }
      v <- list(V = V, Rinv = Rinv, Q = Q)
      return(v)
    },
    
    get.compute.r = function(t, v, G, Y, E, q, c) {
      # data
      tvec <- .self$get.data(param.name = 't')
      v <- .self$get.data(param.name = 'v')
      G <- .self$get.data(param.name = 'G')
      Y <- .self$get.data(param.name = 'Y')
      E <- .self$get.data(param.name = 'E')
      q <- .self$get.data(param.name = 'q')
      nc <- .self$get.data(long.name = 'nC')
      # computing
      nt <- length(tvec)
      echi <- exp(chi)
      intercept <- rep(rho + chi, nC)
      m <- a <- e <- vector('list', length = nt)
      tt <- 1
      a[[tt]] <- G %*% rep(0, nC) - zeta * U - echi * q * E[1]
      f <- intercept + a[[tt]]
      e[[tt]] <- Y[, tt] - f
      m[[tt]] <- v$V[[tt]] %*% (v$Rinv[[tt]] %*% a[[tt]] + Y[, tt] / o)
      for (tt in 2:nt) {
        a[[tt]] <- G %*% rep(0, nC) - zeta * U - echi * q * E[tt - 1]
        f <- intercept + a[[tt]]
        e[[tt]] <- Y[, tt] - f
        m[[tt]] <- v$V[[tt]] %*% (v$Rinv[[tt]] %*% a[[tt]] + Y[, tt] / o)
      }
      r <- list(a = a, e = e, m = m)
      return(r)
    },
    
    get.compute.l = function(r, v, c) {
      # data
      r <- .self$get.data(param.name = 'r')
      v <- .self$get.data(param.name = 'v')
      nc <- .self$get.data(long.name = 'nC')
      # computing
      clog2pi <- nC * log(2 * pi)
      partial.llik <- mapply(r$e, v$Q, FUN = function(et, Qt) {
        -0.5 * (clog2pi + det(Qt) + as.numeric(crossprod(et, solve(Qt, et))))
      })
      l <- list(partial.llik = partial.llik, llik = llik)
      return(l)
    }
    
    # --------------------------------------------------------------------------
    # - Is methods -------------------------------------------------------------
    # --------------------------------------------------------------------------
    
    # --------------------------------------------------------------------------
    # - Internal methods -------------------------------------------------------
    # --------------------------------------------------------------------------
    
  )
)
