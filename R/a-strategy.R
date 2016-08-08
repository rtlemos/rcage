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
    
    get.compute.theta = function(r, v) {
      'Log-abundance'
      
      theta <- list(type = 'unknown', mean = r, var = v)
      return(phi)
    },
    
    get.compute.phi = function() {
      'Elasticity'
      
      phi <- list(type = 'unknown', lb = 0, ub = 1)
      return(phi)
    },
    
    get.compute.rho = function() {
      'Log-carrying capacity'
      
      rho <- list(type = 'unknown', lb = 0, ub = 1)
      return(rho)
    },
    
    get.compute.xi = function() {
      'Natural mortality rate'
      
      xi <- list(type = 'unknown', lb = 0, ub = 0.4)
      return(xi)
    },
    
    get.compute.chi = function() {
      'Log-catchability'
      
      chi <- list(type = 'unknown', lb = -7, ub = -1)
      return(chi)
    },
    
    get.compute.Y = function() {
      'Log-CPUE data'
      
      return(2) #RTL TODO
    },
    
    get.compute.E = function() {
      'Effort data'
      
      return(1) #RTL TODO
    },
    
    get.compute.u = function() {
      'Species longevity (aU)'
      
      return(7)
    },
    
    get.compute.d = function(u) {
      'Number of cohorts (nC)'
      
      d <- u + 1
      return(d)
    },
    
    get.compute.m = function() {
      'Age at maturity (aM)'
      
      return(2)
    },
    
    get.compute.s = function(m, u) {
      'Mean spawning age (aS)'
      
      s <- 0.5 * (m + u)
      return(s)
    },
    
    get.compute.h = function() {
      'Initial system variance'
      
      h <- 1
      return(h)
    },
    
    get.compute.o = function() {
      'Measurement error variance'
      
      o <- 1
      return(o)
    },
    
    get.compute.t = function() {
      'Time bounds'
      
      t.posix <- as.POSIXlt(c('1964/01/01', '2004/01/01'), tz = 'GMT')
      t.bnd <- as.numeric(t.posix)
      return(t.bnd)
    },
    
    get.compute.U = function(d, s, phi) {
      'System intercept vector'
      
      U <- c(-s * (1 - phi), rep(1, d - 1))
      return(U)
    },
    
    get.compute.G = function(d, m, s, phi) {
      'Transition matrix'
      
      G <- matrix(nrow = d, ncol = d, 0)
      for(i in 2:d) G[i, i - 1] <- 1
      G[1, (m + 1):d] <- (1 - phi) / s
      return(G)
    },
    
    get.compute.i = function(d) {
      'Initial mean of state vector'
      
      i <- rep(0, d)
      return(i)
    },
    
    get.compute.q = function(u, eta, alpha) {
      'Selectivity function'

      age <- c(0:u)
      q <- 1 / (1 + exp(-eta * (age - alpha)))
      return(q)
    },
    
    get.compute.eta = function() {
      'Selectivity function slope'
      
      eta <- list(type = 'unknown', mean = 0, var = 1)
    },
    
    get.compute.omega = function() {
      'System error variance'
      
      omega <- list(type = 'unknown', mean = 0, var = 1)
    },
    
    get.compute.alpha = function() {
      'Selectivity function intercept'
      
      alpha <- list(type = 'unknown', mean = 0, var = 1)
    },
    
    get.compute.v = function(t, h, o, d, G, omega) {
      'System variance matrices'
      
      nt <- length(t)
      V <- Rinv <- Q <- vector('list', length = nt)
      O <- diag(o, d, d)
      Oinv <- diag(1 / o, d, d)
      eps <- diag(0.0001, d, d)
      W <- diag(omega, d, d)
      tG <- t(G)
      R <- G %*% diag(h, d, d) %*% tG + W
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
    
    get.compute.r = function(t, v, G, Y, E, U, q, d, chi, rho, xi) {
      'System mean vectors'
      
      nt <- length(t)
      echi <- exp(chi)
      intercept <- rep(rho + chi, d)
      m <- a <- e <- vector('list', length = nt)
      tt <- 1
      a[[tt]] <- G %*% rep(0, d) - xi * U - echi * q * E[1]
      f <- intercept + a[[tt]]
      e[[tt]] <- Y[, tt] - f
      m[[tt]] <- v$V[[tt]] %*% (v$Rinv[[tt]] %*% a[[tt]] + Y[, tt] / o)
      for (tt in 2:nt) {
        a[[tt]] <- G %*% rep(0, d) - xi * U - echi * q * E[tt - 1]
        f <- intercept + a[[tt]]
        e[[tt]] <- Y[, tt] - f
        m[[tt]] <- v$V[[tt]] %*% (v$Rinv[[tt]] %*% a[[tt]] + Y[, tt] / o)
      }
      r <- list(a = a, e = e, m = m)
      return(r)
    },
    
    get.compute.l = function(r, v, d) {
      'Likelihood'
      
      clog2pi <- d * log(2 * pi)
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
