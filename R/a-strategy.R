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
    
    get.compute.theta = function(r, v, t, G) {
      'Log-abundance'
      
      theta <- GaussianState(name = 'log-abundance', 
                             full.conditional.mean = r$m, 
                             full.conditional.var = v$V,
                             prior.mean = r$a,
                             prior.precision = v$R.inv,
                             evolution.matrix = G,
                             n.instants = length(t))
      return(theta)
    },
    
    get.compute.phi = function() {
      'Elasticity'
      
      phi <- Uniform(name = 'elasticity', lb = 0, ub = 1)
      return(phi)
    },
    
    get.compute.rho = function() {
      'Log-carrying capacity'
      
      rho <- Uniform(name = 'log-K', lb = 0, ub = 1)
      return(rho)
    },
    
    get.compute.xi = function() {
      'Natural mortality rate'
      
      xi <- Uniform(name = 'natural mortality rate', lb = 0, ub = 0.4)
      return(xi)
    },
    
    get.compute.chi = function() {
      'Log-catchability'
      
      chi <- Uniform(name = 'log-catchability', lb = -7, ub = -1)
      return(chi)
    },
    
    get.compute.Y = function() {
      'Log-CPUE data'
      
      Y <- Constant(name = 'log-cpue', path = paste0(tempdir(), '/Y.RData'))
      return(Y)
    },
    
    get.compute.E = function() {
      'Effort data'
      
      E <- Constant(name = 'effort', path = paste0(tempdir(), '/E.RData'))
      return(E)
    },
    
    get.compute.aU = function() {
      'Species longevity'
      
      aU <- Constant(6)
      return(aU)
    },
    
    get.compute.nC = function(aU) {
      'Number of cohorts (nC)'
      
      nC <- aU + 1
      return(nC)
    },
    
    get.compute.aM = function() {
      'Age at maturity (aM)'
      
      aM <- Constant(2)
      return(aM)
    },
    
    get.compute.aS = function(aM, aU) {
      'Mean spawning age'
      
      aS <- 0.5 * (aM + aU)
      return(aS)
    },
    
    get.compute.h = function() {
      'Initial system variance'
      
      h <- Constant(1)
      return(h)
    },
    
    get.compute.o = function() {
      'Measurement error variance'
      
      o <- Constant(1)
      return(o)
    },
    
    get.compute.t = function() {
      'Time'
      
      t <- Time(lb = '1964/01/01', ub = '2004/01/01', step = 'year')
      return(t)
    },
    
    get.compute.U = function(nC, aS, phi) {
      'System intercept vector'
      
      U <- c(-aS * (1 - phi), rep(1, nC - 1))
      return(U)
    },
    
    get.compute.G = function(nC, aM, aS, phi) {
      'Transition matrix'
      
      G <- matrix(nrow = nC, ncol = nC, 0)
      for(i in 2:nC) G[i, i - 1] <- 1
      G[1, (aM + 1):nC] <- (1 - phi) / aS
      return(G)
    },
    
    get.compute.i = function(nC) {
      'Initial mean of state vector'
      
      i <- rep(0, nC)
      return(i)
    },
    
    get.compute.q = function(aU, eta, alpha) {
      'Selectivity function'

      age <- c(0:aU)
      q <- 1 / (1 + exp(-eta * (age - alpha)))
      return(q)
    },
    
    get.compute.eta = function() {
      'Selectivity function slope'
      
      eta <- Normal(name = 'selectivity slope', mean = 0, var = 1)
    },
    
    get.compute.omega = function() {
      'System error variance'
      
      omega <- InverseGamma(name = 'system variance', rate = 3, shape = 1)
    },
    
    get.compute.alpha = function() {
      'Selectivity function intercept'
      
      alpha <- Normal(name = 'selectivity intercept', mean = 0, var = 1)
    },
    
    get.compute.v = function(t, h, o, nC, G, omega) {
      'System variance matrices'
      
      nt <- length(t)
      V <- R.inv <- Q <- vector('list', length = nt)
      O <- diag(o, nC, nC)
      Oinv <- diag(1 / o, nC, nC)
      eps <- diag(0.0001, nC, nC)
      W <- diag(omega, nC, nC)
      tG <- t(G)
      R <- G %*% diag(h, nC, nC) %*% tG + W
      tt <- 1
      Q[[tt]] <- R + O
      R.inv[[tt]] <- solve(R)
      vm <- solve(R.inv[[tt]] + Oinv)
      V[[tt]] <- 0.5 * vm + 0.5 * t(vm) + eps
      for (tt in (2:nt)) {
        R <- G %*% V[[tt - 1]] %*% tG + W
        Q[[tt]] <- R + O
        R.inv[[tt]] <- solve(R)
        vm <- solve(R.inv[[tt]] + Oinv)
        V[[tt]] <- 0.5 * vm + 0.5 * t(vm) + eps
      }
      v <- list(V = V, R.inv = R.inv, Q = Q)
      return(v)
    },
    
    get.compute.r = function(t, v, G, Y, E, U, q, nC, o, chi, rho, xi) {
      'System mean vectors'
      
      nt <- length(t)
      echi <- exp(chi)
      intercept <- rep(rho + chi, nC)
      m <- a <- e <- Y
      tt <- 1
      a[tt,] <- G %*% rep(0, nC) - xi * U - echi * q * E[1, 1]
      f <- intercept + a[tt,]
      e[tt,] <- Y[tt, ] - f
      m[tt,] <- v$V[[tt]] %*% (v$R.inv[[tt]] %*% a[tt,] + Y[tt, ] / o)
      for (tt in 2:nt) {
        a[tt,] <- G %*% rep(0, nC) - xi * U - echi * q * E[tt - 1, 1]
        f <- intercept + a[tt,]
        e[tt,] <- Y[tt, ] - f
        m[tt,] <- v$V[[tt]] %*% (v$R.inv[[tt]] %*% a[tt,] + Y[tt, ] / o)
      }
      r <- list(a = a, e = e, m = m)
      return(r)
    },
    
    get.compute.l = function(r, v, nC) {
      'Likelihood'
      
      clog2pi <- nC * log(2 * pi)
      partial.llik <- mapply(1:length(v$Q), FUN = function(tt) {
        -0.5 * (clog2pi + determinant(v$Q[[tt]], logarithm = TRUE)$modulus + 
                  as.numeric(crossprod(r$e[tt,], solve(v$Q[[tt]], r$e[tt,]))))
      })
      l <- sum(partial.llik)
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
