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
    },
    
    get.compute.phi = function() {
      'Elasticity'
      
      phi <- Uniform(name = 'elasticity', lb = 0.02, ub = 0.98)
    },
    
    get.compute.rho = function() {
      'Log-carrying capacity'
      
      # rho <- TruncatedNormal(name = 'log-K', lb = 3, ub = 30,
      #                        mean = 9, var = 9)
      rho <- Uniform(name = 'log-K', lb = 0, ub = 35)
      # rho <- Constant(9)
    },
    
    get.compute.xi = function() {
      'Natural mortality rate'
      
      #xi <- Uniform(name = 'natural mortality rate', lb = 0.01, ub = 2.5)
      xi <- TruncatedNormal(name = 'natural mortality rate', 
                            mean = 1.1, var = 0.1, lb = 0.7, ub = 3.06)
    },
    
    get.compute.chi = function() {
      'Log-catchability'
      
      #chi <- Uniform(name = 'log-catchability', lb = -20, ub = 0)
      chi <- TruncatedNormal(name = 'log-catchability', lb = -20, ub = 0,
                             mean = -10, var = 1)
      # chi <- Constant(-0.01)
    },
    
    get.compute.C = function() {
      'Catch data'
      
      C <- Constant(name = 'catches', path = paste0(tempdir(), '/C.RData'))
    },
    
    get.compute.E = function() {
      'Effort data'
      
      E <- Constant(name = 'effort', path = paste0(tempdir(), '/E.RData'))
    },
    
    get.compute.aR = function() {
      'Age at recruitment'
      
      aR <- Constant(0)
    },
    
    get.compute.aU = function() {
      'Species longevity'
      
      aU <- Constant(6)
    },
    
    get.compute.nC = function(aU) {
      'Number of cohorts (nC)'
      
      nC <- aU + 1
    },
    
    get.compute.aM = function() {
      'Age at maturity (aM)'
      
      aM <- Constant(2)
    },
    
    get.compute.aS = function(aM, aU) {
      'Mean spawning age'
      
      aS <- 0.5 * (aM + aU)
    },
    
    get.compute.nS = function(aM, aU) {
      'Mean spawning age'
      
      nS <- aU - aM + 1
    },
    
    get.compute.h = function() {
      'Initial system variance'
      
      h <- Constant(16)
    },
    
    get.compute.o = function() {
      'Measurement error variance'
      
      # log(1+0.04^2), where CV[x]=4% is the CV reported by Vauhgan et al
      # and the expression comes from the use of cpue~LogNormal (see wikipedia)
      o <- TruncatedNormal(name = 'cpue error variance', lb = -7, ub = 1,
                           mean = -6.44, var = 9)
    },
    
    get.compute.t = function() {
      'Time'
      
      t <- Time(lb = '1964/01/01', ub = '2004/01/01', step = 'year')
    },
    
    get.compute.U = function(nC, aS, phi) {
      'System intercept vector'
      
      U <- c(-aS * (1 - phi), rep(1, nC - 1))
    },
    
    get.compute.G = function(nC, aM, nS, phi) {
      'Transition matrix'
      
      G <- matrix(nrow = nC, ncol = nC, 0)
      for(i in 2:nC) G[i, i - 1] <- 1
      G[1, (aM + 1):nC] <- (1 - phi) / nS
      return(G)
    },
    
    get.compute.iniM = function(aR, aS, aU, phi, xi, chi) {
      'Initial mean of state vector'
      
      k <- (1 - phi) / phi * (aS - aR)
      echi <- exp(chi)
      iniM <- mapply(0:aU, FUN = function(age) {
        x <- -age * xi - k * echi * E[1, 1]
        out <- if (age < aR) x else x - (age - aR) * echi * E[1, 1]
        return(out)
      })
      #iniM <- mapply(0:aU, FUN = function(age) -xi * age )
      return(iniM)
    },
    
    get.compute.alpha = function() {
      'Selectivity function intercept'
      
      # alpha <- TruncatedNormal(name = 'selectivity intercept', 
      #                          lb = 0.5, ub = 2.5, mean = 2, var = 1)
      alpha <- Uniform(name = 'selectivity intercept', lb = 0.5, ub = 2.5)
      #alpha <- Constant(0.9)
    },
    
    get.compute.eta = function() {
      'Selectivity function slope'
      
      # eta <- TruncatedNormal(name = 'selectivity slope', lb = 4, ub = 20,
      #                        mean = 10, var = 4)
      eta <- Uniform(name = 'selectivity slope', lb = 2, ub = 10)
      #eta <- Constant(6)
    },
    
    get.compute.s = function(aU, eta, alpha) {
      'Selectivity function'

      s <- 1 / (1 + exp(-eta * (c(0:aU) - alpha)))
    },
    
    get.compute.omega = function() {
      'System error variance'
      
      #omega <- InverseGamma(name = 'system variance', rate = 3, shape = 1)
      #omega <- TruncatedNormal(name = 'system variance', 
      #                         mean = 0, var = 9, lb = -20, ub = 20)
      omega <- Uniform(name = 'system variance', lb = -20, ub = 20)
    },
    
    get.compute.v = function(t, h, o, nC, G, omega, C) {
      'System variance matrices'
      
      nt <- length(t)
      eo <- exp(o)
      V <- R.inv <- Q <- vector('list', length = nt)
      O <- diag(eo, nC, nC)
      eps <- diag(0.0001, nC, nC)
      W <- diag(exp(omega), nC, nC)
      tG <- t(G)
      R <- G %*% diag(h, nC, nC) %*% tG + W
      tt <- 1
      ok <- !is.na(C[tt, ])
      Q[[tt]] <- (R + O)[ok, ok]
      R.inv[[tt]] <- solve(R)
      vm <- solve(R.inv[[tt]] + diag(ok / eo, nC, nC))
      V[[tt]] <- 0.5 * vm + 0.5 * t(vm) + eps
      for (tt in (2:nt)) {
        R <- G %*% V[[tt - 1]] %*% tG + W
        ok <- !is.na(C[tt, ])
        Q[[tt]] <- (R + O)[ok, ok]
        R.inv[[tt]] <- solve(R)
        vm <- solve(R.inv[[tt]] + diag(ok / eo, nC, nC))
        V[[tt]] <- 0.5 * vm + 0.5 * t(vm) + eps
      }
      v <- list(V = V, R.inv = R.inv, Q = Q)
      return(v)
    },
    
    get.compute.r = function(t, v, G, C, E, U, s, nC, iniM, o, 
                             chi, rho, xi) {
      'System mean vectors'
      
      Y <- log(C / E[, 1])
      nt <- length(t)
      echi <- exp(chi)
      eo <- exp(o)
      intercept <- rep(rho + chi, nC)
      m <- a <- err <- f <- matrix(nrow = nrow(Y), ncol = ncol(Y))
      tt <- 1
      a[tt,] <- G %*% iniM - xi * U - echi * s * E[1, 1]
      f[tt,] <- intercept + a[tt,] + log(s)
      ok <- !is.na(Y[tt, ])
      err[tt, ] <- Y[tt, ] - f[tt, ]
      m[tt,] <- v$V[[tt]] %*% v$R.inv[[tt]] %*% a[tt,]
      m[tt, ok] <- m[tt, ok] + 
        v$V[[tt]][ok, ok] %*% ((Y[tt, ok] - intercept[ok] - log(s[ok])) / eo)
      for (tt in 2:nt) {
        a[tt,] <- G %*% m[tt - 1,] - xi * U - echi * s * E[tt - 1, 1]
        f[tt,] <- intercept + a[tt,] + log(s)
        ok <- !is.na(Y[tt, ])
        err[tt, ] <- Y[tt, ] - f[tt, ]
        m[tt,] <- v$V[[tt]] %*% v$R.inv[[tt]] %*% a[tt,]
        m[tt, ok] <- m[tt, ok] + 
          v$V[[tt]][ok, ok] %*% ((Y[tt, ok] - intercept[ok] - log(s[ok])) / eo)
      }
      r <- list(a = a, f = f, err = err, m = m)
      return(r)
    },
    
    get.compute.l = function(r, v, nC) {
      'Log-likelihood'
      
      clog2pi <- nC * log(2 * pi)
      partial.llik <- as.numeric(mapply(1:length(v$Q), FUN = function(tt) {
        ok <- !is.na(r$err[tt, ])
        -0.5 * (clog2pi + determinant(v$Q[[tt]], logarithm = TRUE)$modulus + 
                  crossprod(r$err[tt, ok], solve(v$Q[[tt]], r$err[tt, ok])))
      }))
      l <- LogLikelihood(sum(partial.llik))
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
