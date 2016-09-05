gg = rcage(autoconstruct = TRUE); gg$fit(mle = FALSE)
gg$strategy$model.res
gg$plotter$selectivity.plot()
gg$plotter$obs.vs.pred(TRUE)
gg$plotter$obs.vs.pred(FALSE)
gg$strategy$set.independence.sampler(nits = 10000)
gg$plotter$mcmc.traceplots()
gg$plotter$posterior.density.plots(par.names=c('alpha', 'eta'))
gg$plotter$posterior.density.plots(par.names=c('chi', 'rho'))
gg$plotter$posterior.density.plots(par.names=c('phi', 'xi'))
gg$plotter$posterior.density.plots(par.names=c('xi', 'rho'))
gg$strategy$get.posterior.correl()
gg$strategy$get.posterior.correl(method = 'mcmc')


.self = gg$strategy
aU <- 6
aR <- 0
aM <- 2
aS <- 0.5 * (aM + aU)
aB <- 0.856 * aR + 0.144 * aS
nR <- aU - aR + 1
phi <- .self$mcmc.res['phi',]
rho <- .self$mcmc.res['rho',]
xi <- .self$mcmc.res['xi',]
chi <- .self$mcmc.res['chi',]

R0 <- exp(rho - aR * xi)
N0 <- R0 * (1 - exp(-nR * (xi + 0))) / (1 - exp(-(xi + 0)))

H.MER <- (1 - phi) ^ ((1 - phi) / phi)
MER <- H.MER * phi * exp(rho)
F.MER <- -log(1 - phi) / (aS - aR)
R.MER <- exp(rho - aR * xi) * H.MER
N.MER <- R.MER * (1 - exp(-nR * (xi + F.MER))) / (1 - exp(-(xi + F.MER)))
E.MER <- F.MER / exp(xi) 
F.MSY <- F.MER * H.MER * (exp(1) - phi)
H.MSY <- (1 - phi) ^ (-(1 - exp(1) / phi) * (1 - phi) ^ (1 / phi))
R.MSY <- exp(rho - aR * xi) * H.MSY
N.MSY <- R.MSY * (1 - exp(-nR * (xi + F.MSY))) / (1 - exp(-(xi + F.MSY)))
E.MSY <- F.MSY / exp(xi)
BAR.MER <- (1 - phi) ^ (1 / phi - 0.856)

# F time series plot
tt <- .self$get.data(param.name = 't', field.name = 'year')
qtl <- quantile(chi, probs = c(0.025, 0.5, 0.975))
Ft <- outer(exp(qtl), E[, 1])
Fm <- quantile(F.MSY, probs = c(0.025, 0.5, 0.975))
plot(tt, rep(NA, length(tt)), ylim = c(min(Ft, Fm), max(Ft, Fm)),
     ylab = 'F', xlab = '')
rect(tt[1], Fm[1], tt[length(tt)], Fm[3], col = 'light grey', border = NA)
abline(h = Fm[2], lty = 2, col = grey(0.4))
lines(tt, Ft[1,], lty = 1, col = grey(0.5))
lines(tt, Ft[2,], lty = 1, col = 'black')
lines(tt, Ft[3,], lty = 1, col = grey(0.5))
