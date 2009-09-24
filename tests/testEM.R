require(Genominator)

##
## Simulation Tests.
##
params <- list(mu = 4, pi = .7, phi = 2)

res.p <- t(replicate(1000, emZeroInflated(simData(1000, params$mu, params$pi),
                                          dfx = function(x, mu) {
                                            dpois(x, mu)
                                          })))
res.nb <- t(replicate(1000, emZeroInflated(simData(1000, params$mu,params$pi, phi = params$phi),
                                           dfx = function(x, mu) {
                                             dnbinom(x, mu = mu, size = 1/params$phi)
                                          })))

## I think that this makes sense mostly if I have a fair number
## of replicates. 
if (t.test(res.p[,"pi"], mu = params$pi)$p.value < .05)
  stop("Failed em Tests for poisson")
if (t.test(res.nb[,"pi"], mu = params$pi)$p.value < .05)
  stop("Failed em Tests for poisson")
if (t.test(res.p[,"mu"], mu = params$mu)$p.value < .05)
  stop("Failed em Tests for poisson")
if (t.test(res.nb[,"mu"], mu = params$mu)$p.value < .05)
  stop("Failed em Tests for poisson")

doTest <- function(v) {
  ss <- replicate(1000, {
    a1 <- sample(1:1000, size = 500)
    a2 <- setdiff(1:1000, a1)
    wilcox.test(v[a1], v[a2])$p.value
  })
  return(ss)
}
if(ks.test(doTest(res.p[,"pi"]), "punif")$p.value < .05)
  stop("Problem in test.")
if(ks.test(doTest(res.nb[,"pi"]), "punif")$p.value < .05)
  stop("Problem in test.")
if(ks.test(doTest(res.p[,"mu"]), "punif")$p.value < .05)
  stop("Problem in test.")
if(ks.test(doTest(res.nb[,"mu"]), "punif")$p.value < .05)
  stop("Problem in test.")

## par(mfrow=c(1,2))
## plot(density(res.nb[,1]), ylim = c(0,30))
## lines(density(res.p[,1]), col = "blue")
## abline(v = params$pi)
## plot(density(res.nb[,2]), ylim = c(0,10))
## lines(density(res.p[,2]), col = "blue")
## abline(v = params$mu)


