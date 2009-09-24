## 
## XXX: There is a problem currently in the reference distribution for the
##      tabular goodness of fit test. The problem is related to binning it
##      seems and I don't know currently how to solve it. 
##
require(Genominator)

##
## Testing the Goodness of fit tools
##
L <- 8
N <- 100000
G <- 400

dta.1 <- matrix(rpois(N*L / 2, 10), nrow = N)
dta.2 <- matrix(rpois(N*L / 2, 12), nrow = N)
colnames(dta.1) <- paste("a", 1:(L/2), sep = "_")
colnames(dta.2) <- paste("b", 1:(L/2), sep = "_")
dta <- data.frame(chr = rep(1, N), strand = rep(1, N), location = 1:N, dta.1, dta.2)
ed <- importToExpData(dta, "/tmp/dd.db", "ttt", overwrite = TRUE)
genes <- matrix(sort(sample(1:N, size = G*2, replace = FALSE)), ncol = 2, byrow = TRUE)
genes <- data.frame(chr = 1, strand = 1, start = genes[,1], end = genes[,2])

nms <- c(rep("a", 4), rep("b", 4))
cols <- getColnames(ed, all = FALSE)
xx <- regionGoodnessOfFit(ed, genes, groups = nms)
yy <- regionGoodnessOfFit(ed, genes)
zz <- regionGoodnessOfFit(ed, genes, what = cols[1:4])
ww <- regionGoodnessOfFit(ed, genes, groups = nms[1:4]) # should be a warning

if (!all.equal(xx[[1]], zz[[1]]))
    stop("error!")

##
## Goodness of fit for just the regular old Poisson and NB
##
p.gf <- t(replicate(5000, { 
    tabularGoodnessOfFit(rpois(1000, 2))
}))

if (ks.test(p.gf, "punif")$p.value < .05)
    stop("something didn't work in the tabular goodness of fit with KS test.")

nb.gf <- t(replicate(1000, { 
    x <- rnbinom(1000, mu = 2, size = 1/1.5)
    unlist(tabularGoodnessOfFit(x, dfx = function(y)
                                dnbinom(y, mu = mean(x), size = 1/1.5), levels = 10))
}))

if (ks.test(nb.gf, "punif")$p.value < .01)
    stop("something didn't work in the tabular goodness of fit with KS test.")

##
## Some plotting code. 
##
## par(mfrow=c(1,2))
## hist(p.gf[,2], breaks = 20, prob = TRUE)
## boxplot(p.gf[,2] ~ p.gf[,3], varwidth = TRUE, notch = TRUE)
## qqplot(runif(3000), p.gf[,2]); abline(0,1,col="red")
## table(p.gf[,3], exclude = NULL)
## hist(p.nb[,2], breaks = 20, prob = TRUE)
## qqplot(runif(3000), p.nb[,2]); abline(0,1,col="red")
## boxplot(p.nb[,2] ~ p.nb[,3], varwidth = TRUE)

