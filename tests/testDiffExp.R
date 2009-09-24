require(Genominator)

R <- 2
N <- 1000
D <- 50
AL <- 3e6
EFFECT <- 2
groups <- rep(c("a", "b"), each = R)

iter <- 50
RES <- array(dim=c(2,2,iter))

for (jj in 1:iter) {
  BETAS <- rnorm(R*2, AL, AL^(1/1.2))
  LAMBDAS.1 <- LAMBDAS.2 <- runif(N, 0, 1/(10*N))
  LAMBDAS.2[1:D] <- LAMBDAS.2[1:D]*EFFECT
  LAMBDAS.1 <- c(LAMBDAS.1, 1 - sum(LAMBDAS.1))
  LAMBDAS.2 <- c(LAMBDAS.2, 1 - sum(LAMBDAS.2))
  if (any(LAMBDAS.1 < 0 | LAMBDAS.2 < 0))
    browser()
  DTA <- apply(cbind(outer(LAMBDAS.1, BETAS[1:R]),
                     outer(LAMBDAS.2, BETAS[(R+1):(R*2)])), c(1,2),
               function(mm) { rpois(1, mm) })
  a <- rowSums(DTA[,1:4])/sum(DTA[,1:4])
  b <- rowSums(DTA[,5:8])/sum(DTA[,5:8])
  num <- log(a/b)
  den <- sqrt(1/rowSums(DTA[,1:4]) + 1/rowSums(DTA[,5:8]))
  stat <- num/den
  gg <- sapply(1:nrow(DTA), function(i) {
    coefficients(summary(glm(DTA[i,] ~ colSums(DTA) + as.factor(groups),
                             family = "poisson")))[3,4]
  })
  
  RES[,,jj] <- rbind("glm" = c(TP = sum(gg[1:D] < .05)/D, FP = sum(gg[(D+1):N] < .05)/(N-D)),
                     "mp" = c(TP = sum(abs(stat[1:D]) > 1.96, na.rm = TRUE)/D,
                       FP = sum(abs(stat[(D+1):N]) > 1.96, na.rm = TRUE)/(N-D)))
}
dimnames(RES) <- list(c("glm", "mp"), c("TP", "FP"), 1:(dim(RES)[3]))
boxplot(as.data.frame(t(RES[1:2,2,])))
boxplot(as.data.frame(t(RES[1:2,1,])))


dta <- c(rep(0, 1000), rep(1, 1000))
Y <- sapply(exp(1 + dta*2), function(x) rpois(1, x))
coefficients(summary(glm(Y ~ factor(dta), family = "poisson")))

Yc <- tapply(Y, dta, sum)
coefficients(summary(glm(Yc ~ factor(c(0,1)), family = "poisson")))

