##
## Simple implementation of a zero inflated Poisson / NB EM algorithm.
##
simData <- function(l, mu, pi, phi) {
    if(missing(phi)) 
        f <- function(n, mu) rpois(n, mu) 
    else 
        f <- function(n, mu) rnbinom(n, mu = mu, size = 1/phi)
    
    res <- numeric(l)
    zeros <- rbinom(1, prob = 1 - pi, size = l)
    nr <- (zeros + 1):length(res)
    res[nr] <- f(length(nr), mu)
    return(res)
}

emZeroInflated <- function(v, dfx, maxiter = 300, tol = 1e-5,
                           verbose = getOption("verbose")) {
    mu.start <- mean(v[v != 0]) 
    pi.start <- 1 - mean(v == 0)
    if (verbose) cat("mu:", mu.start, "pi:", pi.start, "\n")

    loop <- TRUE
    iter <- 1
    
    while (loop) {
        dd <- pi.start*dfx(v, mu.start)
        
        tau <- dd/(dd + (1-pi.start)*(v==0))
        tmp.mu <- sum(tau*v)/sum(tau)
        tmp.pi <- mean(tau)
        
        if (!(is.finite(tmp.mu) && is.finite(tmp.pi))) {
            mu.start <- NA
            pi.start <- NA
            break 
        }
        
        if ((abs(tmp.mu - mu.start) + abs(tmp.pi - pi.start)) < tol) {
            loop <- FALSE
        }
        pi.start <- tmp.pi
        mu.start <- tmp.mu
        iter <- 1 + iter 
        
        if (verbose) 
            if (iter %% 20 == 0) cat("mu:", mu.start, "pi:", pi.start, "\n")
        
        if (iter >= maxiter) {
            warning("Failed to converge in EM.")
            loop <- FALSE
        }
    }
    return(c(pi = pi.start, mu = mu.start, iter = iter))
}

zNB <- function(v, mu, phi, pi) {
    (v == 0)*(1 - pi) + pi*dnbinom(v, mu = mu, size = 1/phi)
}

zPois <- function(v, mu, pi) {
    (v == 0)*(1 - pi) + pi*dpois(v, mu)
}

