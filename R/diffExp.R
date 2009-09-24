computeDiffExpT <- function(expData, annoData, groups, what = getColnames(expData, all = FALSE), 
                            verbose = FALSE, addTo = 0, ...) {
    if (length(groups) != length(what)) 
        stop("groups must be a factor with the same length as what")
    
    totals <- summarizeByAnnotation(expData, annoData, what = what, verbose = verbose, ...)

    idx <- split(1:length(groups), groups)
    levs <- groups[!duplicated(groups)]
 
    diffE <- function(l1, l2) {
        le1 <- sum(totals[, idx[[l1]]] + addTo)
        le2 <- sum(totals[, idx[[l2]]] + addTo)
        lam1 <- rowSums(totals[,idx[[l1]], drop = FALSE] + addTo)
        lam2 <- rowSums(totals[,idx[[l2]], drop = FALSE] + addTo)
        (log(lam1/le1) - log(lam2/le2))/sqrt(1/lam1 + 1/lam2)
    }
    
    res <- list()
    for (i in 1:(length(levs) - 1)) {
        for (j in (i+1):length(levs)) {
            res[[paste(levs[i], levs[j], sep = "-")]] <- diffE(levs[i], levs[j])
        }
    }

    lapply(res, function(a) {
      lst <- list("statistic" = a, "p.value" = 1 - pnorm(abs(a)))
    })
}



