##
## This is Sequencing Specific.
##

##
## A simple function to compute coverage in a relatively naive way.
##

computeCoverage <- function (expData, annoData,
                             cutoff  = function(x, anno, group) x > 10, 
                             effort  = seq(1e+05, 5e+07, length = 20),
                             smooth  = function(probs) probs, 
                             groups  = rep("ALL", length(what)),
                             what    = getColnames(expData, all = FALSE),
                             totals  = summarizeExpData(expData, what = what, verbose = verbose),
                             ignoreStrand = FALSE,
                             verbose = FALSE, ...) 
{
    if (!is.list(annoData))
        stop("annoData must either be a data.frame, or a list of data.frames.")
        
    if (!is.data.frame(annoData) && length(groups) != 1) {
        stop("If looking across groups then only 1 sample at a time is allowed.")
    }

    if (length(groups) > 1) {
        if (length(groups) != length(what)) 
            stop("groups needs to be a factor/character of same length of columns.")
        groups <- factor(groups)
    }
    else {
        groups <- factor(rep(groups, length(what)))
    }
    groupIdxs <- split(1:length(what), groups)

    computeC <- function(annoData) {
        A <- summarizeByAnnotation(expData, annoData, what = what, verbose = verbose,
                                   ignoreStrand = ignoreStrand)
        L <- totals
        res <- mapply(function(idx, name) {
            A <- rowSums(A[, idx, drop = FALSE])
            if (is.null(dim(L)))
                L <- sum(L[idx])
            else 
                L <- sum(L[, idx, drop = FALSE])
            probs <- A/L
            probs <- c(probs, 1 - sum(probs))
            buckets <- 1:(length(A) + 1)
            probs <- smooth(probs)
            E <- outer(probs, effort)
            E <- E[-nrow(E), ]
            res <- matrix(NA, ncol = ncol(E), nrow = nrow(E))
            for (i in 1:nrow(E)) {
                res[i, ] <- cutoff(E[i, ], annoData[i, ], name, ...)
            }
            colSums(res)/nrow(res)
            }, groupIdxs, names(groupIdxs), SIMPLIFY = FALSE)

        return(do.call(cbind, res))
    }

    if (!is.data.frame(annoData)) {
        cc <- lapply(annoData, computeC)
        coverage <- do.call(cbind, cc)
        colnames(coverage) <- names(cc)
    }
    else {
        coverage <- computeC(annoData)
    }
    
    collapsedTotals <- tapply(totals, groups, sum)
    res <- list("coverage" = coverage, "effort" = effort, "totals" = collapsedTotals)
    class(res) <- "genominator.coverage"
    return(res)

}

plot.genominator.coverage <- function(x, type = "l", col = NULL,
                                      draw.totals = TRUE, draw.legend = TRUE, legend.location = NULL,
                                      ...) {
    ngroups <- ncol(x$coverage)
    if (is.null(col)) {
        col <- 1:ngroups
    }
    res <- matplot(x$effort, x$coverage, type = type, col = col, ...)

    if (draw.totals) {
        segments(x$totals, rep(0, ngroups), x$totals, rep(1, ngroups), col = col)
    }
    if (draw.legend) {
        if (is.null(legend.location))
            legend.location <- c(min(x$effort), max(x$coverage))
        legend(legend.location[1], legend.location[2], legend = colnames(x$coverage), fill = col)
    }

    return(invisible(res))
}
