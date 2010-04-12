computePrimingWeights <- function(aln, biasedIndex = 1:2, unbiasedIndex = 24:29, weightsLength = 7L){
    if(!isTRUE(class(aln) == "AlignedRead"))
        stop("argument 'aln' needs to be of class 'AlignedRead'")
    if(!isTRUE(is.integer(weightsLength) && weightsLength > 0))
        stop("argument 'weightsLength' needs to be a single, positive integer")
    if(!isTRUE(is.integer(biasedIndex) && all(biasedIndex > 0)))
        stop("argument 'biasedIndex' needs to be positive integers")
    if(!isTRUE(is.integer(unbiasedIndex) && all(unbiasedIndex > 0)))
        stop("argument 'unbiasedIndex' needs to be positive integers")
    ## More argument checking
    p_base <- rep(0, 4^weightsLength)
    names(p_base) <- mkAllStrings(c("A", "C", "G", "T"), width = weightsLength)
    p.compute <- function(i) {
        tab <- tables(subseq(sread(aln), start = i, width = weightsLength), n = NULL)$top
        out <- p_base
        common.names <- intersect(names(out), names(tab))
        out[common.names] <- tab[common.names]
        out
    }
    p_biased <- rowMeans(do.call(cbind, lapply(biasedIndex, p.compute)))
    p_unbiased <- rowMeans(do.call(cbind, lapply(unbiasedIndex, p.compute)))
    if(any(is.na(p_biased) | p_biased < sqrt(.Machine$double.eps)))
        warning("weights might be degenerate")
    if(any(is.na(p_unbiased) | p_unbiased < sqrt(.Machine$double.eps)))
        warning("weights might be degenerate")
    weights <- p_unbiased / p_biased
    weights
}

addPrimingWeights <- function(aln, weights = NULL, overwrite = FALSE, ...) {
    if("weights" %in% varLabels(alignData(aln)) && !overwrite)
        stop("argument 'aln' already has weights, use overwrite = TRUE")
    if(is.null(weights)) {
        weights <- computePrimingWeights(aln, ...)
    }
    weightsLength <- nchar(names(weights[1]))
    readWeights <- as.numeric(weights[as.character(subseq(sread(aln), start = 1, width = weightsLength))])
    if("weights" %in% varLabels(alignData(aln))) {
        alignData(aln)$"weights" <- readWeights
        varMetadata(alignData(aln))["weights", "labelDescription"] <- "Priming weights"
    } else {
        aln@alignData <- AlignedDataFrame(data = cbind(pData(alignData(aln)), weights = readWeights),
                                          metadata = rbind(varMetadata(alignData(aln)),
                                          data.frame(labelDescription = c(weights = "Priming weights"))))
    }
    aln
}

fixWeights <- function(wt, verbose = FALSE) {
    if(verbose) cat("number of weights", length(wt),
                    "\n", fill = TRUE)
    if(any(is.na(wt))) {
        wh.na <- which(is.na(wt))
        if(verbose) cat("number of na weights", length(wh.na),
                        "\n", fill = TRUE)
        wt[wh.na] <- 0
        if(any(wt[-wh.na] == 0)) {
            if(verbose) cat("number of zero weights", sum(wt[-wh.na] == 0),
                            "\n", fill = TRUE)
            wt[wt == 0] <- 1/100
            wt[wh.na] <- 0
        }
    } else {
        if(any(wt == 0)) {
            if(verbose) cat("number of zero weights", sum(wt == 0),
                            "\n", fill = TRUE)
            wt[wt == 0] <- 1/100
        }
    }
    if(any(wt == Inf)) {
        if(verbose) cat("number of Inf weights", sum(wt == Inf),
                        "\n", fill = TRUE)
        wt[wt == Inf] <- 100
    }
    wt
}
