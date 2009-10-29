##
## Compute goodness of fit statistics.
##
setGeneric("regionGoodnessOfFit", function(obj, ...) standardGeneric("regionGoodnessOfFit"))

setMethod("regionGoodnessOfFit", "data.frame",
          definition = function(obj, denominator = colSums(obj), groups = rep("A", ncol(obj))) {
              if (length(groups) != ncol(obj) & length(groups) != length(denominator))
                  stop("Arguments groups must be same length as obj and denominator.")
              sameSample <- split(1:ncol(obj), groups)
              denominator <- split(denominator, groups)

              res <- lapply(1:length(sameSample), function(i) {
                  idx <- sameSample[[i]]
                  lc <- denominator[[i]]
                  x <- obj[,idx]
                  lambda <- rowSums(x)/sum(lc)
                  E <- outer(lambda, lc)
                  chi <- rowSums((x - E)^2/E)
                  pval <- 1 - pchisq(chi, ncol(x) - 1)
                  cbind(chi = chi, pval = pval)
              })
              names(res) <- names(sameSample)
              
              if (length(groups) > 1) {
                  res <- res[groups[!duplicated(groups)]]
              }
              
              res <- list(stats = res, dfs = sapply(sameSample, length) - 1)
              class(res) <- "genominator.goodness.of.fit"
              return(res)
          })

setMethod("regionGoodnessOfFit", "ExpData",
          definition = function(obj, annoData, groups = rep("A", length(what)), what = getColnames(obj, all = FALSE),
                                denominator = c("regions", "lanes"), verbose = FALSE) {

              if (missing(groups)) {
                  groups <- rep("group", length(what))
              }
              
              if (length(groups) != length(what) & length(groups) != 1)
                  warning("Arguments groups and what are of different lengths.")
              regionSums <- summarizeByAnnotation(expData = obj, annoData = annoData,
                                                  what = what, verbose = verbose)
              if (!is.numeric(denominator)) {
                denominator <- match.arg(denominator)
                denominator <- switch(denominator,
                                      "regions" = colSums(regionSums),
                                      "lanes" = summarizeExpData(expData = obj,
                                        what = what, verbose = verbose))
              } else {
                stopifnot(length(denominator) == groups)
              }
              regionGoodnessOfFit(regionSums, denominator = denominator, groups = groups)
          })

plot.genominator.goodness.of.fit <- function(x, chisq = FALSE, plotCol = TRUE, sample = FALSE,
                                             nsamples = 5000, xlab = "theoretical quantiles",
                                             ylab = "observed quantiles", main = names(x),
                                             pch = 16, cex = .75, ...) {
    op <- par(no.readonly = TRUE)

    dfs <- x$dfs
    x <- x$stats
    args <- list(...)
    
    a <- sqrt(length(x))
    b <- a
    if (!is.integer(a)) {
        a <- floor(a)
        b <- length(x)/a
    }
    par(mfrow=c(a,b))
    
    g <- mapply(function(mat, name) {
        if (chisq) {
            y <- sort(mat[, "chi"])
            xx <- sort(rchisq(length(y), df = dfs[name]))
        } else {
            y <- sort(mat[, "pval"])
            xx <- seq(0, 1, length = length(y))
        }
        if (sample & length(y) > nsamples) {
            idx <- sort(sample(1:length(y), size = nsamples))
            y <- y[idx]
            xx <- xx[idx]
        }
        if (plotCol) {
            colors <- rep("black", length(y))
            colors[floor(length(colors)*.95):length(colors)] <- "red"
            colors[floor(length(colors)*.99):length(colors)] <- "violet"
            colors[floor(length(colors)*.999):length(colors)] <- "orange"
            plot(xx, y, xlab = xlab , ylab = ylab, main = name, col = colors, pch = pch, cex = cex, ...)
        }
        else {
            plot(xx, y, xlab = xlab, ylab = ylab, main = name, pch = pch, cex = cex, ...)
        }
        abline(0, 1, col = "blue")
    }, x, names(x))

    par(op)
}

##
## This function operates on a vector
##
tabularGoodnessOfFit <- function(v, dfx = function(y) dpois(y, mean(v)),
                                 levels = 10, bins = NULL, verbose = FALSE) {
    levs <- unique(floor(quantile(v, seq(0, 1, length = levels))))

    if (length(levs) <= 2)
        return(rep(NA, 3))

    levs <- c(-1, levs, max(v) + 1)
    tbl <- table(cut(v, levs))

    # collapse low counts.
    if (tbl[length(tbl)] == 0) {
        levs <- levs[-length(levs)]
        tbl <- table(cut(v, levs))
    }

    if (!is.null(bins)) {
      tbl <- table(cut(v, bins))
      levs <- bins
    }

    if (any(tbl < 5)) {
        warning("Low counts in table, possibly bad approximation with chi squared.")
    }

    dlevs <- diff(levs) - 1
    dlevs[length(dlevs)] <- dlevs[length(dlevs)] + 1
    probs <- numeric(length(tbl))

    for (i in 1:length(dlevs)) {
        start <- levs[i] + 1
        end <- start + dlevs[i]
        probs[i] <- sum(dfx(start:end))
    }
    probs[i] <- probs[i] + (1 - sum(probs))

    if (sum(tbl) != length(v))
        stop("Unforseen problem in tabularGoodnessOfFit: Programming Error!")

    stat <- sum((tbl - probs*sum(tbl))^2/(probs*sum(tbl)))

    return(list(statistic = stat, p.value = 1 - pchisq(stat, df = length(probs) - 2),
                bins = length(probs)))
}


