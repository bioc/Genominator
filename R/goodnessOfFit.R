##
## Compute goodness of fit statistics.
##
setGeneric("regionGoodnessOfFit", function(obj, ...) standardGeneric("regionGoodnessOfFit"))

setMethod("regionGoodnessOfFit", "data.frame",
          definition = function(obj, denominator = colSums(obj), groups = rep("A", ncol(obj))) {
              if (length(groups) != ncol(obj) & length(groups) != length(denominator))
                  stop("Arguments 'groups', 'obj' and 'denominator' must all have the same length.")
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
          definition = function(obj, annoData, groups = rep("A", length(what)),
                                what = getColnames(obj, all = FALSE),
                                denominator = c("regions", "lanes"),
                                verbose = getOption("verbose")) {

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

plot.genominator.goodness.of.fit <- function(x, chisq = FALSE, plotCol = TRUE, qqline = FALSE,
                                             xlab = "theoretical quantiles",
                                             ylab = "observed quantiles", main,
                                             pch = 16, cex = .75, ...) {
    stats <- x$stats
    dfs <- x$dfs
    if(missing(main))
        main <- names(stats)
        
    a <- sqrt(length(stats))
    b <- a
    if (!is.integer(a)) {
        a <- floor(a)
        b <- length(stats)/a
    }
    op <- par(mfrow=c(a,b), no.readonly = TRUE)
    on.exit(par(op))
    
    mapply(function(mat, name, mainplot) {
        if (chisq) {
            yy <- mat[, "chi"]
            xx <- qchisq(ppoints(yy), df = dfs[name])
        } else {
            yy <- sort(mat[, "pval"])
            xx <- qunif(ppoints(yy), min = 0, max = 1)
        }
        qq <- qqplot(xx, yy, plot.it = FALSE)
        if (plotCol) {
            colors <- rep("black", length(yy))
            colors[floor(length(colors)*.95):length(colors)] <- "red"
            colors[floor(length(colors)*.99):length(colors)] <- "violet"
            colors[floor(length(colors)*.999):length(colors)] <- "orange"
        } else {
            colors <- "black"
        }
        plot(qq$x, qq$y, col = colors, xlab = xlab, ylab = ylab, main = mainplot, pch = pch, cex = cex, ...)
        if(qqline) {
            yyqt <- quantile(qq$y, c(0.25, 0.75))
            xxqt <- quantile(qq$x, c(0.25, 0.75))
            slope <- diff(yyqt) / diff(xxqt)
            int <- yyqt[1L] - slope * xxqt[1L]
            abline(int, slope, col = "blue")
        }
        abline(0, 1, col = "blue", lty = 2)
    }, stats, names(stats), main)
    invisible(NULL)
}

##
## This function operates on a vector
##
tabularGoodnessOfFit <- function(v, dfx = function(y) dpois(y, mean(v)),
                                 levels = 10, bins = NULL,
                                 verbose = getOption("verbose")) {
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


