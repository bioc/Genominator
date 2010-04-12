toAnnoData <- function(annoTrack) {
    cbind(chr = rep(annoTrack@chr, nrow(annoTrack@regions)), start = annoTrack@regions[, "start"],
          end = annoTrack@regions[, "end"], strand = rep(annoTrack@strand, nrow(annoTrack@regions)),
          annoTrack@regions[,-match(c("start", "end"), colnames(annoTrack@regions))])
}

makeAnnoFactory <- function(x,
                            dp = DisplayPars(),
                            featureColumnName = "feature",
                            groupColumnName = NULL,
                            idColumnName = NULL,
                            expandAround = 1000,
                            chrFunction = function(yy) yy,
                            strandFunction = function(yy) yy) {

    if (!require(GenomeGraphs))
        stop("GenomeGraphs must be installed to use this functionality.")

    cChrFunction <- function(y) as.character(chrFunction(y))
    cStrandFunction <- function(y) as.character(strandFunction(y))
    
    regionAnno.Biomart <- function(chr, start, end, strand) { 
        geneRegionBiomart(chr, start - expandAround, end + expandAround,
                          strand, biomart, dp,
                          chrFunction = cChrFunction,
                          strandFunction = cStrandFunction)
    }
    
    regionAnno.AnnoData <- function(chr, start, end, strand) {
        mat <- annoData[annoData[,"start"] >= start - expandAround &
                        annoData[,"end"] <= end + expandAround &
                        annoData[,"chr"] == chr &
                        annoData[,"strand"] == strand, ]
        if (nrow(mat) == 0) {
            return(NULL)
        }

        mat <- makeGeneRepresentation(mat, type = "Ugene", gene.id = NULL)
        

        ## XXX: This is bad - it implies that I have to know some things about
        ##      how genome graphs works *internally*
        groups <- if (is.null(groupColumnName))
            rownames(mat)
        else
            mat[,groupColumnName]

        IDs <- if (is.null(idColumnName)) 
            rownames(mat)
        else
            mat[,idColumnName]

        features <- if(! featureColumnName %in% colnames(mat))
            "protein_coding"
        else
            mat[, featureColumnName]
        
        regions <- data.frame(start = mat[,"start"], end = mat[, "end"],
                              feature = features,
                              group = groups, ID = IDs)
        ## Fix for transcript based annotation
        regions <- regions[! duplicated(regions[, c("start", "end", "feature", "group")]), ]
        regions <- regions[order(regions$start),]
        makeAnnotationTrack(chr = chr, strand = strand, regions =
                            regions, dp = dp)
    }


    if(class(x) == "Mart") {
        biomart <- x
        return(regionAnno.Biomart)
    }
    if(validAnnotation(x)) {
        annoData <- x
        return(regionAnno.AnnoData)
    }
    stop("Unknown object as the value of 'x'")
}


makeRegionPlotter <- function(expDataLists, annoFactory = NULL) {
    if (!require(GenomeGraphs)) {
        stop("GenomeGraphs must be installed to use this functionality")
    }

    ##-- return a function which gets called on a region to produce a
    ##   GenomeGraphs object.
    fxs <- lapply(expDataLists, function(lst) {
        dp <- if (is.null(lst$dp)) {
            DisplayPars()
        } else {
            lst$dp
        }
        fx <- if (is.null(lst$fx)) {
            function(x) { return(x) } 
        } else {
            lst$fx
        }
        segFactory <- if (is.null(lst$segFactory)) {
            function(chr, start, end) { return(NULL) }
        } else {
            lst$segFactory
        }

        if (!is.null(lst$expandBy))
            expandBy <- as.numeric(lst$expandBy)
        else
            expandBy <- 0
        
        function(chr, start, end) {
            strand <- if (is.null(lst$strand)) {
                0
            } else {
                lst$strand
            }
            
            trackOverlay <- segFactory(chr, start, end)
            
            res <- getRegion(lst$expData, chr, start - expandBy, end + expandBy, strand, what = c("location", lst$what))

            if (is.null(res) || nrow(res) == 0) {
                xx <- integer(1)
                yy <- numeric(1)
                xx[1] <- 1
                yy[1] <- NA
                return(makeBaseTrack(base = xx, value = yy, dp = dp, trackOverlay = trackOverlay))
            }

            if (is.null(lst$class)) {
                klass <- "BaseTrack"
            } else {
                klass <- lst$class
            }
            
            
            switch (klass,
                    BaseTrack = {
                        makeBaseTrack(base = res[, "location"], value = fx(res[, lst$what]), dp = dp,
                                      trackOverlay = trackOverlay)
                    },
                    GenericArray = {
                        makeGenericArray(probeStart = res[, "location"],
                                         intensity = fx(as.matrix(res[, lst$what, drop = FALSE])),
                                         dp = dp, trackOverlay = trackOverlay)
                    },
                {
                    stop("Unknown class specified.")
                })
        }
    })

    function(chr, start, end, overlays = NULL, title = NULL, ...) {
        if (!is.null(annoFactory)) {
          gr <- lapply(c(-1, 0, 1), function(str) {
              annoFactory(chr, start, end, strand = str)
          })
          names(gr) <- c("+", "0", "-")
        } else {
          gr <- NULL
        }
        lst <- lapply(fxs, function(fx) {
            fx(chr, start, end)
        })
        unstranded <- if(is.null(gr[["0"]])) {
            NULL
        } else {
            list(" " = gr[["0"]], makeGenomeAxis())
        }
        ann <- c("+" = gr[["+"]], makeGenomeAxis(), " " = unstranded, "-" = gr[["-"]])
        lst <- Filter(function(x) !is.null(x), c(ann, lst))

        gdPlot(c(lst, title), minBase = start, maxBase = end, overlays = overlays, ...)
    }
}

## makeCoveragePlotter <- function(eData, what, readLength = 1, dp = DisplayPars()) {
##     ## use both anno as char and data.frame and also accept chr, start, end
##     makeTracks <- function(anno, extend = 0) {
##         if(nrow(anno) != 1) stop("")
##         anno$start <- anno$start - extend
##         anno$end <- anno$end + extend
##         if(class(eData) == "ExpData")
##             res <- splitByAnnotation(eData, anno = anno, ignoreStrand = TRUE, , what = c("location", what), expand = TRUE)[[1]]
##         else {
##             res <- lapply(eData, function(x) {
##                 subset(x, x[, "location"] >= anno$start &
##                           x[, "location"] <= anno$end)
##             })
##             what <- setdiff(colnames(res[[1]]), c("strand", "chr","location"))
##         }
##         tracks <- lapply(what, function(wh) {
##             counts.pos <- res[["1"]][, wh]
##             counts.neg <- res[["-1"]][, wh]
##             if(readLength > 1) {
##                 weights <- rep(1, readLength)
##                 pad <- rep(0, readLength - 1)
##                 res.pos <- convolve(c(pad, counts.pos), weights, type = "filter")
##                 res.neg <- rev(convolve(c(pad, rev(counts.neg)), weights, type = "filter"))
##                 values <- res.pos + res.neg
##             } else {
##                 values <- counts.pos + counts.neg
##             }
##             makeBaseTrack(base = res[["1"]][, "location"], value = values, dp = dp)
##         })
##         names(tracks) <- what
##         gdPlot(c(makeGenomeAxis(), tracks), minBase = anno$start, maxBase = anno$end)
##         return(invisible(list(tracks = tracks, results = res)))
##     }
## }



## resCoverage <- function(res, start, end, readLength) {
##     resnames <- setdiff(colnames(res), c("chr", "location", "strand"))
##     res[is.na(res)] <- 0
##     if(length(readLength) == 1)
##         readLength <- rep(readLength, length(resnames))
##     out <- lapply(seq(along = resnames), function(ii) {
##         ir <- IRanges(start = ifelse(res[, "strand"] == 1L,
##                       res[,"location"], NA),
##                       end = ifelse(res[, "strand"] == -1L,
##                       res[,"location"], NA),
##                       width = readLength[ii])
##         cov <- coverage(ir, weight = res[, resnames[ii]],
##                         shift = - start, width = end - start + 1)
##     })
##     names(out) <- resnames
##     out
## }
                  
    
## setClass("BaseTrackList", contains = c("gdObject", "ImplementsTrackOverlay"),
##          representation(basetracks = "list"),
##          prototype(basetracks = list(new("BaseTrack")),
##                    dp = DisplayPars(size = 5),
##                    trackOverlay = NULL)
## )

## setMethod("getGenomicRange", signature("BaseTrackList"), function(obj){
##     c(min(sapply(obj@basetracks, function(xx) min(xx@base))),
##       max(sapply(obj@basetracks, function(xx) max(xx@base))))
## })



## setMethod("drawGD", signature("BaseTrackList"), function(gdObject, minBase, maxBase, vpPosition) {
##     xlims <- sapply(gdObject@basetracks, function(xx) {
##         out <- getPar(xx, "xlim")
##         if(is.null(out))
##             out <- c(minBase, maxBase)
##         out
##     })
##     xlim <- c(min(xlims[,1]), max(xlims[,2]))

##     ylims <- sapply(gdObject@basetracks, function(xx) {
##         out <- getPar(xx, "ylim")
##         if(is.null(out))
##             out <- range(GenomeGraphs:::getBaseValue(xx), na.rm = TRUE, finite = TRUE)
##         out
##     })
##     ylim <- c(min(ylims[,1]), max(ylims[,2]))

##     drawAxis <- sapply(gdObject@basetracks, function(xx) {
##         out <- getPar(xx, "drawAxis")
##         if(is.null(out))
##             out <- TRUE
##     })
##     if(any(drawAxis))
##         drawAxis <- TRUE
##     else
##         drawAxis <- FALSE

##     if (drawAxis) {
##         pushViewport(dataViewport(xData = xlim, yData = ylim, extension = 0,
##                                   layout.pos.col = 1, layout.pos.row = vpPosition))
##     } else {
##         ## XXX: THESE ARE TOTAL HACKS FOR NOW.
##         if(is.na(ylim) || diff(ylim) == 0)
##             ylim <- c(0,.1)
##         if(any(!is.finite(ylim))) {
##             ylim <- c(0,.1)
##         }
##         pushViewport(dataViewport(xData = xlim, yscale = ylim, extension = 0, clip = TRUE,
##                                   layout.pos.col = 1, layout.pos.row = vpPosition))
##     }

##     ## here i probably want to vectorize these two.

##     lapply(gdObject@basetracks, function(xx) {
##         baseValue <- GenomeGraphs:::getBaseValue(xx)
##         pos <- GenomeGraphs:::getBase(xx)
##         lwd <- GenomeGraphs:::getPar(xx, "lwd")
##         lty <- GenomeGraphs:::getPar(xx, "lty")
##         pty <- GenomeGraphs:::getPar(xx, "type")
##         col <- if (length(color <- getPar(xx, "color")) == length(pos)) {
##             color
##         }
##         else {
##             rep(color, length(pos))[1:length(pos)]
##         }
##         whBase <- (pos > minBase & pos < maxBase)
##         col <- col[whBase]
##         baseValue <- baseValue[whBase]
##         pos <- pos[whBase]
##         if (sum(whBase) > 0) {
##         dP <- function() {
##             grid.points(pos, baseValue, default.units = "native", gp = gpar(col=col),
##                         size = unit(lwd, "char"), pch = 16)
##         }
##         dVL <- function() {
##             mapply(function(a, b, c) {
##                 grid.lines(x = c(a, a),
##                            y = c(min(ylim), b),
##                            default.units = "native", gp = gpar(col=c, lwd = unit(lwd, "char")))
##             }, pos, baseValue, col)
##         }
##         dCL <- function() {
##             grid.lines(x = pos, y = baseValue, default.units = "native", 
##                        gp = gpar(col = col, lwd = unit(lwd, "char")))
##         }

##         if (pty == "p") {
##             dP()
##         }
##         else if (pty == "h") {
##             dVL()
##         }
##         else if (pty == "l") {
##             dCL()
##         }
##         grid.yaxis()
##     }
##     })

##     GenomeGraphs:::.drawTrackOverlays(gdObject, minBase, maxBase)

##     popViewport()
## })
          

## makeCoveragePlotter2 <- function(expDataList, annoFactory, readLength = 1, dp = DisplayPars()) {
    
##     function(chr, start, end, overlays = NULL, title = NULL, ...) {
##         if (!is.null(annoFactory)) {
##           geneTracks <- lapply(c(-1, 0, 1), function(str) {
##               annoFactory(chr, start, end, strand = str)
##           })
##           names(geneTracks) <- c("+", "*", "-")
##         } else {
##           geneTracks <- NULL
##         }
##         if(is.null(geneTracks[["0"]])) {
##             annoTracks <- c("+" = geneTracks[["+"]],
##                             makeGenomeAxis(), 
##                             "-" = geneTracks[["-"]])
##         } else {
##             annoTracks <- c("+" = geneTracks[["+"]],
##                             makeGenomeAxis(), 
##                             "*" = geneTracks[["*"]],
##                             makeGenomeAxis(), 
##                             "-" = geneTracks[["-"]])
##         }
##         res <- lapply(expDataList, function(expData) {
##             res <- getRegion(expData = expData, chr = chr, start = start, end = end, strand = 0)
##         })
##         res <- Reduce(function(xx, yy) {
##             merge(xx, yy, by = getIndexColumns(expDataList[[1]]),
##                   all = TRUE)
##         }, res)
##         res[is.na(res)] <- 0
##         cov <- resCoverage(res, start = start, end = end, readLength = readLength)
        
##         res.1 <- res[res[,"strand"] == 1L,]
##         res.2 <- res[res[,"strand"] == -1L,]
##         cov.1 <- resCoverage(res.1, start = start, end = end, readLength = readLength)
##         cov.2 <- resCoverage(res.2, start = start, end = end, readLength = readLength)
                
##         dataTracks <- lapply(cov, function(cv) {
##             value <- rep(runValue(cv), each = 2)
##             base <- cumsum(runLength(cv))
##             base <- sort(c(base, base - runLength(cv) + 1))
##             makeBaseTrack(value = value, base = base + start - 1, dp = dp)
##         })

##         dataTracks2 <- mapply(function(cv1,cv2) {
##             value <- rep(runValue(cv1), each = 2)
##             base <- cumsum(runLength(cv1))
##             base <- sort(c(base, base - runLength(cv1) + 1))
##             bt1 <- makeBaseTrack(value = value, base = base + start - 1, dp = DisplayPars(color = "orange", type = "l"))
##             value <- rep(runValue(cv2), each = 2)
##             base <- cumsum(runLength(cv2))
##             base <- sort(c(base, base - runLength(cv2) + 1))
##             bt2 <- makeBaseTrack(value = value, base = base + start - 1, dp = DisplayPars(color = "red", type = "l"))
##             new("BaseTrackList", basetracks = list(bt1, bt2))
##         }, cov.1, cov.2)

##         allTracks <- Filter(function(x) !is.null(x), c(annoTracks, dataTracks))
##         gdPlot(c(allTracks, title), minBase = start, maxBase = end,
##                overlays = overlays, ...)
##         invisible(allTracks)
##     }
## }


##
## Useful for converting to pileup representation.
##
makeConvolver <- function(readLength = 32, direction = 1, weight = 1) {
  pad <- rep(0, readLength - 1)
  weights <- rep(weight, readLength)
  
  function(dta) {
    if (direction == 1) {
      dta <- c(pad, dta)
    } else {
      dta <- c(dta, pad)
    }
    convolve(dta, weights, type = "f")
  }
}
