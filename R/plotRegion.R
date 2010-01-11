toAnnoData <- function(annoTrack) {
    cbind(chr = rep(annoTrack@chr, nrow(annoTrack@regions)), start = annoTrack@regions[, "start"],
          end = annoTrack@regions[, "end"], strand = rep(annoTrack@strand, nrow(annoTrack@regions)),
          annoTrack@regions[,-match(c("start", "end"), colnames(annoTrack@regions))])
}

makeAnnoFactory.Biomart <- function(biomart, dp = DisplayPars(),
                                    chrFunction = function(x) x,
                                    strandFunction = function(x) x) {
    if (!require(GenomeGraphs))
        stop("GenomeGraphs must be installed to use this functionality.")

    cChrFunction <- function(y) as.character(chrFunction(y))
    cStrandFunction <- function(y) as.character(strandFunction(y))
    
    function(chr, start, end, strand) { 
        geneRegionBiomart(chr, start, end, strand, biomart, dp, chrFunction = cChrFunction,
                          strandFunction = cStrandFunction)
    }
}

makeAnnoFactory.AnnoData <- function(annoData,
                                     dp = DisplayPars(),
                                     featureColumnName = "feature",
                                     groupColumnName = NULL,
                                     idColumnName = NULL,
                                     expandAround = 1000) {
    if (!require(GenomeGraphs))
        stop("GenomeGraphs must be installed to use this functionality.")
    
    function(chr, start, end, strand) {
        mat <- annoData[annoData[,"start"] >= start - expandAround & annoData[,"end"] <= end + expandAround &
                        annoData[,"chr"] == chr & annoData[,"strand"] == strand , ]
        
        if (nrow(mat) == 0) {
            return(NULL)
        }
        
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
        
        regions <- data.frame(start = mat[,"start"], end = mat[, "end"],
                              feature = mat[, featureColumnName], group = groups, ID = IDs)

        makeAnnotationTrack(chr = chr, strand = strand, regions =
                            regions, dp = dp)
    }
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

makeCoveragePlotter <- function(eData, what, readLength = 1, dp = DisplayPars()) {
    ## use both anno as char and data.frame and also accept chr, start, end
    makeTracks <- function(anno, extend = 0) {
        if(nrow(anno) != 1) stop("")
        anno$start <- anno$start - extend
        anno$end <- anno$end + extend
        if(class(eData) == "ExpData")
            res <- splitByAnnotation(eData, anno = anno, ignoreStrand = TRUE, , what = c("location", what), expand = TRUE)[[1]]
        else {
            res <- lapply(eData, function(x) {
                subset(x, x[, "location"] >= anno$start &
                          x[, "location"] <= anno$end)
            })
            what <- setdiff(colnames(res[[1]]), c("strand", "chr","location"))
        }
        tracks <- lapply(what, function(wh) {
            counts.pos <- res[["1"]][, wh]
            counts.neg <- res[["-1"]][, wh]
            if(readLength > 1) {
                weights <- rep(1, readLength)
                pad <- rep(0, readLength - 1)
                res.pos <- convolve(c(pad, counts.pos), weights, type = "filter")
                res.neg <- rev(convolve(c(pad, rev(counts.neg)), weights, type = "filter"))
                values <- res.pos + res.neg
            } else {
                values <- counts.pos + counts.neg
            }
            makeBaseTrack(base = res[["1"]][, "location"], value = values, dp = dp)
        })
        names(tracks) <- what
        gdPlot(c(makeGenomeAxis(), tracks), minBase = anno$start, maxBase = anno$end)
        return(invisible(list(tracks = tracks, results = res)))
    }
}


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
