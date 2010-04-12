validAnnotation <- function(annoData) {
    if(!is.data.frame(annoData))
        stop("An annotation object should be a data.frame")
    if(length(intersect(c("chr", "strand", "start", "end"), names(annoData))) != 4)
        stop("An annotation object should have chr, strand, start, stop information")
    if(!is.numeric(annoData$chr) || !is.numeric(annoData$strand) ||
       !is.numeric(annoData$start) || !is.numeric(annoData$end))
        stop("An annotation object should have numeric genomic locations")
    if(any(is.na(annoData$chr)) || any(is.na(annoData$strand)) ||
       any(is.na(annoData$start)) || any(is.na(annoData$end)))
        stop("An annotation object should have non-missing genomic locations")
    if(!all(is.element(annoData$strand, c(-1,0,1))))
        stop("An annotation object should have strand values in {-1,0,1}")
    if(!(all(annoData$start <= annoData$end) || 
         (annoData$start[annoData$strand != -1] <= annoData$end[annoData$strand != -1] &&
          annoData$start[annoData$strand == -1] >= annoData$end[annoData$strand == -1])))
        stop("An annotation object should have either all start <= end or only start > end if strand == -1")
    return(TRUE)
}

makeGeneRepresentation <- function(annoData, type = c("UIgene", "Ugene", "ROCE", "background"),
                                   gene.id = "ensembl_gene_id", transcript.id = "ensembl_transcript_id",
                                   bind.columns, ignoreStrand = TRUE, verbose = getOption("verbose")) {
    ## Want to augument to compute ROCEs, have annotation that only is part of the subtraction phase.
    findRelevantGeneOverlap <- function(geneRanges) {
        ## This function takes a names IRangesList and find the overlap between any of the components and
        ## the rest of the IRanges
        reduced <- reduce(geneRanges)
        allReduced <- unlist(reduced)
        overlaps <- matchMatrix(findOverlaps(allReduced, ignoreSelf = TRUE))
        geneNamesIdx <- seq(along = allReduced)
        names(geneNamesIdx) <- rep(names(reduced), times = sapply(reduced, length))
        relevantGeneOverlap <- do.call(IRangesList, lapply(names(reduced), function(nam) {
            geneIdx <- geneNamesIdx[names(geneNamesIdx) == nam]
            relevantOverlaps <- overlaps[overlaps[, "query"] %in% geneIdx &
                                         !overlaps[, "subject"] %in% geneIdx, "subject"]
            out <- reduce(allReduced[relevantOverlaps])
            out
        }))
        names(relevantGeneOverlap) <- names(reduced)
        relevantGeneOverlap
    }
    
    formFinalDataFrame <- function(geneSplit, rangesList, bindColumns){
        reps <- sapply(rangesList, length)
        newSplit <- geneSplit[names(rangesList)]
        getColumn <- function(nam){
            rep(sapply(newSplit, function(x) x[1, nam]), times = reps)
        }
                
        df <- data.frame(chr = getColumn("chr"),
                         start = unlist(lapply(rangesList, function(x) start(x))),
                         end = unlist(lapply(rangesList, function(x) end(x))),
                         strand = getColumn("strand"),
                         rep(names(rangesList), times = reps),
                         stringsAsFactors = FALSE)
        names(df)[5] <- gene.id
        if(length(bindColumns) > 0) {
            names(bindColumns) <- bindColumns
            bindList <- lapply(bindColumns, getColumn)
            df <- cbind(df, bindList)
        }
        rownames(df) <- NULL
        df
    }
        
    computeUorUI <- function(chrAnno) {
        if(verbose) cat("computing on chromosome", chrAnno$chr[1], "\n")
        geneSplit <- split(chrAnno, chrAnno[, gene.id])
        geneRanges <- do.call(IRangesList, lapply(geneSplit, function(gene) {
            IRanges(gene$start, gene$end)
        }))
        reducedGeneRanges <- reduce(geneRanges)
        relevantGeneOverlap <- findRelevantGeneOverlap(reducedGeneRanges)
        geneModel <- switch(type, "UIgene" = {
            intersectGeneRanges <- do.call(IRangesList, lapply(geneSplit, function(gene) {
                transcriptList <- lapply(split(gene, gene[, transcript.id]), function(exons) {
                    reduce(IRanges(start = exons$start, end = exons$end))
                })
                Reduce(intersect, transcriptList)
            }))
            intersectGeneRanges
        }, "Ugene" = {
            reducedGeneRanges
        })

        nameList <- names(geneModel)
        names(nameList) <- names(geneModel)
        geneModel <- do.call(IRangesList, lapply(nameList, function(nam) {
            geneOverlap <- relevantGeneOverlap[[nam]]
            model <- geneModel[[nam]]
            if(length(geneOverlap) > 0) {
                model <- setdiff(model, geneOverlap)
            }
            model
        }))
        df <- formFinalDataFrame(geneSplit, geneModel, bindColumns = bind.columns)
        df
    }
    
    computeROCE <- function(chrAnno) {
        if(verbose) cat("computing on chromosome", chrAnno$chr[1], "\n")
        geneSplit <- split(chrAnno, chrAnno[, gene.id])
        geneRanges <- do.call(IRangesList, lapply(geneSplit, function(gene) {
            IRanges(gene$start, end = gene$end)
        }))
        relevantGeneOverlap <- findRelevantGeneOverlap(geneRanges)
        disjointGeneRanges <- disjoin(geneRanges)
        
        nameList <- names(disjointGeneRanges)
        names(nameList) <- names(disjointGeneRanges)
        disjointGeneRanges <- do.call(IRangesList, lapply(nameList, function(nam) {
            disjRange <- disjointGeneRanges[[nam]]
            geneOverlap <- relevantGeneOverlap[[nam]]
            if(length(geneOverlap) > 0) {
                disjRange <- do.call(c, lapply(seq(along = disjRange), function(ii) {
                    setdiff(disjRange[ii], geneOverlap)
                }))
            }
            disjRange
        }))
        df <- formFinalDataFrame(geneSplit, disjointGeneRanges, bindColumns = bind.columns)
        df
    }
    
    computeBackground <- function(chrAnno) {
        if(verbose) cat("computing on chromosome", chrAnno$chr[1], "\n")
        allRanges <- reduce(IRanges(start = chrAnno$start, end = chrAnno$end))
        background <- setdiff(IRanges(start = 1, end = max(chrAnno$end)), allRanges)
        if(max(chrAnno$strand) > 0.5) {
            strand <- 1L
        } else {
            strand <- -1L
        }
        df <- data.frame(chr = chrAnno$chr[1],
                         strand = strand,
                         start = start(background),
                         end = end(background), stringsAsFactors = FALSE)
        rownames(df) <- NULL
        df
    }

    if(!require(IRanges))
        stop("This functionality requirees 'IRanges'")
    if(class(annoData) != "data.frame")
        stop("Argument 'annoData' needs to be a data.frame")
    type <- match.arg(type)
    if(!missing(bind.columns)) {
        bindOK <- all(bind.columns %in% names(annoData))
    } else {
        bindOK <- TRUE
        bind.columns <- character(0)
    }
    
    switch(type, "UIgene" = {
        stopifnot(gene.id %in% names(annoData),
                  transcript.id %in% names(annoData),
                  bindOK)
        representation <- computeUorUI
    }, "Ugene" = {
        if(is.null(gene.id)) {
            gene.id <- "FAKE_GENE_ID"
            annoData$"FAKE_GENE_ID" <- "A"
        }
        stopifnot(gene.id %in% names(annoData),
                  bindOK)
        representation <- computeUorUI
    }, "background" = {
        representation <- computeBackground
    }, "ROCE" = {
        stopifnot(gene.id %in% names(annoData),
                  bindOK)
        representation <- computeROCE
    })

    chrSplit <- split(annoData, annoData$chr)

    if(!ignoreStrand) {
        chrSplit <- do.call(c, lapply(chrSplit, function(chrData) {
            list(subset(chrData, strand %in% c(-1L, 0L)),
                 subset(chrData, strand %in% c(1L, 0L)))
        }))
    }

    df <- do.call(rbind, lapply(chrSplit, representation))
    rownames(df) <- NULL
    df
}

