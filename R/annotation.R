
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

makeUIgenes <- function(annoData) {
    makeUIlocally <- function(chrAnno) {
        geneSplit <- split(chrAnno, chrAnno$gene_id)
        geneRanges <- do.call(IRangesList, lapply(geneSplit, function(gene) {
            reduce(IRanges(gene$start, gene$end))
        }))
        unionComplementToGene <- lapply(seq(along = geneRanges), function(i) {
            reduce(unlist(geneRanges[-i]))
        })
        consecutivePartsOfGene <- do.call(IRangesList, lapply(geneSplit, function(gene) {
            transcriptList <- lapply(split(gene, gene$transcript_id), function(exons) {
                reduce(IRanges(start = exons$start, end = exons$end))
            })
            transcriptIntersection <- Reduce(intersect, transcriptList)
        }))
        cleanedUIgene <- lapply(seq(along = consecutivePartsOfGene), function(i) {
            setdiff(consecutivePartsOfGene[[i]], unionComplementToGene[[i]])
        })
        reps <- sapply(cleanedUIgene, length)
        data.frame(chr = chrAnno$chr[1],
                   strand = rep(sapply(geneSplit, function(x) x$chr[1]), reps),
                   start = do.call(c, lapply(cleanedUIgene, start)),
                   end = do.call(c, lapply(cleanedUIgene, end)),
                   gene_id = rep(names(geneSplit), reps))
    }
    chrSplit <- split(annoData, annoData$chr)
    do.call(rbind, lapply(chrSplit, makeUIlocally))
}

## system.time(tmp <- makeUIgenes(chr1Anno))
