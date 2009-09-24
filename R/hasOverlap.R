isContained <- function(anno1, anno2, extend = 0, ignoreStrand = FALSE) {
    if(!is.data.frame(anno1) || !is.data.frame(anno2))
        stop("anno1 and anno2 needs to be data.frame")
    s1 <- anno1[,.ANNO.COLS]
    s2 <- anno2[,.ANNO.COLS]
    s2$start <- s2$start - extend
    s2$end <- s2$end + extend
    outerComp <- (outer(s1$chr, s2$chr, "==") &
                  outer(s1$start, s2$start, ">=") &
                  outer(s1$end, s2$end, "<="))
    if(!ignoreStrand)
        outerComp <- outerComp & outer(s1$strand, s2$strand, "==")
    lapply(1:nrow(s1), function(i) {
        which(outerComp[i,])
    })
}

hasOverlap <- function (anno1, anno2, extend = 0, ignoreStrand = FALSE){ 
    if(!is.data.frame(anno1) || !is.data.frame(anno2))
        stop("anno1 and anno2 needs to be data.frame")
    s1 <- anno1[,.ANNO.COLS]
    s2 <- anno2[,.ANNO.COLS]
    s2$start <- s2$start - extend
    s2$end <- s2$end + extend
    sameChr <- outer(s1$chr, s2$chr, "==")
    if(!ignoreStrand)
        sameChr <- sameChr & outer(s1$strand, s2$strand, "==")
    isContained <- (outer(s1$start, s2$start, ">=") &
                    outer(s1$end, s2$end, "<="))
    isPartial <- ( (outer(s1$start, s2$start, "<=") &
                    outer(s1$end, s2$start, ">=")) |
                   (outer(s1$start, s2$end, "<=") &
                    outer(s1$end, s2$end, ">=")) )
    outerComp <- sameChr & (isContained | isPartial)
    lapply(1:nrow(s1), function(i) {
        which(outerComp[i,])
    })
}

## hasOverlap <- function(segment, anno, ignoreStrand = FALSE, pad = 0){
##   if(nrow(segment) > 1)
##     stop("too many segments")
##   if(ignoreStrand)
##     checkStrand <- TRUE
##   else
##     checkStrand <- (segment$strand == anno$strand)
##   partial <- which(segment$chr == anno$chr & checkStrand &  (
##                       (segment$start >= anno$start - pad & segment$start <= anno$end + pad) |
##                       (segment$end >= anno$start - pad & segment$end <= anno$end + pad)
##                       ))
##   within <- which(segment$chr == anno$chr & checkStrand &
##                   segment$start <= anno$start - pad & segment$end >= anno$end + pad)
##   overlap <- sort(unique(c(partial, within)))
##   return(overlap)
## }

