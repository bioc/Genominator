##
## Here I am testing whether or not various functions do what
## I expect. 
##
require(Genominator)

N <- 1000
M <- 100

rawData <- data.frame(chr = sample(1:2, size = N, replace = TRUE),
                      location = floor(runif(N, 1, 1e3)),
                      strand = sample(c(-1,1), size = N, replace = TRUE))
eData <- importToExpData(rawData, tablename = "test", filename = "test.db",
                         overwrite = TRUE, verbose = TRUE)
eData <- aggregateExpData(eData, overwrite = TRUE, verbose = TRUE)

r1 <- data.frame(chr = 1, start = 10, end = 100, strand = 1)
r2 <- r1; r2$strand == 0

##
## testing splitByAnnotation
##
a1 <- rawData[rawData$chr == r1$chr & rawData$location >= r1$start & rawData$location <= r1$end &
              rawData$strand == r1$strand, ]
a2 <- splitByAnnotation(eData, r1)[[1]]
if (sum(a2[,"counts"]) != nrow(a1))
  stop("Counts should add up!")

a1 <- rawData[rawData$chr == r1$chr & rawData$location >= r1$start & rawData$location <= r1$end, ]
a2 <- splitByAnnotation(eData, r1, ignoreStrand = TRUE)[[1]]
if (sum(a2[,"counts"]) != nrow(a1))
  stop("Counts should add up!")

r1$strand <- 0
a1 <- rawData[rawData$chr == r1$chr & rawData$location >= r1$start & rawData$location <= r1$end, ]
a2 <- splitByAnnotation(eData, r1)[[1]]
if (sum(a2[,"counts"]) != nrow(a1))
  stop("Counts should add up!")

a2 <- splitByAnnotation(eData, r1, what = "counts")[[1]]
if (sum(a2) != nrow(a1))
  stop("Counts should ad up!")

a2 <- splitByAnnotation(eData, r1, what = "counts", expand = TRUE, ignoreStrand = TRUE,
                        addOverStrands = TRUE)[[1]]
if (sum(a2) != nrow(a1))
  stop("Counts should ad up!")

if(length(a2) != (r1$end - r1$start + 1))
  stop("Expand should have same length of returned object.")

a2 <- splitByAnnotation(eData, r1, what = "counts", expand = TRUE, ignoreStrand = TRUE,
                        addOverStrands = FALSE)[[1]]

a1m <- rawData[rawData$chr == r1$chr & rawData$location >= r1$start & rawData$location <= r1$end &
              rawData$strand == -1, ]
a1p <- rawData[rawData$chr == r1$chr & rawData$location >= r1$start & rawData$location <= r1$end &
               rawData$strand == 1, ]

if (nrow(a1m) != sum(a2[["-1"]]) ||
    nrow(a1p) != sum(a2[["1"]]))
  stop("Strands problem.")


if (summarizeByAnnotation(eData, r1) !=
  sum(splitByAnnotation(eData, r1, what = "counts", expand = TRUE, ignoreStrand = TRUE,
                    addOverStrands = TRUE)[[1]]))
  stop("Summarize and split disagreement.")



