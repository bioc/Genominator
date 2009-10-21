##
## Test the hasOverlap and isContained functions. 
##

test1 <- data.frame(chr = c(1,1,1), strand = c(1,1,1), start = c(2,2,2), end = c(5,5,5))
a <- isContained(test1, test1)
stopifnot(all(sapply(a, function(b) all(b == c(1,2,3)))))

hasOverlap(test1, test1)
test2 <- data.frame(chr = 1, strand = 1, start = c(3,3,1,7,1), end = c(4,8,4,10,10))
a <- isContained(test2, test1, ignoreStrand = FALSE)
stopifnot(all(a[[1]] == c(1,2,3)))

a <- hasOverlap(test2, test1, ignoreStrand = FALSE)
stopifnot(all(sapply(a[c(1,2,3,5)], function(b) all(b == c(1,2,3)))))

test2 <- data.frame(chr = 1, strand = -1, start = c(3,3,1,7,1), end = c(4,8,4,10,10))
a <- isContained(test2, test1, ignoreStrand = TRUE)
stopifnot(all(a[[1]] == c(1,2,3)))

a <- hasOverlap(test2, test1, ignoreStrand = TRUE)
stopifnot(all(sapply(a[c(1,2,3,5)], function(b) all(b == c(1,2,3)))))
