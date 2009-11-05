require(Genominator)

options(verbose = TRUE)
N <- 1000000 # the number of observations. 
K <- 1000    # the number of annotation regions, not less than 10

annoData <- data.frame(chr = sample(1:16, size = K, replace = TRUE),
                       strand = sample(c(1, -1), size = K, replace = TRUE),
                       start = (st <- sample(1:1000, size = K, replace = TRUE)),
                       end = st + rpois(K, 75))

init <- function() {
  df <- data.frame(chr = sample(1:16, size = N, replace = TRUE),
                   location = sample(1:1000, size = N, replace = TRUE),
                   strand = sample(c(1L,-1L), size = N, replace = TRUE))
  head(df)
  eDataRaw <- importToExpData(df, filename = "my.db", tablename = "ex_tbl", 
                              overwrite = TRUE)
  eData <- aggregateExpData(eDataRaw, tablename = "counts_tbl", deleteOriginal = FALSE, 
                            overwrite = TRUE)
  
  return(eData)
}

timeQ <- function(q, prepend = "") {
  tm <- system.time({r <- dbGetQuery(getDB(eData), paste(prepend, q))})
  return(list(tm, r))
}

IDXS <- list(c("chr", "strand", "location"),
             c("chr", "location", "strand"),
             c("location", "strand", "chr"),
             c("location", "chr", "strand"),
             c("strand", "location", "chr"),
             c("strand", "chr", "location"))
names(IDXS) <- sapply(IDXS, paste, collapse = ", ")

L <- list(q0 = paste("SELECT __regions__.id, TOTAL(counts) FROM __regions__",
            "LEFT OUTER JOIN counts_tbl ON __regions__.chr = counts_tbl.chr",
            "AND counts_tbl.location BETWEEN __regions__.start AND __regions__.end__1",
            "AND (counts_tbl.strand = __regions__.strand OR __regions__.strand = 0",
            "OR counts_tbl.strand = 0) GROUP BY __regions__.id ORDER BY __regions__.id"),
          q1 = paste("SELECT __regions__.id, TOTAL(counts) FROM __regions__",
            "LEFT OUTER JOIN counts_tbl ON __regions__.chr = counts_tbl.chr",
            "AND (counts_tbl.strand = __regions__.strand OR __regions__.strand = 0 OR counts_tbl.strand = 0)",
            "AND counts_tbl.location BETWEEN __regions__.start AND __regions__.end__1",
            "GROUP BY __regions__.id ORDER BY __regions__.id"),
          q3 = paste("SELECT counts_tbl.* , __regions__.id FROM counts_tbl",
            "INNER JOIN __regions__ ON __regions__.chr = counts_tbl.chr",
            "AND counts_tbl.location BETWEEN __regions__.start AND __regions__.end__1",
            "AND (counts_tbl.strand = __regions__.strand OR __regions__.strand = 0 OR counts_tbl.strand = 0)",
            "ORDER BY __regions__.id, counts_tbl.chr, counts_tbl.location, counts_tbl.strand"),
          q4 = paste("SELECT counts_tbl.* , __regions__.id FROM counts_tbl",
            "INNER JOIN __regions__ ON __regions__.chr = counts_tbl.chr",
            "AND (counts_tbl.strand = __regions__.strand OR __regions__.strand = 0 OR counts_tbl.strand = 0)",
            "AND counts_tbl.location BETWEEN __regions__.start AND __regions__.end__1",
            "ORDER BY __regions__.id, counts_tbl.chr, counts_tbl.location, counts_tbl.strand"))

          
eData <- init()

S <- sapply(IDXS, function(idx) {
  sapply(L, function(q) {
    dbGetQuery(getDB(eData), "DROP INDEX counts_tblIDX")
    dbGetQuery(getDB(eData), paste("CREATE INDEX counts_tblIDX ON counts_tbl(", paste(idx, collapse = ","), ")"))
    timeQ(q)[[1]][3]
  })
})

