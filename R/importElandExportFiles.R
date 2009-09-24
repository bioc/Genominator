## Here is an example call: 
##
## importElandExportFiles(inFiles = fl, database = DBH, tablename = "exp_2",
##                        overwrite = TRUE, verbose = TRUE, preprocess = TRUE,
##                        chrMap = c(paste("c", c(1:22, "X", "Y"), ".fa", sep = ""), "splice_sites-31.fa"),
##                        nrows = -1)

## inFiles: is a named vector of files to read and import where the
## names correspond the columns in the final joined table and thus
## they are quite relevant.

## chrMap: is an ordered vector of names which will correspond to
## chromosome 1,2, ... it is important that all of the chromosomes be
## included, otherwise NAs will arise. The name of the chromosome is
## generally the fasta file which it mapped to.

## preprocess: essentially does the initial preprocessing the export
## files removing columns, sometimes it is useful for debugging to do
## this once and then set the argument to FALSE.

## nrows: says how many rows of the resulting data.frame from each
## export file should be included. This is useful for debugging as
## reading only 100000 rows makes a huge difference from 5 million.

importElandExportFiles <- function(inFiles, database, tablename, chrMap,
                                   purityFilter = TRUE,
                                   matchPerfectly = NULL,
                                   overwrite = FALSE,
                                   nrows = -1, preprocess = TRUE, verbose = FALSE) {
    if (is.null(names(inFiles))) {
        stop("inFiles must have names which corresponsd to eventual column names.")
    }
    if (missing(database) || missing(tablename)) {
        stop("Must specify both tablename and database.")
    }
    
    cleanFiles <- paste(inFiles, "-clean", sep = "")
    
    if (preprocess) {
        sapply(inFiles, function(nm) {
            if (purityFilter) 
                pf <- "grep -v N$ |"
            else
                pf <- ""

            if (!is.null(matchPerfectly)) {
                matchPerfectly <- paste("| grep ", matchPerfectly, "$", sep = "")
            }
            matchPerfectly <- paste(matchPerfectly, "| cut -f 1-3 -d\" \"")

            tmp <- sprintf(paste("cat %s | cut -f 11,13,14,15,22 |", pf, "grep -v NM |",
                                 "grep -v RM | grep -v QC | cut -f 1-4 | sed 's/\\\t/ /g'",
                                 matchPerfectly, "> %s-clean"), nm, nm)
            if (verbose) cat("processing sed:\n", tmp, "\n")
            system(tmp)
        })
    }

  eDatas <- mapply(function(file, name) {
      if (verbose)
          print(paste("processing:", file, "with name:", name))
      tbl <- read.table(file, stringsAsFactors = FALSE, header = FALSE, nrows = nrows)
      colnames(tbl) <- c("chr", "location", "strand")
      
      tbl$chr <- match(tbl$chr, chrMap)
      tbl$strand <- ifelse(tbl$strand == "F", 1, -1)
      
      nnn <- (is.na(tbl$chr) | is.na(tbl$strand) | is.na(tbl$location))
      if (any(nnn)) {
          warning(paste(sum(nnn), "NAs found after doing chromosome/strand mapping, dropping."))
          print(tbl[nnn,])
          tbl <- tbl[!nnn, ] 
      }
      ed <- importToExpData(tbl, database, tablename = name, overwrite = TRUE, verbose = verbose)

      if (verbose) {
          print(head(ed))
      }
      
      ad <- aggregateExpData(ed, deleteOriginal = TRUE, overwrite = TRUE, verbose = verbose, colname = name)
      if (verbose) {
          print(head(ad))
      }
    return(ad)
  }, cleanFiles, names(inFiles))
  
  jd <- joinExpData(eDatas, tablename = tablename, verbose = verbose, overwrite = overwrite, deleteOriginals = TRUE)
  
  if (preprocess) {
      if (any(wh <- (!sapply(cleanFiles, file.remove))))
          warning(paste("Unable to delete some intermediate files.", paste(cleanFiles[wh], collapse = ", ")))
  }
  return(jd)
}
