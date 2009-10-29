
##-- These are the columns which are necessary.
.ANNO.COLS <- c("chr", "start", "end", "strand")
.REGION.TABLE.NAME <- '__regions__'
.REGION.TABLE.TMP.NAME <- "__tmp_regions__"

.makeIndexName <- function(tablename) {
    paste(tablename, "IDX", sep = "")
}

importFromAlignedReads <- function(alignedReads, chrMap, filename, tablename,
                                   overwrite = TRUE, deleteIntermediates = TRUE,
                                   verbose = FALSE, ...) {
    if (!require(ShortRead))
        stop("ShortRead package must be installed.")
    if (!all(sapply(alignedReads, class) == "AlignedRead"))
        stop("alignedReads must be a list of AlignedRead objects.")

    laneNames <- names(alignedReads)
    
    if (is.null(laneNames) || any(laneNames == "") || length(unique(laneNames)) != length(laneNames))
        stop("alignedReads must be a named list with unique names.")

    expDatas <- mapply(function(name, aln) {
        loc <- position(aln)
        str <- c(-1,1)[match(strand(aln), c("-", "+"))] 
        chr <- match(chromosome(aln), chrMap)
        ed <- importToExpData(data.frame(chr = chr, location = loc, strand = str), filename = filename,
                              tablename = name, overwrite = TRUE, verbose = verbose)
        aggregateExpData(ed, colname = name, verbose = verbose, overwrite = TRUE)
    }, laneNames, alignedReads)
    joinExpData(expDatas, tablename = tablename, overwrite = overwrite, verbose = verbose,
                deleteOriginals = deleteIntermediates)
}

##
## Imports data from a data.frame.
##
importToExpData <- function(df, filename, tablename, overwrite = FALSE, verbose = FALSE,
                            columns = NULL) {
  db <- dbConnect(dbDriver("SQLite"), filename)
  allCols <- colnames(df)
  
  if (missing(columns))
    COLS <- attr(getClass("ExpData")@prototype, "indexColumns")
  else
    COLS <- columns
  
  if (!all(COLS %in% allCols))
    stop(paste("data.frame must have:", paste(COLS, collapse = ", "), "columns"))

  ## order them.
  df <- df[, c(COLS, setdiff(allCols, COLS))]

  ## now make the required columns integers.
  for (i in length(COLS)) {
    if(!is.numeric(df[,i]))
      stop("df needs only numeric columns")
    if (!is.integer(df[,i]))
      df[,i] <- as.integer(df[,i])
  }
  
  if (ncol(df) > length(COLS)) {
    for (i in (1 + length(COLS)):ncol(df)) {
      if (!extends(class(df[,i]), "numeric"))
        stop("Currently only numeric columns are supported!")
    }
  }

  ## now order them while you have the table in memory because
  ## constructing the index is infinitely faster this way.
  df <- df[order(df[,COLS[1]], df[,COLS[2]], df[,COLS[3]]), ]
  
  .timeAndPrint( { if (!dbWriteTable(db, tablename, df, row.names = FALSE, overwrite = overwrite))
                     stop("Unable to write the table, delete table or specify overwrite.")}, "Writing table",
                print = verbose)

  idxName <- .makeIndexName(tablename)
  q <- paste("CREATE INDEX ",  idxName, " ON ", tablename, "(",
             paste(COLS, collapse = ", "),");", sep = "")

  if (overwrite) {
      tryCatch(dbGetQuery(db, paste("DROP INDEX", idxName)), error = function(x) {})
  }

  .timeAndPrint(tryCatch(x <- dbGetQuery(db, q), error = print, q), "Creating index",
                print = verbose)
  
  ## clean up -- we will reconnect in the line below.
  dbDisconnect(db)
  
  return(ExpData(filename, tablename, indexColumns = COLS, mode = 'w'))
}

aggregateExpData <- function(expData, by = getIndexColumns(expData), tablename = NULL, deleteOriginal = FALSE,
                             overwrite = FALSE, verbose = FALSE, colname = "counts",
                             aggregator = paste("count(", by[1], ")", sep = ""))
{
    moveTable <- FALSE

    if (is.null(tablename)) {
        tablename <- sprintf("__tmp_%d", floor(runif(1, 1000, 9999)))
        deleteOriginal <- TRUE
        moveTable <- TRUE
    }

    if (overwrite) {
        tryCatch(dbGetQuery(getDB(expData), paste("DROP table", tablename)),
                 error = function(x) {})
    }

    cols <- paste(paste(c(by, colname), "INTEGER"), collapse = ",")
    .timeAndPrint(dbGetQuery(getDB(expData), sprintf("CREATE TABLE %s (%s)", tablename, cols)),
                  paste("Creating table:", tablename), print = verbose)

    statement <- sprintf("INSERT INTO %s SELECT %s FROM %s GROUP BY %s",
                         tablename,
                         paste(c(by, aggregator), collapse = ","),
                         getTablename(expData),
                         paste(by, collapse = ","))
    .timeAndPrint(dbGetQuery(getDB(expData), statement),
                  "inserting", print = verbose)

    if (deleteOriginal) {
        .timeAndPrint(dbGetQuery(getDB(expData), paste("DROP TABLE", getTablename(expData))),
                      "droping original table", print = verbose)
    }

    if (moveTable) {
        .timeAndPrint(dbGetQuery(getDB(expData), sprintf("ALTER TABLE %s RENAME TO %s",
                                                         tablename, getTablename(expData))),
                      "renaming table", print = verbose)

        tablename <- getTablename(expData)
    }

    ## now create the index on the correct tablename.
    statement <-  sprintf("CREATE INDEX %s ON %s (%s);", .makeIndexName(tablename),
                          tablename, paste(by, collapse = ","))
    .timeAndPrint(dbGetQuery(getDB(expData), statement), "creating index", print = verbose)

    ## return a new expData.
    return(ExpData(getDBName(expData), tablename, indexColumns = by, mode = 'w'))
}

collapseExpData <- function(expData, tablename = NULL, what = getColnames(expData, all = FALSE),
                            groups = "COL", collapse = c("sum", "avg", "weighted.avg"),
                            overwrite = FALSE, deleteOriginal = FALSE, verbose = FALSE) {
    ## XXX: Code Duplication!
    moveTable <- FALSE
    if (is.null(tablename)) {
        tablename <- sprintf("__tmp_%d", floor(runif(1, 1000, 9999)))
        deleteOriginal <- TRUE
        moveTable <- TRUE
    }
    if (overwrite) {
        tryCatch(dbGetQuery(getDB(expData), paste("DROP table", tablename)),
                 error = function(x) {})
    }

    if (length(groups) == 1)
        groups <- rep(groups, length(what))

    newCols <- groups[!duplicated(groups)]

    if (length(groups) != length(what))
        stop("Groups and what must match up in length.")

    collapse <- match.arg(collapse)

    sel <- tapply(what, groups, function(cols) {
        groupTotals <- paste("TOTAL(", cols, ")", sep = "")
        switch(collapse,
               sum = {
                   types <- getSchema(expData)[cols]
                   if(all(toupper(types) == "INTEGER")) {
                       groupTotals <- paste("CAST(", groupTotals, "AS INTEGER)", sep = "")
                   }
                   paste(groupTotals, collapse = "+") },
               avg = sprintf("(%s)/%s", paste(groupTotals, collapse = "+"), length(groupTotals)),
               weighted.avg = {
                   totals <- summarizeExpData(expData, what = cols)
                   totals <- totals/sum(totals)
                   paste(paste(groupTotals, "*", totals), collapse = "+")
               })
    })[newCols]

    ## Unfortunately, the TOTAL function in SQLite always returns a REAL,
    ## even when applied to INTEGERS. Operating on REALS are "slow" in my
    ## experience, so if we use "collapse = sum" I want to possibly insert
    ## a cast statement

    ## XXX : Here, i have enforced something about index columns that I don't probably want
    ##       to enforce, INTEGERness
    statement <- sprintf("CREATE TABLE %s (%s)", tablename, paste(c(paste(getIndexColumns(expData), "INTEGER"), newCols),
                                                                  collapse = ", "))
    .timeAndPrint(dbGetQuery(getDB(expData), statement), "creating table", print = verbose, statement)

    statement <- sprintf("INSERT INTO %s SELECT %s FROM %s GROUP BY %s",
                         tablename,
                         paste(c(getIndexColumns(expData), sel), collapse = ", "), getTablename(expData),
                         paste(getIndexColumns(expData), collapse = ","), paste(getIndexColumns(expData), collapse = ","))
    .timeAndPrint(dbGetQuery(getDB(expData), statement), "inserting data", print = verbose, statement)

    ## XXX: Code duplication
    if (deleteOriginal) {
        .timeAndPrint(dbGetQuery(getDB(expData), paste("DROP TABLE", getTablename(expData))),
                      "droping original table", print = verbose)
    }
    if (moveTable) {
        .timeAndPrint(dbGetQuery(getDB(expData), sprintf("ALTER TABLE %s RENAME TO %s",
                                                         tablename, getTablename(expData))),
                      "renaming table", print = verbose)
        tablename <- getTablename(expData)
    }

    ## now create the index on the correct tablename.
    statement <-  sprintf("CREATE INDEX %s ON %s (%s);", .makeIndexName(tablename),
                          tablename, paste(getIndexColumns(expData), collapse = ","))
    .timeAndPrint(dbGetQuery(getDB(expData), statement), "creating index", print = verbose)

    ## return a new expData.
    return(ExpData(getDBName(expData), tablename, mode = 'w'))
}

##-- Infer the column types. pretty lame.
dbListFieldsAndTypes <- function(con, tablename) {
  TYPE.MAP <- c("INTEGER", "REAL", "VARCHAR")
  names(TYPE.MAP) <- c("integer", "numeric", "character")

  q <- paste("SELECT * FROM", tablename, "LIMIT 1;")
  df <- dbGetQuery(con, q)
  res <- character(ncol(df))

  for (i in 1:ncol(df)) {
    res[i] <- TYPE.MAP[class(df[,i])]
  }
  names(res) <- colnames(df)

  return(res)
}

joinExpData <- function(expDataList, fields = NULL, tablename = "aggtable",
                        overwrite = TRUE, deleteOriginals = FALSE, verbose = FALSE){

  ## Does this work if fields is empty?
  ## Arguments:
  ## expDataList : list of expData objects
  ## fields : a list, preferably named as the tables in the expDataList,
  ##  each component is a named vector where the names of the vector are
  ##  the original column names of the table and the values of the vector
  ##  are the new columns names.
  ##  Example fields = list(c("a" = "b"), c("a" = "c"))
  ## tablename : the new tablename
  ## overwrite : overwrite the new table?
  ## deleteOriginals : delete all tables in the expDataList
  ## verbose : timing
  if(class(expDataList) != "list" || length(expDataList) < 2)
    stop("argument 'expDataList' must be a list of at least 2 components")
  if(!length(unique((sapply(expDataList, getDBName)))))
    stop("All expData must point to the same database!")
  if(!is.null(names(fields)) && !(setequal(sapply(expDataList, getTablename), names(fields))))
    stop("fields must be named appropriately: as the tables in expDataList")
  if(!is.null(names(fields)))
    fields <- fields[sapply(expDataList, getTablename)]

  ## enforce that the indices are the same.
  if (!all(length(Reduce(intersect, lapply(expDataList, getIndexColumns))) ==
           sapply(expDataList, function(a) length(getIndexColumns(a))))) {
    stop("Index columns must be the same for expData joining.")
  }
  .COLS <- getIndexColumns(expDataList[[1]])
    
  efields <- lapply(expDataList, function(x) {
    setdiff(getColnames(x), .COLS)
  })
  if(is.null(fields))
    nfields <- efields
  else {
    nfields <- mapply(function(x,y) {
      if(is.null(y))
        return(x)
      if(length(y) > 0 && is.null(names(y)))
        stop("'fields' must be properly named")
      m <- match(x, names(y))
      x[!is.na(m)] <- y[na.omit(m)]
      return(x)
    }, efields, fields, SIMPLIFY = FALSE)
  }
  if(any(duplicated(unlist(nfields))))
    stop("The merged table will have duplicated column names, use 'fields' to correct.")

  ## currently I am not using types for anything, probably wrong...
  types <- mapply(function(x,y) {
    getSchema(x)[y]
  }, expDataList, efields, SIMPLIFY = FALSE)

  nExpData <- length(expDataList)
  expdatatables <- sapply(expDataList, getTablename)
  db <- getDB(expDataList[[1]])

  if(overwrite) {
    dbGetQuery(db, paste("DROP TABLE IF EXISTS", tablename, ";"))
  }

  ## We now make nExpData __TEMP__ tables and we let the last one be
  ## the new table.
  temptables <- c(paste("__TEMP__", 1:nExpData, sep = ""), tablename)
  dbGetQuery(db, sprintf("CREATE TABLE %s (%s);", temptables[1],
                         paste(.COLS, "INTEGER", collapse = ", ")))
  COMMA.COLS <- paste(.COLS, collapse = ", ")
  select <- paste("SELECT", COMMA.COLS, "FROM", expdatatables, collapse = " UNION ")
  statement <- sprintf("INSERT INTO %s (%s) %s ORDER BY %s;",
                       temptables[1], COMMA.COLS, select, COMMA.COLS)
  .timeAndPrint(dbGetQuery(db, statement), "Creating union", print = verbose)

  ## Now we start a for loop where we add one expData at a time.
  ## This is a loop of sequential outer left joins.
  ## There is no need for a full outer join because of the way we have
  ## constructed __TEMP__1
  currentCols <- .COLS
  currentTypes <- rep("INTEGER", length(.COLS))
  joinstatement <- paste(paste("a", .COLS, sep = "."), "=",
                         paste("b", .COLS, sep = "."), collapse = " AND ")

  for(i in 1:nExpData) {
    nextCols <- c(currentCols, nfields[[i]])
    nextTypes <- c(currentTypes, types[[i]])
    statement <- sprintf("CREATE TABLE %s (%s);", temptables[i+1],
                         paste(nextCols, nextTypes, collapse = ", "))
    dbGetQuery(db, statement)
    statement <- sprintf("INSERT INTO %s (%s) SELECT %s, %s FROM %s OUTER LEFT JOIN %s ON %s;",
                         temptables[i+1],
                         paste(nextCols, collapse = ", "),
                         paste("a", currentCols, sep = ".", collapse = ", "),
                         paste("b", efields[[i]], sep = ".", collapse = ","),
                         paste(temptables[i], "AS a"),
                         paste(expdatatables[i], "AS b"),
                         joinstatement)
    currentCols <- nextCols
    currentTypes <- nextTypes
    .timeAndPrint(dbGetQuery(db, statement),
                  paste("Left outer join with table", expdatatables[i]),
                  print = verbose)
    dbGetQuery(db, sprintf("DROP TABLE %s;", temptables[i]))
  }

  ## Optional cleanup
  if(deleteOriginals) {
    sapply(expdatatables, function(x) {
      dbGetQuery(db, sprintf("DROP TABLE %s;", x))
    })
  }

  statement <- sprintf("CREATE INDEX %s ON %s (%s);",
                       paste(tablename, "IDX", sep = ""),
                       tablename, paste(.COLS, collapse = ", "))
  .timeAndPrint(dbGetQuery(db, statement), "Indexing", print = verbose)

  ## Return pointer to new expData
  return(ExpData(db = getDBName(expDataList[[1]]), tablename = tablename, mode = 'w'))
}

