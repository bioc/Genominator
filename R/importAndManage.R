.makeIndexName <- function(tablename) {
    paste(tablename, "IDX", sep = "")
}

importFromAlignedReads <- function(x, chrMap, dbFilename, tablename,
                                   overwrite = TRUE, deleteIntermediates = TRUE,
                                   readPosition = c("5prime", "left", "center"),
                                   verbose = getOption("verbose"), ...) {
    if (!require(ShortRead))
        stop("ShortRead package must be installed.")
    if(is.null(names(x)) || any(names(x) == ""))
        stop("'x' must have existing, non-empty names")
    readPosition <- match.arg(readPosition)
    weights.AlignedRead <- function(object, ...) {
        if("weights" %in% varLabels(alignData(object))) {
            alignData(object)$"weights"
        } else {
            NULL
        }
    }
    switch(class(x),
           list = {
               if (!all(sapply(x, class) == "AlignedRead") || 
                   anyDuplicated(names(x)))
                   stop("'x' must be a list of 'AlignedRead' objects with unique names.")
               method <- "AlignedReadList"
           },
           character = {
               if (!all(file.exists(x)))
                   stop("'x' must be a named character vector of filenames (of existing files).")
               method <- "filenames"
           },
           stop("seems the 'x' argument is wrong")
           )
    
    importObject <- function(name, aln, verbose) {
        switch(readPosition,
               "5prime" = {
                   loc <- position(aln) + ifelse(strand(aln) == "-", width(aln) - 1, 0)
               },
               "left" = {
                   loc <- position(aln)
               },
               "center" = {
                   loc <- position(aln) + ifelse(strand(aln) == "+", floor(width(aln) / 2) - 1,
                                                 ceiling(width(aln) / 2) - 1)
               })
        str <- c(-1L, 0L, 1L)[match(strand(aln), c("-", "*", "+"))] 
        chr <- match(chromosome(aln), chrMap)
        if("weights" %in% varLabels(alignData(aln))) {
ed <- importToExpData(data.frame(chr = chr, location = loc, 
                    strand = str, weights = alignData(aln)$weights),
                                  dbFilename = dbFilename, tablename = name,
                                  overwrite = TRUE, verbose = verbose)
            aggregateExpData(ed, colname = name, overwrite = TRUE, 
                             verbose = verbose, aggregator = "total(weights)")
        } else {
            ed <- importToExpData(data.frame(chr = chr, location = loc, 
                                             strand = str), dbFilename = dbFilename, tablename = name, 
                                  overwrite = TRUE, verbose = verbose)
            aggregateExpData(ed, colname = name, overwrite = TRUE, 
                             verbose = verbose)
        }
    }
    
    switch(method,
           AlignedReadList = {
               laneNames <- names(x)
               expDataList <- mapply(importObject, laneNames, x, MoreArgs = list(verbose = TRUE))
           },
           filenames = {
               filenames <- x
               columnList <- split(filenames, names(filenames))
               columnNames <- names(columnList)
               if(verbose)
                   cat("filename : column",
                       paste(sapply(names(columnList), function(xx) {
                           paste(paste(columnList[[xx]], ":", xx), collapse = "\n")
                       }), collapse = "\n"), fill = TRUE)
               expDataList <- as.list(columnNames)
               names(expDataList) <- columnNames
               for(name in columnNames) {
                   if(verbose)
                       cat("process column", name, "in ...", fill = TRUE)
                   files <- columnList[[name]]
                   etime <- round(system.time({
                       aln <- readAligned(dirPath = files, ...)
                   })[3], 4)
                   if(verbose)
                       cat("  parsing using ShortRead in", etime, "secs", fill = TRUE)
                   etime <- round(system.time({
                       expDataList[[name]] <- importObject(name = name, aln = aln, verbose = FALSE)
                   })[3], 4)
                   if(verbose)
                       cat("  importing in", etime, "secs", fill = TRUE)
               }
           })
    joinExpData(expDataList, tablename = tablename, overwrite = overwrite, verbose = verbose,
                deleteOriginals = deleteIntermediates)
}

##
## Imports data from a data.frame.
##
importToExpData <- function(df, dbFilename, tablename, overwrite = FALSE,
                            verbose = getOption("verbose"), columns = NULL) {
  db <- dbConnect(dbDriver("SQLite"), dbFilename)
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
  for (i in 1:length(COLS)) {
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
  df <- df[ stats::complete.cases(df[, COLS]), ]
  if(nrow(df) == 0)
      stop("After removing missing locations, df has no rows.")
  df <- df[order(df[,COLS[1]], df[,COLS[2]], df[,COLS[3]]), ]
  
  .timeAndPrint( { if (!dbWriteTable(db, tablename, df, row.names = FALSE, overwrite = overwrite))
                     stop("Unable to write the table, delete table or specify overwrite.")},
                txt = "Writing table", print = verbose)

  idxName <- .makeIndexName(tablename)

  if (overwrite) {
      dbGetQuery(db, paste("DROP INDEX IF EXISTS", idxName))
  }
  q <- sprintf("CREATE INDEX %s ON %s (%s)", idxName, tablename,
               paste(COLS, collapse = ", "))
  .timeAndPrint(dbGetQuery(db, q), txt = "Creating index", print = verbose)
  
  ## clean up -- we will reconnect in the line below.
  dbDisconnect(db)
  
  return(ExpData(dbFilename = dbFilename, tablename, indexColumns = COLS, mode = 'w'))
}

aggregateExpData <- function(expData, by = getIndexColumns(expData), tablename = NULL,
                             deleteOriginal = FALSE, overwrite = FALSE,
                             verbose = getOption("verbose"), colname = "counts",
                             aggregator = paste("count(", by[1], ")", sep = ""))
{
    .checkWrite(expData)

    moveTable <- FALSE

    if (is.null(tablename)) {
        tablename <- sprintf("__tmp_%d", floor(runif(1, 1000, 9999)))
        deleteOriginal <- TRUE
        moveTable <- TRUE
    }

    if (overwrite) {
        dbGetQuery(getDBConnection(expData), paste("DROP TABLE IF EXISTS", tablename))
    }

    cols <- paste(paste(c(by, colname), "INTEGER"), collapse = ",")
    .timeAndPrint(dbGetQuery(getDBConnection(expData), sprintf("CREATE TABLE %s (%s)", tablename, cols)),
                  txt = paste("Creating table:", tablename), print = verbose)

    statement <- sprintf("INSERT INTO %s SELECT %s FROM %s GROUP BY %s",
                         tablename,
                         paste(c(by, aggregator), collapse = ","),
                         getTablename(expData),
                         paste(by, collapse = ","))
    .timeAndPrint(dbGetQuery(getDBConnection(expData), statement),
                  txt = "inserting", print = verbose)

    if (deleteOriginal) {
        .timeAndPrint(dbGetQuery(getDBConnection(expData), paste("DROP TABLE", getTablename(expData))),
                      txt = "droping original table", print = verbose)
    }

    if (moveTable) {
        .timeAndPrint(dbGetQuery(getDBConnection(expData), sprintf("ALTER TABLE %s RENAME TO %s",
                                                         tablename, getTablename(expData))),
                      txt = "renaming table", print = verbose)

        tablename <- getTablename(expData)
    }

    ## now create the index on the correct tablename.
    statement <-  sprintf("CREATE INDEX %s ON %s (%s);", .makeIndexName(tablename),
                          tablename, paste(by, collapse = ","))
    .timeAndPrint(dbGetQuery(getDBConnection(expData), statement),
                  txt = "creating index", print = verbose)

    ## return a new expData.
    return(ExpData(dbFilename = getDBFilename(expData), tablename = tablename,
                   indexColumns = by, mode = 'w'))
}

collapseExpData <- function(expData, tablename = NULL, what = getColnames(expData, all = FALSE),
                            groups = "COL", collapse = c("sum", "avg", "weighted.avg"),
                            overwrite = FALSE, deleteOriginal = FALSE,
                            verbose = getOption("verbose")) {
    .checkWrite(expData)

    ## XXX: Code Duplication!
    moveTable <- FALSE
    if (is.null(tablename)) {
        tablename <- sprintf("__tmp_%d", floor(runif(1, 1000, 9999)))
        deleteOriginal <- TRUE
        moveTable <- TRUE
    }
    if (overwrite) {
        dbGetQuery(getDBConnection(expData), paste("DROP TABLE IF EXISTS", tablename))
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
    .timeAndPrint(dbGetQuery(getDBConnection(expData), statement),
                  txt = "creating table", print = verbose, query = statement)

    statement <- sprintf("INSERT INTO %s SELECT %s FROM %s GROUP BY %s",
                         tablename,
                         paste(c(getIndexColumns(expData), sel), collapse = ", "), getTablename(expData),
                         paste(getIndexColumns(expData), collapse = ","), paste(getIndexColumns(expData), collapse = ","))
    .timeAndPrint(dbGetQuery(getDBConnection(expData), statement),
                  txt = "inserting data", print = verbose, query = statement)

    ## XXX: Code duplication
    if (deleteOriginal) {
        .timeAndPrint(dbGetQuery(getDBConnection(expData), paste("DROP TABLE", getTablename(expData))),
                      txt = "droping original table", print = verbose)
    }
    if (moveTable) {
        .timeAndPrint(dbGetQuery(getDBConnection(expData), sprintf("ALTER TABLE %s RENAME TO %s",
                                                         tablename, getTablename(expData))),
                      txt = "renaming table", print = verbose)
        tablename <- getTablename(expData)
    }

    ## now create the index on the correct tablename.
    statement <-  sprintf("CREATE INDEX %s ON %s (%s);", .makeIndexName(tablename),
                          tablename, paste(getIndexColumns(expData), collapse = ","))
    .timeAndPrint(dbGetQuery(getDBConnection(expData), statement),
                  txt = "creating index", print = verbose)
    
    ## return a new ExpData
    return(ExpData(dbFilename = getDBFilename(expData), tablename = tablename, mode = 'w')) 
}

## Infer the column types. pretty lame.
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
                        overwrite = TRUE, deleteOriginals = FALSE,
                        verbose = getOption("verbose")){
    .checkWrite(expDataList)
    
  if(class(expDataList) != "list" || length(expDataList) < 2)
    stop("argument 'expDataList' must be a list of at least 2 components")
  if(!length(unique((sapply(expDataList, getDBFilename)))))
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
  db <- getDBConnection(expDataList[[1]])

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
  .timeAndPrint(dbGetQuery(db, statement),
                txt = "Creating union", print = verbose)

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
                  txt = paste("Left outer join with table", expdatatables[i]),
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
  .timeAndPrint(dbGetQuery(db, statement), txt = "Indexing", print = verbose)

  ## Return pointer to new expData
  return(ExpData(dbFilename = getDBFilename(expDataList[[1]]), tablename = tablename, mode = 'w'))
}

