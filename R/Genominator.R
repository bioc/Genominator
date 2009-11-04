
setClass("ExpData",
         representation(db           = "character",
                        tablename    = "character",
                        tableSchema  = "character",
                        indexColumns = "character",
                        mode         = "character",
                        chrMap       = "character",
                        .tmpFile     = "character",
                        .pool        = "environment"
                        ),
         prototype(indexColumns = c("chr", "location", "strand"),
                   .pool = new.env(hash=TRUE, parent = emptyenv())))

setMethod("show", "ExpData", function(object) {
    cat("table:", getTablename(object), "\n")
    cat("database file:", getDBName(object), "\n")
    cat("index columns:", getIndexColumns(object), "\n")
    cat("mode:", getMode(object), "\n")
    cat("schema:\n")
    print(getSchema(object))
})

setMethod("head", "ExpData", function(x, ...) {
    args <- list(...)
    if (length(args) == 1 && is.numeric(args[[1]]))
        n <- args[[1]]
    else
        n <- 10
    dbGetQuery(getDB(x), sprintf("SELECT * FROM %s LIMIT %s;", getTablename(x), n))
})

setMethod("initialize", signature(.Object = "ExpData"),
          function(.Object, db, tablename = "", pragmas = NULL, indexColumns = NULL,
                   mode = c('r', 'w')) {
            if (missing(db))
              stop("Must specify db argument.")
            
            db <- Sys.glob(db)
            
            .Object@mode <- match.arg(mode)
            .Object@db <- db
            .Object@tablename <- tablename
             
            if (.Object@mode == 'r')
              .Object@.tmpFile <- tempfile()
            
            if (!is.null(indexColumns))
              .Object@indexColumns <- indexColumns

            ##
            ## this is a simple sharing mechanism so that we don't have to think about
            ## instanteating expData objects and overloading the number of connections
            ## a database can have.

            ## XXX : think about the fact that you can't see the schema
            ##       and the fact that something you do might make the db unsafe w/multiple
            ##       threads accessing it. 
            ##
            if (length(.Object@db) != 0 && .Object@db != "") {
              if (exists(.Object@db, envir = .Object@.pool)) {
                ## refresh each time we call this function, could have stale
                ## handles if you are borrowing the connection from an object. 
                dbDisconnect(getDB(.Object))
              }
            }
            
            if (.Object@mode == 'r') 
              v <- .Object@.tmpFile
            else
              v <- .Object@db

            ## connect.
            assign(.Object@db, dbConnect(dbDriver("SQLite"), v), .Object@.pool)

            ## now set the pragmas
            if (!is.null(pragmas)) {
              dbGetQuery(getDB(.Object), pragmas)
            }

            ## now attach the data if we are in read mode. I'll also
            ## want to register a finalizer somehow.
            if (.Object@mode == 'r') {
              dbGetQuery(getDB(.Object), sprintf("attach \"%s\" as \"\";", .Object@db))
            }
            
            reg.finalizer(.Object@.pool, function(o) {
              dbDisconnect(getDB(.Object))
              if (.Object@mode == 'r') unlink(.Object@.tmpFile)
            })
            
              
            .Object@tableSchema <- dbListFieldsAndTypes(getDB(.Object), .Object@tablename)
            return(.Object)
          })


ExpData <- function(db, tablename, mode = c('r', 'w'), indexColumns = NULL, pragmas = NULL) {
  if (is.null(indexColumns)) 
    new("ExpData", db = db, tablename = tablename, pragmas = pragmas,
        mode = mode)
  else
    new("ExpData", db = db, tablename = tablename, pragmas = pragmas,
        indexColumns = indexColumns, mode = mode)
}

getIndexColumns <- function(expData) {
    return(expData@indexColumns)
}

getMode <- function(expData) {
  return(expData@mode)
}

getDB <- function(expData) {
  get(expData@db, expData@.pool)
}

getDBName <- function(expData) {
    expData@db
}

getTablename <- function(expData) {
    expData@tablename
}

getSchema <- function(expData) {
  expData@tableSchema
}

getColnames <- function(expData, all = TRUE) {
  cc <- names(getSchema(expData))
  if (all) {
    return(cc)
  } else {
    return(setdiff(cc, getIndexColumns(expData)))
  }
}

listTables <- function(db) { 
  return(dbListTables(dbConnect(dbDriver("SQLite"), db)))
}

setMethod("$", signature = "ExpData", definition = function(x, name) {
    dbGetQuery(getDB(x), sprintf("SElECT %s from %s", name, getTablename(x)))
})

##
## The semantics here is slightly off because there are really no
## "rownames", however there should be a way to make this more
## general.
##
setMethod("[", signature = "ExpData", definition = function (x, i, j, ..., drop = TRUE) {
  if (!missing(j)) {
      cols <- switch(class(j),
                     "character" = j,
                     "integer" = getColnames(x, all = TRUE)[j],
                     "numeric" = getColnames(x, all = TRUE)[j],
                     "logical" = getColnames(x, all = TRUE)[j],
                     stop("Column subsetting object needs to be either a character, a numeric/integer or a logical."))
      whichClause <- paste(cols, collapse = ", ")
  } else {
      whichClause <- "*"
  }
  if (!missing(i)) {
      whereClause <- switch(class(i),
                            "character" = stop("ExpData objects do not have rownames, subsetting makes no sense."),
                            "integer" = sprintf("WHERE _ROWID_ in (%s)", paste(i, collapse = ", ")),
                            "numeric" = sprintf("WHERE _ROWID_ in (%s)", paste(i, collapse = ", ")),
                            "logical" = sprintf("WHERE _ROWID_ in (%s)", paste(which(i), collapse = ", ")),
                            stop("Row subsetting object needs to be either a numeric/integer or a logical.")
                            )
  } else {
      whereClause <- ""
  }
  q <- sprintf("SELECT %s FROM %s %s", whichClause, getTablename(x), whereClause)
  if ("verbose" %in% names(list(...)))
      print(q)
  dbGetQuery(getDB(x), q)
})

getRegion <- function(expData, chr, start, end, strand, what = "*",
                      whereClause = "", verbose = FALSE) {
    if (what[1] != "*") {
        what <- paste(what, collapse = ",")
    }
    if (missing(strand) || strand == 0) {
        strand <- c(-1, 0, 1)
    }
    
    if (missing(start))
      start <- 0   ## this is better than excluding it because of the index.
    if (missing(end))
      end <- 1e12  ## i imagine this is longer than anything around. 

    if (whereClause != "")
      whereClause <- paste("AND", whereClause)
    
    q <- sprintf(paste("SELECT %s FROM %s WHERE chr = %s AND (strand IN (%s) OR strand = 0) AND",
                       "location between %s AND %s %s ORDER BY %s"),
                 what, getTablename(expData), chr, paste(strand, collapse = ","),
                 start, end, whereClause, paste(getIndexColumns(expData), collapse = ","))
    
    .timeAndPrint(res <- dbGetQuery(getDB(expData), q),
                  "fetching region query", print = verbose, q)
    return(res)
}

splitByAnnotation <- function(expData, annoData, what = "*", ignoreStrand = FALSE,
                              expand = FALSE, addOverStrands = FALSE, verbose = FALSE)
{
  if (is.null(rownames(annoData)))
    stop("need uniqe rownames for the annoData object.")
  
  if (expand) {
    ## In this case, I prepend the annotation columns. This adds
    ## additional size to the database query, but significantly
    ## clarifys the code.
    what <- c(getIndexColumns(expData), what)
  }

  if (what[1] != "*") {
    what <- paste(paste(getTablename(expData), what, sep = "."), collapse = ", ")
  } else {
    what <- paste(getTablename(expData), ".*", sep = "")
  }
  
  .timeAndPrint(.writeRegionsTable(expData, annoData), "writing region table", print = verbose)
  
  ## here we need to order things a certain way so we are guaranteed to have things
  ## come out in a consistent way.
  regionID <- paste(.REGION.TABLE.NAME, ".id", sep = "")
  
  oby <- paste(regionID, paste(paste(getTablename(expData), getIndexColumns(expData), sep = "."),
                               collapse = ", "),
               sep = ", ")
  
  q <- .formRegionsSQL(paste(what, ",", regionID), getTablename(expData),
                       ignoreStrand = ignoreStrand, oby = oby)
  .timeAndPrint(tbl <- dbGetQuery(getDB(expData), q),
                "fetching splits table", print = verbose, q = q)
  
  if (nrow(tbl) == 0 || ncol(tbl) == 0) {
    return(NULL)
  }
  
  ## count query to determine size.
  q <- .formRegionsSQL(sprintf("count(%s), %s", regionID, regionID), getTablename(expData),
                       gby = regionID, ignoreStrand = ignoreStrand,
                       oby = regionID)
  .timeAndPrint(cdb <- dbGetQuery(getDB(expData), q),
                "count query", print = verbose, q = q)
  
  ## this makes things worlds faster, however it also means that what is in
  ## the database must be numbers or all become characters. 
  tbl <- as.matrix(tbl)
  
  res <- vector("list", nrow(cdb))
  names(res) <- rownames(annoData)[cdb[,2]]
  lens <- cdb[,1]; clens <- cumsum(lens)
  bounds <- cbind(c(1, 1 + clens[-length(lens)]), clens)
  
  .timeAndPrint({ for (i in seq_len(nrow(bounds))) {
    ## I am dropping the ID column.
    res[[i]] <- tbl[bounds[i,1]:bounds[i,2], -ncol(tbl), drop = FALSE]
  }}, "performing split", print = verbose)
  
  ## here we expand the region with 0s -- this can considerably decrease
  ## performance.
  if (expand) {
    expandF <- function(oMat, reg) {
      strands <- as.character(sort(unique(oMat[,"strand"])))
      byStrand <- vector("list", length(strands))
      names(byStrand) <- strands
      
      for (strand in strands) {
        cols <- oMat[strand == oMat[,"strand"], 1:length(getIndexColumns(expData)), drop = FALSE]
        mat <- oMat[strand == oMat[,"strand"], -(1:length(getIndexColumns(expData))), drop = FALSE]
        
        len <- reg[,"end"] - reg[,"start"] + 1
        mm <- matrix(0, nrow = len, ncol = ncol(mat))
        colnames(mm) <- colnames(mat)                
        mm[cols[, "location"] - reg[,"start"] + 1, ] <- mat[cols[,"strand"] == strand, ]
        
        ## convert NAs to 0, maybe not right?
        mm[is.na(mm)] <- 0
        
        if ("location" %in% colnames(mm))
          mm[,"location"] <- as.integer(reg[,"start"]:reg[,"end"])
        if ("strand" %in% colnames(mm))
          mm[,"strand"] <- as.integer(strand)
        if ("chr" %in% colnames(mm)) {
          mm[,"chr"] <- as.integer(reg[,"chr"])
        }
        byStrand[[strand]] <- mm
      }
      retVal <- byStrand
      
      ## you expand and you only have one strand.
      ## if (length(byStrand) == 1 && ("strand" %in% colnames(byStrand[[1]])))
      ##    retVal <- byStrand[[1]][, -which("strand" == colnames(byStrand[[1]]))]
      
      if (addOverStrands) {
        oCols <- setdiff(colnames(byStrand[[1]]), getIndexColumns(expData))
        cCols <- getIndexColumns(expData)[getIndexColumns(expData) %in% colnames(byStrand[[1]])]
        
        if (length(oCols) >= 1) {
          common <- byStrand[[1]][, cCols]
        }
        tmp <- Reduce('+', lapply(byStrand, function(a)  a[, oCols]))
        
        byStrand <- cbind(common, tmp)
        colnames(byStrand) <- c(cCols, oCols)
        
        if ("strand" %in% colnames(byStrand)) {
          retVal <- byStrand[, -which("strand" == colnames(byStrand))]
        }
        else {
          retVal <- byStrand
        }
      }
      return(retVal)
    }
    res <- applyMapped(res, annoData, expandF)
  }
  return(res)
}

mergeWithAnnotation <- function(expData, annoData, what = "*", ignoreStrand = FALSE, splitBy = NULL,
                                verbose = FALSE) {
  
  .timeAndPrint(.writeRegionsTable(expData, annoData, dropCols = FALSE),
                "writing regions table", print = verbose)
  
  if (what[1] != "*") {
    what <- paste(what, collapse = ",")
  }
  if (!is.null(splitBy) & (what != "*")) {
    ##-- here i need to add the splitBy to the what and then remove it.
    what <- paste(c(what, splitBy), collapse = ",")
  }
  q <- .formRegionsSQL(what, getTablename(expData), ignoreStrand = ignoreStrand)
  .timeAndPrint(tbl <- dbGetQuery(getDB(expData), q), "fetching merge table", print = verbose, q)
  if (!is.null(splitBy)) {
    .timeAndPrint(tbl <- split(tbl[,-ncol(tbl)], tbl[, splitBy]), paste("splitting by:", splitBy),
                  print = verbose)
  }
  return(tbl)
}

summarizeExpData <- function(expData, what = getColnames(expData, all = FALSE), fxs = c("TOTAL"),
                             preserveColnames = TRUE, whereClause = "",
                             verbose = FALSE) {
    originalWhat <- what
    
    what <- paste(sapply(fxs, function(f) {
        paste(f, "(", what, ")", sep = "")
    }), collapse = ", ")

    if (whereClause != "") {
        whereClause <- paste("WHERE", whereClause)
    }
    
    q <- sprintf("SELECT %s FROM %s %s;", what, getTablename(expData), whereClause)
    .timeAndPrint(tbl <- dbGetQuery(getDB(expData), q), "fetching summary", print = verbose)

    if (preserveColnames & length(originalWhat) == ncol(tbl)) {
        colnames(tbl) <- originalWhat
    }

    return(drop(as.matrix(tbl)))
}

summarizeByAnnotation <- function(expData, annoData, what = getColnames(expData, all = FALSE),
                                  fxs = c("TOTAL"), ignoreStrand = FALSE, splitBy = NULL,
                                  bindAnno = FALSE, preserveColnames = TRUE, verbose = FALSE)
{
    .timeAndPrint(ids <- .writeRegionsTable(expData, annoData),
                  "writing regions table", print = verbose)

    if (length(fxs) > 1 & preserveColnames) {
        if (!missing(preserveColnames))
            warning("Cannot preserve column names when you are applying more than one function to the columns.")
        preserveColnames <- FALSE
    }
    
    originalWhat <- what
    what <- paste(sapply(fxs, function(f) {
        paste(f, "(", what, ")", sep = "")
    }), collapse = ", ")

    what <- paste("__regions__.id", what, sep = ", ")
    
    q <- sprintf(paste("SELECT %s FROM __regions__ LEFT OUTER JOIN %s ON __regions__.chr = %s.chr",
                       "AND %s.location BETWEEN __regions__.start AND __regions__.end__1",
                       if(!ignoreStrand) {
                           sprintf("AND (%s.strand = __regions__.strand OR __regions__.strand = 0 OR %s.strand = 0)",
                                   getTablename(expData), getTablename(expData))
                       }, "GROUP BY __regions__.id ORDER BY __regions__.id"),
                 what, getTablename(expData), getTablename(expData),
                 getTablename(expData))
    
    .timeAndPrint(tbl <- dbGetQuery(getDB(expData), q), "fetching summary table", print = verbose, q)
    rownames(tbl) <- rownames(annoData)
    tbl <- tbl[, -1, drop = FALSE]
    
    if (preserveColnames)
        colnames(tbl) <- originalWhat

    if (bindAnno) {
        tbl <- cbind(annoData, tbl)
    }


    if (is.null(splitBy)) {
        return(tbl)
    } else {
        return(split(tbl, annoData[, splitBy], drop = TRUE))
    }
    
}

applyMapped <- function(mapped, annoData, FUN, bindAnno = FALSE) {
    annoData <- annoData[names(mapped),]

    x <- lapply(1:length(mapped), function(i) {
        FUN(mapped[[i]], annoData[i,])
    })
    names(x) <- rownames(annoData)
    
    if (bindAnno) {
        x <- do.call(rbind, x)
        x <- cbind(annoData, x)
    }
    else {
        names(x) <- names(mapped)
    }

    return(x)
}

.timeAndPrint <- function(exp, txt, print = TRUE, q = NULL) {
    if (print) {
        if (!is.null(q))
            cat(q, "\n")
        cat(txt, ": ", sep = "")
    }
    time <- round(system.time(exp)[3], 4)

    if (print) {
        cat(time, "sec\n")
    }
}

.formRegionsSQL <- function(what, tablename, gby = NULL, ignoreStrand = FALSE,
                            oby = NULL, regionTableName = .REGION.TABLE.NAME) {
  
  s <- paste(sprintf(paste("SELECT %s FROM %s INNER JOIN %s ON %s.chr = %s.chr AND",
                           "%s.location BETWEEN %s.start AND %s.end__1"),
                     what, tablename, regionTableName, regionTableName, tablename,
                     tablename, regionTableName, regionTableName),
             if(!ignoreStrand) {
               sprintf("AND (%s.strand = %s.strand OR %s.strand = 0 OR %s.strand = 0)",
                       tablename, regionTableName, regionTableName, tablename)
             } else {
               ""
             })
    
  if (!is.null(gby)) {
    s <- paste(s, sprintf("GROUP BY %s", gby))
  }
  if (!is.null(oby)) {
    s <- paste(s, sprintf("ORDER BY %s", oby))
  }
  return(s)
}

## XXX: I have to deal with this name munging which turns end into end__1.
.writeRegionsTable <- function(expData, annoData, dropCols = TRUE) {
  id <- as.integer(1:nrow(annoData))

  joinCols <- annoData[, .ANNO.COLS]
  for (i in 1:ncol(joinCols)) {
      if (class(joinCols[,i]) != "integer")
          joinCols[,i] <- as.integer(joinCols[,i])
  }
  
  if (dropCols) {
    regions <- cbind(id = id, joinCols)
  } else {
    regions <- cbind(id = id, joinCols, annoData[, -which(.ANNO.COLS %in% colnames(annoData)), drop = FALSE])
  }
  con <- getDB(expData)
  
  ## This approach uses the dbWriteTable
  dbWriteTable(con, .REGION.TABLE.NAME, regions, overwrite = TRUE, row.names = FALSE)
  
  ## This approach uses the temporary table -- this doesn't appear faster and is
  ## significantly more code, but I haven't done sufficient testing.
  ## dbWriteTable(con, .REGION.TABLE.TMP.NAME, regions, row.names = FALSE, overwrite = TRUE)
  ## tryCatch(dbGetQuery(con, sprintf("DROP TABLE %s", .REGION.TABLE.NAME)), error = function(a) {})
  ## tmp <- colnames(regions)
  ## tmp[tmp == "end"] <- "end__1"
  ## cols <- paste(paste(tmp, "INTEGER"), collapse = ",")
  ## sql <- sprintf("CREATE TEMPORARY TABLE %s (%s)", .REGION.TABLE.NAME, cols)
  ## dbGetQuery(con, sql)
  ## sql <- sprintf("INSERT INTO %s SELECT * FROM %s", .REGION.TABLE.NAME, .REGION.TABLE.TMP.NAME)
  ## dbGetQuery(con, sql)

  return(id)
}
