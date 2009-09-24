library(RSQLite)
source("../../R/sequalizer.R")
source("../../R/processSequencingRun.R")

RunOneDir <- paste(Sys.getenv("COLLAB"), "extdata-solexa/marioni/ForAnalysisRunOne", sep = "/")
RunTwoDir <- paste(Sys.getenv("COLLAB"), "extdata-solexa/marioni/ForAnalysisRunTwo", sep = "/")
EndDir <- paste(Sys.getenv("COLLAB"), "extdata-solexa/marioni", sep = "/")

## First we aggregate
## Comment: bad that we need to disconnect expData, should be fixed somehow

for(database in list.files(RunOneDir, pattern = "s_..db$", full.names = TRUE)) {
    expData <- new("ExpData.db", db = database, tablename = "rawReads")
    aggregateExpData(expData, by = c("chr", "strand", "location"), tablename = "rawCounts")
    dbDisconnect(getDB(expData))
}

for(database in list.files(RunTwoDir, pattern = "s_..db$", full.names = TRUE)) {
    expData <- new("ExpData.db", db = database, tablename = "rawReads")
    aggregateExpData(expData, by = c("chr", "strand", "location"), tablename = "rawCounts")
    dbDisconnect(getDB(expData))
}

## Now we merge. This is painful at the moment, but I am trying to figure out how to do it best...
## The code below is just for ForAnalysisRunOne

databases <- list.files(RunOneDir, pattern = "s_..db$", full.names = TRUE)
expData <- new("ExpData.db", db = databases[grep("s_1", databases)], tablename = "rawCounts")
expData2 <- new("ExpData.db", db = databases[grep("s_2", databases)], tablename = "rawCounts")
expOut <- new("ExpData.db", db = paste(EndDir, "cell1.db", sep = "/"))

expOut <- mergeExpData(exp1 = expData, exp2 = expData2, expout = expOut, tablename = "s1s2", fields1 = c("counts" = "counts_1_s1"),
                       fields2 = c("counts" = "counts_1_s2"))
rm(expData2)

expData <- new("ExpData.db", db = databases[grep("s_3", databases)], tablename = "rawCounts")
expOut <- mergeExpData(exp1 = expOut, exp2 = expData, tablename = "s1s2s3", fields2 = c("counts" = "counts_1_s3"))
.dbGetAndTime(getDB(expOut), "DROP TABLE s1s2;")

expData <- new("ExpData.db", db = databases[grep("s_4", databases)], tablename = "rawCounts")
expOut <- mergeExpData(exp1 = expOut, exp2 = expData, tablename = "s1s2s3s4", fields2 = c("counts" = "counts_1_s4"))
.dbGetAndTime(getDB(expOut), "DROP TABLE s1s2s3;")

expData <- new("ExpData.db", db = databases[grep("s_6", databases)], tablename = "rawCounts")
expOut <- mergeExpData(exp1 = expOut, exp2 = expData, tablename = "s1s2s3s4s6", fields2 = c("counts" = "counts_1_s6"))
.dbGetAndTime(getDB(expOut), "DROP TABLE s1s2s3s4;")

expData <- new("ExpData.db", db = databases[grep("s_7", databases)], tablename = "rawCounts")
expOut <- mergeExpData(exp1 = expOut, exp2 = expData, tablename = "s1s2s3s4s6s7", fields2 = c("counts" = "counts_1_s7"))
.dbGetAndTime(getDB(expOut), "DROP TABLE s1s2s3s4s6;")

expData <- new("ExpData.db", db = databases[grep("s_8", databases)], tablename = "rawCounts")
expOut <- mergeExpData(exp1 = expOut, exp2 = expData, tablename = "s1s2s3s4s6s7s8", fields2 = c("counts" = "counts_1_s8"))
.dbGetAndTime(getDB(expOut), "DROP TABLE s1s2s3s4s6s7;")

## Now ForAnalysisRunTwo

databases <- list.files(RunTwoDir, pattern = "s_..db$", full.names = TRUE)
expData <- new("ExpData.db", db = databases[grep("s_1", databases)], tablename = "rawCounts")
expData2 <- new("ExpData.db", db = databases[grep("s_2", databases)], tablename = "rawCounts")
expOut <- new("ExpData.db", db = paste(EndDir, "cell2.db", sep = "/"))

expOut <- mergeExpData(exp1 = expData, exp2 = expData2, expout = expOut, tablename = "s1s2", fields1 = c("counts" = "counts_2_s1"),
                       fields2 = c("counts" = "counts_2_s2"))
rm(expData2)

expData <- new("ExpData.db", db = databases[grep("s_3", databases)], tablename = "rawCounts")
expOut <- mergeExpData(exp1 = expOut, exp2 = expData, tablename = "s1s2s3", fields2 = c("counts" = "counts_2_s3"))
.dbGetAndTime(getDB(expOut), "DROP TABLE s1s2;")

expData <- new("ExpData.db", db = databases[grep("s_4", databases)], tablename = "rawCounts")
expOut <- mergeExpData(exp1 = expOut, exp2 = expData, tablename = "s1s2s3s4", fields2 = c("counts" = "counts_2_s4"))
.dbGetAndTime(getDB(expOut), "DROP TABLE s1s2s3;")

expData <- new("ExpData.db", db = databases[grep("s_6", databases)], tablename = "rawCounts")
expOut <- mergeExpData(exp1 = expOut, exp2 = expData, tablename = "s1s2s3s4s6", fields2 = c("counts" = "counts_2_s6"))
.dbGetAndTime(getDB(expOut), "DROP TABLE s1s2s3s4;")

expData <- new("ExpData.db", db = databases[grep("s_7", databases)], tablename = "rawCounts")
expOut <- mergeExpData(exp1 = expOut, exp2 = expData, tablename = "s1s2s3s4s6s7", fields2 = c("counts" = "counts_2_s7"))
.dbGetAndTime(getDB(expOut), "DROP TABLE s1s2s3s4s6;")

expData <- new("ExpData.db", db = databases[grep("s_8", databases)], tablename = "rawCounts")
expOut <- mergeExpData(exp1 = expOut, exp2 = expData, tablename = "s1s2s3s4s6s7s8", fields2 = c("counts" = "counts_2_s8"))
.dbGetAndTime(getDB(expOut), "DROP TABLE s1s2s3s4s6s7;")

## Now we merge the two cell databases

dbDisconnect(getDB(expOut)) # This is just to gracefully close the object because we will overwrite it

cellfiles <- list.files(EndDir, pattern = "cell..db", full.names = TRUE)
expCell1 <- new("ExpData.db", db = cellfiles[1], tablename = "s1s2s3s4s6s7s8")
expCell2 <- new("ExpData.db", db = cellfiles[2], tablename = "s1s2s3s4s6s7s8")
expOut <- new("ExpData.db", db = paste(EndDir, "marioni_pivot.db", sep = "/"))
expOut <- mergeExpData(exp1 = expCell1, exp2 = expCell2, expout = expOut)
