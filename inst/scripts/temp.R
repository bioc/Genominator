## For yeast1
library(RSQLite)
dbfile <- paste(Sys.getenv("COLLAB"), "/extdata-yeastTiling/UHTS/yeast1.db", sep = "")
expData <- new("ExpData.db", db = dbConnect(dbDriver("SQLite"), dbfile), tablename = "rawReads")
expData <- aggregateBy(expData)

## For restoring in case something is screwed up
debug(aggregateBy)
dbGetQuery(getDB(expData), "DROP TABLE rawCounts;")
expData@tablename <- "rawReads"

## For marioni
library(RSQLite)
dbfile <- paste(Sys.getenv("COLLAB"), "/extdata-solexa/marioni/marioni.db", sep = "")
expData <- new("ExpData.db", db = dbConnect(dbDriver("SQLite"), dbfile), tablename = "rawReads")
system.time(expData <- aggregateBy(expData, by = c("chr", "strand", "location", "lane"),
                       database = paste(Sys.getenv("COLLAB"),
                       "/extdata-solexa/marioni/marioni_agg.db", sep = "")))



## Some stuff for Marioni

library(biomaRt)
mart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
annot <- getBM(c("structure_gene_stable_id",
                 "structure_transcript_stable_id",
                 "structure_exon_stable_id",
                 "structure_chrom_name",
                 "structure_exon_chrom_start",
                 "structure_exon_chrom_end",
                 "structure_exon_rank",
                 "structure_biotype"),
               filter = "chromosome_name", values = 21, mart = mart)
