#!/bin/bash

## This script assumed that we have a "reshaped" read file, ie.
## a file with four columns separated by "|". The column order is
##   chromosome, strand, location, "other stuff"

INFILE=$1
OUTDIR=$2
DBNAME=$3
TBLNAME=rawReads
OUTRESHAPED=${OUTDIR}/${DBNAME}.reshaped
OUTSORT=${OUTDIR}/${DBNAME}.sort
OUTDB=${OUTDIR}/${DBNAME}


# echo "..sorting"
# time sort -t '|' -k 1,1 -k 2,2 -k 3,3 ${INFILE} > ${OUTSORT}

sqlite3 ${OUTDB} "CREATE TABLE IF NOT EXISTS ${TBLNAME} (chr INTEGER, strand INTEGER, location INTEGER, quality INTEGER);"
echo "..importing"
time sqlite3 ${OUTDB} ".import ${OUTSORT} ${TBLNAME}"
echo "..indexing"
time sqlite3 ${OUTDB} "CREATE INDEX chrstrandlocation${TBLNAME}IDX ON ${TBLNAME} (chr, strand, location);"
echo "..analyzing";
time sqlite3 ${OUTDB} "ANALYZE $TBLNAME;"

rm ${OUTSORT}

exit 0;
