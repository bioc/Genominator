#!/bin/bash

INDIR=${COLLAB}/extdata-yeastTiling/UHTS
OUTDIR=${INDIR}
INFILE=080129_FC2079B_s4__SL152.align_25.scer_Feb08.txt
DBNAME=yeast
TBLNAME=rawReads

SCRIPTS_DIR=`dirname $0`

if [ -f ${OUTDIR}/${DBNAME}.reshaped ] | [ -f ${OUTDIR}/${DBNAME}.sort ] | [ -f ${OUTDIR}/${DBNAME}.db ]
then
  echo "There are old files in ${OUTDIR}, please remove"
  exit 1
fi

echo "..reshaping"
time ${SCRIPTS_DIR}/reshape-yeast1.pl ${INDIR}/${INFILE} >> ${OUTDIR}/${DBNAME}.reshaped
${SCRIPTS_DIR}/sqlimport.sh ${OUTDIR}/${DBNAME}.reshaped ${OUTDIR} ${DBNAME}.db ${TBLNAME}

rm ${OUTDIR}/${DBNAME}.reshaped
exit 0
