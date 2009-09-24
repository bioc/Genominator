#!/bin/bash

INDIR=$1
OUTDIR=$2
DBNAME=$3
TBLNAME=$4

SCRIPTS_DIR=`dirname $0`

if [ "$INDIR" = "" ] | [ "$OUTDIR" = "" ] | [ "$DBNAME" = "" ] | [ "$TBLNAME" = "" ]
then
    echo "$0 <INDIR> <OUTDIR> <DBNAME> <TBLNAME>"
    exit 1
fi

rm -f $DBNAME.reshaped

echo "reshaping"
for f in `ls $INDIR`;
do
  echo "\t $f"
  ${SCRIPTS_DIR}/reshape.pl $INDIR/$f >> $DBNAME.reshaped
done

echo "sorting"
${SCRIPTS_DIR}/sortreshape.sh $DBNAME.reshaped $DBNAME.sorted

echo "importing"
${SCRIPTS_DIR}/sqlimport.sh $DBNAME.sorted $OUTDIR $DBNAME.db $TBLNAME

# rm $DBNAME.reshaped $DBNAME.sorted
exit 0
