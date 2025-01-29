#!/bin/bash

QUALIMAP_DIRECTORY=$1
PREFIX=$2
OUTFILE=$3

INFILE=${QUALIMAP_DIRECTORY}/${PREFIX}/genome_results.txt

if [ ! -f $INFILE ]; then printf "%s\tNO DATA AVAILABLE AT %s\n" $PREFIX $INFILE >> $OUTFILE; exit; fi

#PREFIX
printf "%s\t" $PREFIX >> $OUTFILE

#MAPPED READS
grep "number of reads = " $INFILE | tr -s ' ' | cut -d' ' -f6 | tr -d ',' | tr -d '\n' >> $OUTFILE
printf "\t" >> $OUTFILE

#MEAN COVERAGE
grep "mean coverageData = " $INFILE | tr -s ' ' | cut -d' ' -f5 | tr -d 'X' | tr -d '\n' >> $OUTFILE
printf "\t" >> $OUTFILE

#STD COVERAGE
grep "std coverageData = " $INFILE | tr -s ' ' | cut -d' ' -f5 | tr -d 'X' | tr -d '\n' >> $OUTFILE
printf "\t" >> $OUTFILE

#BREADTH OF COVERAGE
grep "coverageData >= 1X" $INFILE | tr -s ' ' | cut -d' ' -f5 | tr -d '\n' >> $OUTFILE
printf "\t" >> $OUTFILE
grep "coverageData >= 2X" $INFILE | tr -s ' ' | cut -d' ' -f5 | tr -d '\n' >> $OUTFILE
printf "\t" >> $OUTFILE
grep "coverageData >= 3X" $INFILE | tr -s ' ' | cut -d' ' -f5 | tr -d '\n' >> $OUTFILE
printf "\t" >> $OUTFILE
grep "coverageData >= 5X" $INFILE | tr -s ' ' | cut -d' ' -f5 | tr -d '\n' >> $OUTFILE
printf "\t" >> $OUTFILE
grep "coverageData >= 10X" $INFILE | tr -s ' ' | cut -d' ' -f5 | tr -d '\n' >> $OUTFILE
printf "\t" >> $OUTFILE

#FINISH LINE
printf "\n" >> $OUTFILE
