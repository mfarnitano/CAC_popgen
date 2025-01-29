#!/bin/bash

###Pulls two samples from a ancestrycalls.txt.chrpos file (columns are samples, rows are sites, values are 0,1,2,NA)
###Counts number of sites with each combination of genotypes
###Summarizes counts to test maternal-offspring relationship and selfing

###Usage: ./compare_samples_ancestry.sh calls_file sample1 sample2 sample_suffix
###Prints single line to standard out
###Option ./compare_samples_ancestry.sh --header ###Prints header line instead

if [ $1 == "--header" ]; then
	printf "sample1\tsample2\ttotal_sites\ttotal_bothcalled\tNAs\t"
	printf "00\t01\t02\t10\t11\t12\t20\t21\t22\t"
	printf "opposites\tnewhets\topposites_prop\tnewhets_perhom\n"
	exit 0
fi

CALLFILE=$1
SAMPLE1=$2
SAMPLE2=$3
SUFFIX=$4
printf "called variables: $1 $2 $3 $4 \n" >&2
#determine which columns to grab
counter=0
S1_COL=""
S2_COL=""
while read -r line; do
	counter=$(( $counter + 1 ))
	if [ $line == ${SAMPLE1}${SUFFIX} ]; then
		S1_COL=$counter
	fi
	if [ $line == ${SAMPLE2}${SUFFIX} ]; then
		S2_COL=$counter
	fi
done <<<$(head -n1 $CALLFILE | tr '\t' '\n')

if [ "$S1_COL" == "" ]; then printf "Sample %s not found in table\n" $SAMPLE1 ; fi
if [ "$S2_COL" == "" ]; then printf "Sample %s not found in table\n" $SAMPLE2 ; fi

#pull columns and count
paste <(cut -f ${S1_COL} ${CALLFILE}) <(cut -f ${S2_COL} ${CALLFILE}) |
	awk -v x=$SAMPLE1 -v y=$SAMPLE2 \
		'BEGIN{OFS="\t";t=0;u=0;n=0;a=0;b=0;c=0;d=0;e=0;f=0;g=0;h=0;i=0}
		{
			if ($1=="NA" || $2=="NA") {t+=1;n+=1}
			else if ($1==0 && $2==0) {t+=1;u+=1;a+=1}
			else if ($1==0 && $2==1) {t+=1;u+=1;b+=1}
			else if ($1==0 && $2==2) {t+=1;u+=1;c+=1}
			else if ($1==1 && $2==0) {t+=1;u+=1;d+=1}
			else if ($1==1 && $2==1) {t+=1;u+=1;e+=1}
			else if ($1==1 && $2==2) {t+=1;u+=1;f+=1}
			else if ($1==2 && $2==0) {t+=1;u+=1;g+=1}
			else if ($1==2 && $2==1) {t+=1;u+=1;h+=1}
			else if ($1==2 && $2==2) {t+=1;;u+=1;i+=1}
		}
		END{print x,y,t,u,n,a,b,c,d,e,f,g,h,i,c+g,b+h,(c+g)/u,(b+h)/(a+b+c+g+h+i)}'
