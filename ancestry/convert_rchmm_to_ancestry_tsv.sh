#!/bin/bash
LIST=$1
CHR=$2

par1_string=""
par2_string=""
counter=0
first2_par1=first2_par1.temp_${CHR}
first2_par2=first2_par2.temp_${CHR}
cut -f1 $LIST | while read -r sample; do

	counter=$(($counter + 1))
	focal_file=${sample}.posterior
	currname=${focal_file%.gz*}

	if [ $counter -eq 1 ]; then
		temp_par1=${focal_file}.par1.results
		#echo '\t$currname' > $temp_par1
		#cut -f 1-2,3 $focal_file | tail -n +2 >> $temp_par1
		cut -f 1-2 $focal_file | tail -n +2 >> $first2_par1

		temp_par2=${focal_file}.par2.results
		#echo '\t$currname' > $temp_par2
		#cut -f 1-2,3 $focal_file | tail -n +2 >> $temp_par2
		cut -f 1-2 $focal_file | tail -n +2 >> $first2_par1

		par1_string=${first2_par1}" "${temp_par1}
		par2_string=${first2_par2}" "${temp_par2}

		#temp_par1=$(sed -e 's/^(.+?)\t/\\1:/g' $temp_par1)
		#temp_par2=$(sed -e 's/^(.+?)\t/\\1:/g' $temp_par2)
	else
		temp_par1=${focal_file}.par1.results
		#echo $currname > $temp_par1
		#cut -f 3 $focal_file | tail -n +2 >> $temp_par1

		temp_par2=${focal_file}.par2.results
		#echo $currname > $temp_par2
		#cut -f 5 $focal_file | tail -n +2 >> $temp_par2

		par1_string=${par1_string}" "${temp_par1}
		par2_string=${par2_string}" "${temp_par2}
	fi
	printf 'done processing sample %s chr %s\n' $sample $CHR | tee >(cat >&2)
done
file1name=ancestry-probs-par1_transposed_${CHR}.tsv
file2name=ancestry-probs-par2_transposed_${CHR}.tsv
printf '%s\n' $par1_string
printf '%s\n' $par2_string
paste $par1_string > $file1name
paste $par2_string > $file2name
printf 'done processing all samples to tsv for Chr %s\n' $CHR | tee >(cat >&2)
