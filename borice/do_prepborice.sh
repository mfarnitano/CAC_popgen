#!/bin/bash
#SBATCH --job-name=prepborice                     # Job name
#SBATCH --partition=batch	                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=9		                            # Number of cores per task
#SBATCH --mem=80gb			                                # Total memory for job
#SBATCH --time=24:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/borice_hq/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/borice_hq/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

###SETUP
WORKING_DIR=$1
BORICE_BAMINFO=$2 #col1=bam full paths, col2=family ID, col3=par or off
BEAGLE_INPUT=$3 #beagle file from angsd
OUTPREFIX=$4 #final file for borice input will be ${OUTPREFIX}.genotypes.txt

cd ${WORKING_DIR}

cut -f1 ${BORICE_BAMINFO} > ${WORKING_DIR}/${OUTPREFIX}_bamlist.txt

BORICE_IDS=${WORKING_DIR}/${OUTPREFIX}_boriceIDs.txt
if [ ! -f $BORICE_IDS ]; then
	printf "SNP\n" > temp.txt
	cut -f2,3 $BORICE_BAMINFO > namestemp.txt
	cat temp.txt namestemp.txt > $BORICE_IDS
	rm temp.txt
	rm namestemp.txt
fi

zcat ${BEAGLE_INPUT} | tail -n +2 |
	awk '{for(i=1;i<=NF;i++)if($i~/^[0-9]+$/&&i==2||i==3){if($i==0) $i="A"; else if($i==1) $i="C"; else if($i==2) $i="G"; else if($i==3) $i="T"}}1' |
	awk '{for(i=1;i<=NF;i+=3) print $i,$(i+1),$(i+2)}' | tr -s ' ' | sed -e 's/ /\t/g' > ${WORKING_DIR}/${OUTPREFIX}.boricelikes.txt

SCRIPTS_DIR=~/Git_repos/Tn5_sequencing/borice_hq
python ${SCRIPTS_DIR}/likes2borice.py ${WORKING_DIR}/${OUTPREFIX}.boricelikes.txt ${BORICE_IDS} > ${OUTPREFIX}.genotypes.txt
