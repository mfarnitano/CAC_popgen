#!/bin/bash
#SBATCH --job-name=compare_samples_ancestry		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=1		                            # Number of cores per task
#SBATCH --mem=10gb			                                # Total memory for job
#SBATCH --time=00:40:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/myancestry/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/myancestry/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

#prep and submission script for bam_to_GLs_persample.sh

###SETUP

WORKING_DIR=$1
SCRIPTS_DIR=$2
BORICE_FAMS=$3 #three columns: bampath, family, and par/off
CALLTABLE=$4 #columns are samples, rows are sites
SUFFIX=$5 #sample name suffix in header of calltable #_read_1.fastq.gz_mergeChrs.sam.hmm.combined.pass.formatted.posterior
OUTPUT_FILENAME=$6 #will overwrite!
printf 'Script called with inputs: %s\t%s\t%s\t%s\t%s\t%s\n' $1 $2 $3 $4 $5 $6 | tee >(cat >&2)

###SETUP DIRS
mkdir -p $WORKING_DIR
cd $WORKING_DIR



#write header
printf "BORICE_FAM\t" > $OUTPUT_FILENAME
bash ${SCRIPTS_DIR}/compare_samples_ancestry.sh --header >> $OUTPUT_FILENAME

#loop through families
MOM=""
while read -r line; do
	BAMPATH=$(echo "$line" | cut -f1)
	BAMFILE=${BAMPATH##*/}
	PREFIX=${BAMFILE%.*.fds.bam}
	FAM=$(echo "$line" | cut -f2)
	TYPE=$(echo "$line" | cut -f3)
	if [ $TYPE == "par" ]; then
		MOM=$PREFIX
	elif [ $TYPE == "off" ]; then
		printf '%s\t' $FAM >> $OUTPUT_FILENAME
		bash ${SCRIPTS_DIR}/compare_samples_ancestry.sh $CALLTABLE $MOM $PREFIX $SUFFIX >> $OUTPUT_FILENAME
	fi
done < $BORICE_FAMS

printf "Completed all pairs\n"
