#!/bin/bash
#SBATCH --job-name=fastq_align_wrapper		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=5		                            # Number of cores per task
#SBATCH --mem=10gb			                                # Total memory for job
#SBATCH --time=3:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/organelles/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/organelles/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

#prep and submission script for fastq_align_each.sh

###SETUP

PROJECT=Tn5_sequencing
BATCH_NAME=organelles

SCRIPTS_DIR=~/Git_repos/${PROJECT}/${BATCH_NAME}
LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

GENOME1=~/Genomes/Organelles/chloroplastIM767_SSC_IRB_LSC.fa
GENOME2=~/Genomes/Organelles/mitochondriaIM62_excludeChloroRegions.fa
G1TAG=Chloroplast
G2TAG=Mitochondria
FASTQ_TABLE=${WORKING_DIR}/samplelists/everything_fastqtable.txt #two columns, directory and prefix

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

###SETUP DIRS
mkdir -p $WORKING_DIR
cd $WORKING_DIR
mkdir -p fastqc
mkdir -p fastq_align_temp
mkdir -p bams

###MODULES #cluster-updated
ml SAMtools/1.16.1-GCC-11.3.0
ml BWA/0.7.17-GCCcore-11.3.0
ml GATK/4.4.0.0-GCCcore-11.3.0-Java-17
module list

for GENOME in $GENOME1 $GENOME2; do
	if [ $GENOME == $GENOME1 ]; then
		GTAG=$G1TAG
	else
		GTAG=$G2TAG
	fi
	###prep genome
	printf "\n...preparing genome %s \n" $GTAG | tee >(cat >&2)
	if [ ! -f ${GENOME}.fai ]; then
		samtools faidx $GENOME
	fi
	if [ ! -f ${GENOME}.amb ]; then
		bwa index $GENOME
	fi
	if [ ! -f ${GENOME%.fasta}.dict ]; then
		gatk CreateSequenceDictionary -R $GENOME
	fi

	###loop
	while read -r line; do
		INPUT_DIR=$(echo -n $line | tr -s ' ' | cut -d' ' -f1)
		INPUT_PREFIX=$(echo -n $line | tr -s ' ' | cut -d' ' -f2)
		if [ ! -f ${WORKING_DIR}/bams/${INPUT_PREFIX}.${GTAG}.fds.bam ]; then
			printf "\n...submitting sample: $INPUT_DIR/$INPUT_PREFIX\n" | tee >(cat >&2)
			sbatch ${SCRIPTS_DIR}/fastq_align_each.sh $PROJECT $BATCH_NAME $GENOME $GTAG $INPUT_DIR $INPUT_PREFIX
		fi
	done < $FASTQ_TABLE

	printf "\n...submitted all samples for %s\n" $GTAG | tee >(cat >&2)
done
