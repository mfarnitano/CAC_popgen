#!/bin/bash
#SBATCH --job-name=fastq_align_wrapper		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=5		                            # Number of cores per task
#SBATCH --mem=10gb			                                # Total memory for job
#SBATCH --time=3:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/gutnas_panel/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/gutnas_panel/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

#prep and submission script for fastq_align_each.sh

###SETUP

PROJECT=Tn5_sequencing
BATCH_NAME=gutnas_panel

SCRIPTS_DIR=~/Git_repos/${PROJECT}/${BATCH_NAME}
LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

GENOME=~/Genomes/IM62_v3/Mimulus_guttatus_var_IM62.mainGenome.fasta
FASTQ_TABLE=${WORKING_DIR}/gutnas_panel_fastq_table_big.txt #two columns, directory and prefix

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

###SETUP DIRS
mkdir -p $WORKING_DIR
cd $WORKING_DIR
mkdir -p fastqc
mkdir -p fastq_align_temp
mkdir -p bams

###MODULES
ml SAMtools/1.10-iccifort-2019.5.281
ml BWA/0.7.17-GCC-8.3.0
ml GATK/4.1.6.0-GCCcore-8.3.0-Java-1.8
module list

###prep genome
printf "\n...preparing genome\n" | tee >(cat >&2)
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
	if [ ! -f ${WORKING_DIR}/bams/${INPUT_PREFIX}.ar.fds.bam ]; then
		printf "\n...submitting sample: $INPUT_DIR/$INPUT_PREFIX\n" | tee >(cat >&2)
		sbatch ${SCRIPTS_DIR}/fastq_align_each.sh $PROJECT $BATCH_NAME $GENOME $INPUT_DIR $INPUT_PREFIX
	fi
done < $FASTQ_TABLE

printf "\n...submitted all samples\n" | tee >(cat >&2)
