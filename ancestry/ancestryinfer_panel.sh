#!/bin/bash
#SBATCH --job-name=ancestryinfer	                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=10		                            # Number of cores per task
#SBATCH --mem=20gb			                                # Total memory for job
#SBATCH --time=48:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/ancestry_panel/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/ancestry_panel/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

PROJECT=Tn5_sequencing
BATCH_NAME=ancestry_panel

SCRIPTS_DIR=~/Git_repos/${PROJECT}/${BATCH_NAME}
LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

NTHREADS=10
MEM=20

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

mkdir -p ${WORKING_DIR}
cd ${WORKING_DIR}


#make readlist


###READ NAMES MUST BE IN FORMAT ${PREFIX}_read_1.fastq.gz ${PREFIX}_read_2.fastq.gz
###GENOME chr names cannot have '_' in them!!
###downsample FASTQs into working dir
ml seqtk/1.3-GCC-11.3.0

FASTQ_TABLE=${WORKING_DIR}/${BATCH_NAME}_fastqtable.txt

READLIST=${WORKING_DIR}/${BATCH_NAME}_readlist.txt
rm -f $READLIST
while read -r line; do
	folder=$(echo -n $line | tr -s ' ' | cut -d' ' -f1)
	prefix=$(echo -n $line | tr -s ' ' | cut -d' ' -f2)
	seqtk sample -s100 ${folder}/${prefix}_read_1.fastq.gz 2000000 > ./${prefix}_sampled_read_1.fastq.gz
	seqtk sample -s100 ${folder}/${prefix}_read_2.fastq.gz 2000000 > ./${prefix}_sampled_read_2.fastq.gz
	printf '%s_sampled_read_1.fastq.gz\t%s_sampled_read_2.fastq.gz\n' $prefix $prefix >> $READLIST
done < ${FASTQ_TABLE}

source ~/ancestryinferv2/sourceme

###Run ancestryinfer pipeline
Ancestry_HMM_parallel_v7.pl ${SCRIPTS_DIR}/ancestryinfer_panel.cfg
