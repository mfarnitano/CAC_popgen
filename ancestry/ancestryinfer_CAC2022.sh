#!/bin/bash
#SBATCH --job-name=ancestry_v3	                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=10		                            # Number of cores per task
#SBATCH --mem=20gb			                                # Total memory for job
#SBATCH --time=12:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/CAC2022/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/CAC2022/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

PROJECT=Tn5_sequencing
BATCH_NAME=CAC2022

SCRIPTS_DIR=~/Git_repos/${PROJECT}/${BATCH_NAME}
LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

NTHREADS=10
MEM=20

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

cd ${WORKING_DIR}
mkdir -p ${WORKING_DIR}/ancestry_v3
cd ${WORKING_DIR}/ancestry_v3

#make readlist


###READ NAMES MUST BE IN FORMAT ${PREFIX}_read_1.fastq.gz ${PREFIX}_read_2.fastq.gz
###GENOME chr names cannot have '_' in them!!
###do FASTQs need to be in working dir??

FASTQ_TABLE=${WORKING_DIR}/ancestry_v3/${BATCH_NAME}_fastqtable.txt

READLIST=${WORKING_DIR}/ancestry_v3/${BATCH_NAME}_readlist.txt
rm -f $READLIST
while read -r line; do
	folder=$(echo -n $line | tr -s ' ' | cut -d' ' -f1)
	prefix=$(echo -n $line | tr -s ' ' | cut -d' ' -f2)
	#cp ${folder}/${prefix}_read_1.fastq.gz ./${prefix}_read_1.fastq.gz
	#cp ${folder}/${prefix}_read_2.fastq.gz ./${prefix}_read_2.fastq.gz
	#printf '%s/%s_read_1.fastq.gz\t%s/%s_read_2.fastq.gz\n' $folder $prefix $folder $prefix >> $READLIST
	printf './%s_read_1.fastq.gz\t./%s_read_2.fastq.gz\n' $prefix $prefix >> $READLIST
done < ${FASTQ_TABLE}


###did this manually
### conda activate ~/ancestryinfer/ngsutil
### conda install -c bioconda samtools openssl=1.0

### conda activate py2
### conda install -c bioconda samtools openssl=1.0

source ~/ancestryinfer/sourceme
ml Miniconda3/4.7.10
source activate ~/ancestryinfer/ngsutil
#source activate py2
python --version
conda info
samtools --version

###Run ancestryinfer pipeline
perl ~/ancestryinfer/ancestryinfer/Ancestry_HMM_parallel_v6.pl ${SCRIPTS_DIR}/ancestryinfer_${BATCH_NAME}.cfg
