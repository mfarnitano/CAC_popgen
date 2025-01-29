#!/bin/bash
#SBATCH --job-name=fastq_align_each		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=10		                            # Number of cores per task
#SBATCH --mem=40gb			                                # Total memory for job
#SBATCH --time=48:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/organelles/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/organelles/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

###Creates filtered bam file for single sample from raw fastq, requires input variables
###See submission script fastq_align_wrapper.sh

###SETUP
PROJECT=$1
BATCH_NAME=$2
GENOME=$3
GTAG=$4
INPUT_DIR=$5
INPUT_PREFIX=$6

SCRIPTS_DIR=~/Git_repos/${PROJECT}/${BATCH_NAME}
LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

NTHREADS=10
MEM=40

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

###SETUP DIRS
mkdir -p $WORKING_DIR
cd $WORKING_DIR
mkdir -p ${WORKING_DIR}/fastqc
mkdir -p ${WORKING_DIR}/fastq_align_temp
mkdir -p ${WORKING_DIR}/bams
mkdir -p ${WORKING_DIR}/qualimap
mkdir -p ${WORKING_DIR}/qualimap/${GTAG}


###MODULES
# ml Trimmomatic/0.39-Java-11
# ml SAMtools/1.16.1-GCC-11.3.0
# ml BWA/0.7.17-GCCcore-11.3.0
# ml picard/2.25.1-Java-11
# ml Qualimap/2.2.1-foss-2021b-R-4.1.2
module list

printf "\nstarting analysis for sample ${INPUT_PREFIX}\n" | tee >(cat >&2)

###Trim
printf "\n...trimming with Trimmomatic\n" | tee >(cat >&2)
if [ ! -f ${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}_1P.fq.gz ]; then
	ml Trimmomatic/0.39-Java-11
	java -jar ${EBROOTTRIMMOMATIC}/trimmomatic-0.39.jar \
		PE -threads $NTHREADS \
		${INPUT_DIR}/${INPUT_PREFIX}_read_1.fastq.gz ${INPUT_DIR}/${INPUT_PREFIX}_read_2.fastq.gz \
		-baseout ${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.fq.gz \
		ILLUMINACLIP:${SCRIPTS_DIR}/Nextera_adapters_for_trimming.fa:1:30:15 \
		SLIDINGWINDOW:5:20 \
		MINLEN:30
	module purge
fi

###map to reference and sort
printf "\n...mapping to reference with bwa mem and sorting with samtools\n" | tee >(cat >&2)
ml BWA/0.7.17-GCCcore-11.3.0
ml SAMtools/1.16.1-GCC-11.3.0
bwa mem -t $NTHREADS $GENOME \
	${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}_1P.fq.gz ${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}_2P.fq.gz |
	samtools sort --threads $NTHREADS -O bam -o ${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.${GTAG}.s.bam -

###index bam file
printf "\n...indexing bam file\n" | tee >(cat >&2)
samtools index ${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.${GTAG}.s.bam
module purge

###mark and remove duplicates
printf "\n...marking and removing PCR duplicates with picard\n" | tee >(cat >&2)
ml picard/2.25.1-Java-11
java -jar ${EBROOTPICARD}/picard.jar \
	MarkDuplicates \
	I=${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.${GTAG}.s.bam \
	O=${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.${GTAG}.ds.bam \
	M=${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.${GTAG}.dedupe-metrics.txt \
	REMOVE_DUPLICATES=true
module purge

###filter bam file and index
printf "\n...filtering bamfile for MAPQ>=29, both reads mapped and properly paired, passes platform QC\n" | tee >(cat >&2)
ml SAMtools/1.16.1-GCC-11.3.0
samtools view --threads $NTHREADS -q 29 -f 2 -F 524 -b ${WORKING_DIR}/fastq_align_temp/${INPUT_PREFIX}.${GTAG}.ds.bam > ${WORKING_DIR}/bams/${INPUT_PREFIX}.${GTAG}.fds.bam
samtools index ${WORKING_DIR}/bams/${INPUT_PREFIX}.${GTAG}.fds.bam
module purge

printf "\n...generating coverage summary stats with qualimap\n" | tee >(cat >&2)
ml Qualimap/2.2.1-foss-2021b-R-4.1.2
qualimap bamqc -bam ${WORKING_DIR}/bams/${INPUT_PREFIX}.${GTAG}.fds.bam -c -outdir ${WORKING_DIR}/qualimap/${GTAG}/${INPUT_PREFIX}
module purge

printf "\n...completed sample %s alignment to %s\n" ${INPUT_PREFIX} ${GTAG} | tee >(cat >&2)
