#!/bin/bash
#SBATCH --job-name=genotype_each		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=10		                            # Number of cores per task
#SBATCH --mem=40gb			                                # Total memory for job
#SBATCH --time=96:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/CAC2022/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/CAC2022/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

###Creates gvcf for single sample from filtered bam, requires input variables
###See submission script genotype_wrapper.sh

###SETUP
PROJECT=$1
BATCH_NAME=$2
GENOME=$3
INPUT_PREFIX=$4
SEQ_BATCH=$5
REF_CODE=$6

SCRIPTS_DIR=~/Git_repos/${PROJECT}/${BATCH_NAME}
LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

NTHREADS=10
MEM=40

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

###SETUP DIRS
cd $WORKING_DIR
mkdir -p ${WORKING_DIR}/qualimap
mkdir -p ${WORKING_DIR}/genotyping_temp
mkdir -p ${WORKING_DIR}/gvcfs

###MODULES #fixed for cluster change 09/23
# ml SAMtools/1.16.1-GCC-11.3.0
# ml GATK/4.4.0.0-GCCcore-11.3.0-Java-17
# ml Qualimap/2.2.1-foss-2021b-R-4.1.2
module list

printf "\nstarting analysis for sample ${INPUT_PREFIX}\n" | tee >(cat >&2)

###individual genotyping with GATK
printf "\n...adding read groups\n" | tee >(cat >&2)
ml GATK/4.4.0.0-GCCcore-11.3.0-Java-17
gatk AddOrReplaceReadGroups \
	-I ${WORKING_DIR}/bams/${INPUT_PREFIX}.${REF_CODE}.fds.bam \
	-O $WORKING_DIR/genotyping_temp/${INPUT_PREFIX}.${REF_CODE}.rg.bam \
	-RGID $INPUT_PREFIX \
	-LB ${SEQ_BATCH} \
	-PL illumina \
	-PU $INPUT_PREFIX \
	-SM $INPUT_PREFIX
module purge

printf "\n...indexing bam file\n" | tee >(cat >&2)
ml SAMtools/1.16.1-GCC-11.3.0
samtools index $WORKING_DIR/genotyping_temp/${INPUT_PREFIX}.${REF_CODE}.rg.bam
module purge

printf "\n...genotyping with gatk HaplotypeCaller\n" | tee >(cat >&2)
ml GATK/4.4.0.0-GCCcore-11.3.0-Java-17
gatk --java-options "-Xmx${MEM}g" HaplotypeCaller  \
	-R $GENOME \
	-I $WORKING_DIR/genotyping_temp/${INPUT_PREFIX}.${REF_CODE}.rg.bam \
	-O $WORKING_DIR/gvcfs/${INPUT_PREFIX}.${REF_CODE}.g.vcf.gz \
	-ERC GVCF
module purge
printf "...Finished sample ${INPUT_PREFIX}\n" | tee >(cat >&2)
