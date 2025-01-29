#!/bin/bash
#SBATCH --job-name=genotype_wrapper		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=5		                            # Number of cores per task
#SBATCH --mem=10gb			                                # Total memory for job
#SBATCH --time=3:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/CAC2022/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/CAC2022/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

#prep and submission script for genotype_each.sh

###SETUP

PROJECT=Tn5_sequencing
BATCH_NAME=CAC2022
GROUP=Census2022

SCRIPTS_DIR=~/Git_repos/${PROJECT}/${BATCH_NAME}
LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

GENOME=~/Genomes/IM62_v3/Mimulus_guttatus_var_IM62.mainGenome.fixed.fasta
BAM_TABLE=${WORKING_DIR}/${GROUP}_prefixes.txt #two columns, directory and prefix
ls ${WORKING_DIR}/bams | grep 'IM62_v3.fds.bam$' | grep Census2022 | sed 's/.IM62_v3.fds.bam//' > $BAM_TABLE

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

###SETUP DIRS
mkdir -p $WORKING_DIR
cd $WORKING_DIR
mkdir -p genotyping_temp
mkdir -p gvcfs

###MODULES #fixed for cluster change 09/23
###need to load individually
# ml SAMtools/1.16.1-GCC-11.3.0
# ml BWA/0.7.17-GCCcore-11.3.0
# ml GATK/4.4.0.0-GCCcore-11.3.0-Java-17
# module list

###prep genome
printf "\n...preparing genome\n" | tee >(cat >&2)
if [ ! -f ${GENOME}.fai ]; then
	ml SAMtools/1.16.1-GCC-11.3.0
	samtools faidx $GENOME
	module purge
fi
if [ ! -f ${GENOME}.amb ]; then
	ml BWA/0.7.17-GCCcore-11.3.0
	bwa index $GENOME
	module purge
fi
if [ ! -f ${GENOME%.fasta}.dict ]; then
	ml GATK/4.4.0.0-GCCcore-11.3.0-Java-17
	gatk CreateSequenceDictionary -R $GENOME
	module purge
fi

if [ -f ${WORKING_DIR}/Census2022_gvcf_map.txt ]; then
	rm ${WORKING_DIR}/Census2022_gvcf_map.txt
fi

REF_CODE=IM62_v3

###loop
while read -r line; do
	INPUT_PREFIX=${line}
	SEQBATCH="Census2022"
	printf "\n...submitting sample: $INPUT_PREFIX\n" | tee >(cat >&2)
	if [ ! -f $WORKING_DIR/gvcfs/${INPUT_PREFIX}.${REF_CODE}.g.vcf.gz ]; then
		sbatch ${SCRIPTS_DIR}/genotype_each.sh $PROJECT $BATCH_NAME $GENOME $INPUT_PREFIX $SEQBATCH $REF_CODE
	fi
	printf "%s\t%s.%s.g.vcf.gz\n" $INPUT_PREFIX $INPUT_PREFIX $REF_CODE >> ${WORKING_DIR}/${GROUP}_gvcf_map.txt
done < $BAM_TABLE

printf "\n...submitted all samples\n" | tee >(cat >&2)
