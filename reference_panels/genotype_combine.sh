#!/bin/bash
#SBATCH --job-name=genotype_combine                     # Job name
#SBATCH --partition=highmem_p		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=32		                            # Number of cores per task
#SBATCH --mem=240gb			                                # Total memory for job
#SBATCH --time=168:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/gutnas_panel/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/gutnas_panel/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

###SETUP
CHR=$1

PROJECT=Tn5_sequencing
BATCH_NAME=gutnas_panel

SCRIPTS_DIR=~/Git_repos/${PROJECT}/${BATCH_NAME}
LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

GENOME=~/Genomes/IM62_v3/Mimulus_guttatus_var_IM62.mainGenome.fasta
FASTQ_TABLE=${WORKING_DIR}/gutnas_panel_fastq_table_bigger.txt #two columns, directory and prefix

NTHREADS=32
MEM=240

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

###SETUP DIRS
mkdir -p $WORKING_DIR
cd $WORKING_DIR
mkdir -p VCFs

###MODULES #fixed for 09/23 cluster update
ml GATK/4.4.0.0-GCCcore-11.3.0-Java-17
module list

###create chr_list
printf "\n...creating chromosome list\n" | tee >(cat >&2)


##merge gvcfs into database
printf "\n...merging gvcfs into database, selecting only major (chromosome) linkage groups\n" | tee >(cat >&2)
cd ${WORKING_DIR}/gvcfs
gatk --java-options "-Xmx${MEM}g" GenomicsDBImport \
	--sample-name-map ${WORKING_DIR}/gvcf_map.txt \
	--genomicsdb-workspace-path ${WORKING_DIR}/${BATCH_NAME}_${CHR}_vcf_db \
	-L $CHR
cd $WORKING_DIR

###joint genotyping
printf "\n...Genotyping VCFs...outputing all-sites VCF with min call confidence = 30" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" GenotypeGVCFs \
	-R $GENOME \
	-L $CHR \
	-V gendb://${WORKING_DIR}/${BATCH_NAME}_${CHR}_vcf_db \
	-stand-call-conf 30 \
	-all-sites \
	-O ${WORKING_DIR}/VCFs/${BATCH_NAME}.${CHR}.allsites.vcf.gz

printf "\n...Finished genotyping\n" | tee >(cat >&2)
