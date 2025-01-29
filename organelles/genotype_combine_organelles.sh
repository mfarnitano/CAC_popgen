#!/bin/bash
#SBATCH --job-name=genotype_combine                     # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=32		                            # Number of cores per task
#SBATCH --mem=80gb			                                # Total memory for job
#SBATCH --time=48:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/organelles/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/organelles/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

###SETUP
GTAG=$1


PROJECT=Tn5_sequencing
BATCH_NAME=organelles

SCRIPTS_DIR=~/Git_repos/${PROJECT}/${BATCH_NAME}
LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

NTHREADS=32
MEM=80

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


GENOME1=~/Genomes/Organelles/chloroplastIM767_SSC_IRB_LSC.fa
GENOME2=~/Genomes/Organelles/mitochondriaIM62_excludeChloroRegions.fa
G1TAG=Chloroplast
G2TAG=Mitochondria

if [ $GTAG == $G1TAG ]; then
	GENOME=$GENOME1
elif [ $GTAG == $G2TAG ]; then
	GENOME=$GENOME2
else
	printf "\n...genome tag not specified, try %s or %s\n" $G1TAG $G2TAG | tee >(cat >&2)
	exit 1
fi

#make interval list
#cp ${GENOME%.fa}.dict ${GENOME%.fa}.interval_list
#awk 'BEGIN{OFS="\t"}{print $1,"1",$2,"+",$1}' ${GENOME}.fai >> ${GENOME%.fa}.interval_list
#printf "Intervals: \n%s\n" $INTERVALS | tee >(cat >&2)

##merge gvcfs into database
#printf "\n...merging gvcfs into database, selecting only major (chromosome) linkage groups\n" | tee >(cat >&2)
#cd ${WORKING_DIR}/gvcfs
#gatk --java-options "-Xmx${MEM}g" GenomicsDBImport \
#	--sample-name-map ${WORKING_DIR}/Everything_${GTAG}_gvcf_map.txt \
#	--genomicsdb-workspace-path ${WORKING_DIR}/Everything_${GTAG}_vcf_db \
#	-L ${GENOME%.fa}.interval_list
#cd ${WORKING_DIR}

###joint genotyping
printf "\n...Genotyping VCFs...outputing all-sites VCF with min call confidence = 30" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" GenotypeGVCFs \
	-R $GENOME \
	-V gendb://${WORKING_DIR}/Everything_${GTAG}_vcf_db \
	-stand-call-conf 30 \
	-all-sites \
	-O ${WORKING_DIR}/VCFs/Everything_${GTAG}.allsites.vcf.gz

printf "\n...Finished genotyping\n" | tee >(cat >&2)
