#!/bin/bash
#SBATCH --job-name=angsd_GLs                     # Job name
#SBATCH --partition=highmem_p	                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=9		                            # Number of cores per task
#SBATCH --mem=480gb			                                # Total memory for job
#SBATCH --time=96:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/borice_hq/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/borice_hq/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

###SETUP
WORKING_DIR=$1
GROUPID=$2
BAMLIST=$3 #single column, full path to bams
SNPLIST=$4 #two columns, chr and pos
SNPLIST_TAG=$5
CHR=$6 #allChrs if blank

printf "Script called with inputs %s %s %s %s %s %s \n" $1 $2 $3 $4 $5 $6 | tee >(cat >&2)


if [ -z $CHR ]; then
	CHR=allChrs
fi
printf "CHR is %s\n" $CHR | tee >(cat >&2)

SCRIPTS_DIR=~/Git_repos/Tn5_sequencing/angsd_structure
LOG_DIR=/home/mcf96392/logs/borice_hq

NTHREADS=8
MEM=480

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt


GENOME=~/Genomes/IM62_v3/Mimulus_guttatus_var_IM62.mainGenome.fixed.fasta
GTAG=IM62_v3

###SETUP DIRS

mkdir -p ${WORKING_DIR}
cd ${WORKING_DIR}

###MODULES
ml SAMtools/1.14-GCC-11.2.0
ml angsd/0.940-GCC-11.2.0


#make hq SNPs list
# if [ $CHR == "allChrs" ]; then
# 	cut -f1,2 $SNPLIST > ${SNPLIST_TAG}.${CHR}.pos.txt
# 	rm -f ${SNPLIST_TAG}.${CHR}.chrs.txt
# 	for i in Chr{01..14}; do printf '%s\n' $i >> ${SNPLIST_TAG}.${CHR}.chrs.txt; done
# else
# 	grep $CHR $SNPLIST | cut -f1,2 > ${SNPLIST_TAG}.${CHR}.pos.txt
# 	echo $CHR > ${SNPLIST_TAG}.${CHR}.chrs.txt
# fi

printf "\n...starting genotyping for %s\n" ${SNPLIST_TAG}.${CHR} | tee >(cat >&2)
###run angsd to get genotype likelihoods
if [ ! -f ${WORKING_DIR}/genolike.${GROUPID}.${GTAG}.${SNPLIST_TAG}.${CHR}.beagle.gz ]; then
	# angsd sites index ${SNPLIST_TAG}.${CHR}.pos.txt
	angsd -GL 2 -out genolike.${GROUPID}.${GTAG}.IM62_v3.${SNPLIST_TAG}.${CHR} -nThreads $NTHREADS \
		-rf ${SNPLIST_TAG}.${CHR}.chrs.txt -sites ${SNPLIST_TAG}.${CHR}.pos.txt \
		-doGlf 2 -doCounts 1 -doMajorMinor 1 -doMaf 1 -bam ${BAMLIST}

fi

printf "\n...finished genotyping for %s\n" ${SNPLIST_TAG}.$CHR | tee >(cat >&2)
