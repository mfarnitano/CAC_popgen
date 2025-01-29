#!/bin/bash
#SBATCH --job-name=phylofasta	                # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=4		                            # Number of cores per task
#SBATCH --mem=40gb			                                # Total memory for job
#SBATCH --time=12:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/organelles/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/organelles/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

###SETUP

GTAG=$1 #Mitochondria, Chloroplast
# DPMIN #Use 10 for Chloroplast, 5 for Mitochondria

PROJECT=Tn5_sequencing
BATCH_NAME=organelles

SCRIPTS_DIR=~/Git_repos/${PROJECT}/${BATCH_NAME}
LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

NTHREADS=4
MEM=40

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

###SETUP DIRS
mkdir -p $WORKING_DIR
cd $WORKING_DIR

ml BCFtools/1.15.1-GCC-11.3.0

#PREFIX_PATH=${WORKING_DIR}/VCFs/CAC_DP40_plusOrgPanel_${GTAG}
#PREFIX_PATH=${WORKING_DIR}/VCFs/allCACarea_${GTAG}
PREFIX_PATH=${WORKING_DIR}/VCFs/DPRarea_DP40_plusOrgPanel_${GTAG}

GENOME1=~/Genomes/Organelles/chloroplastIM767_SSC_IRB_LSC.renamed.fa
GENOME2=~/Genomes/Organelles/mitochondriaIM62_excludeChloroRegions.renamed.fa
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

mkdir -p ${WORKING_DIR}/phylo

printf "\n...making fastas\n" | tee >(cat >&2)
VCF_TYPE=${PREFIX_PATH}.merged.allcalled
#VCF_TYPE=${PREFIX_PATH}.merged.allcalled.momsonly
#VCF_TYPE=${PREFIX_PATH}.SNPs.fspi.pm.allcalled.rn
PREFIX=${PREFIX_PATH##*/}
#PREFIX=${PREFIX_PATH##*/}_momsonly

rm -f ${WORKING_DIR}/phylo/${PREFIX}.fa
rm -f ${WORKING_DIR}/phylo/${PREFIX}.full.fa

for sample in $(bcftools query -l ${VCF_TYPE}.vcf.gz); do
	printf "\n>%s\n" $sample >> ${WORKING_DIR}/phylo/${PREFIX}.fa
	bcftools consensus -f $GENOME -a '?' -M 'N' -s $sample ${VCF_TYPE}.vcf.gz |
		grep -v "^>" | tr -d '?' | tr -d '\n' >> ${WORKING_DIR}/phylo/${PREFIX}.fa

	#printf "\n>%s\n" $sample >> ${WORKING_DIR}/phylo/${PREFIX}.full.fa
	#bcftools consensus -f $GENOME -M 'N' -s $sample ${VCF_TYPE}.vcf.gz |
	#	grep -v "^>" | tr -d '\n' >> ${WORKING_DIR}/phylo/${PREFIX}.full.fa
done

#make fasta of only unduplicated sequences, store info about duplication groups
#python ${SCRIPTS_DIR}/remove_duplicated_fastas.py ${WORKING_DIR}/phylo/${PREFIX}.fa ${WORKING_DIR}/phylo/${PREFIX}.unique

#add numbers to front of fasta sample names
#awk 'BEGIN{i=0}{if(substr($1,1,1)==">" && (substr($1,2,4)=="CAC" || substr($1,2,3)=="LM")){i++;print ">" i substr($1,2)}else{print $1}}' ${WORKING_DIR}/phylo/${PREFIX}.fa > ${WORKING_DIR}/phylo/${PREFIX}.renamed.fa
#awk 'BEGIN{i=0}{if(substr($1,1,1)==">"){i++;print ">" i substr($1,2)}else{print $1}}' ${WORKING_DIR}/phylo/${PREFIX}.unique.fa > ${WORKING_DIR}/phylo/${PREFIX}.unique.renamed.fa

#add numbers for only CAC+LM samples
#awk 'BEGIN{i=0}{if(substr($1,1,1)==">" && (substr($1,2,3)=="CAC" || substr($1,2,3)=="LM2" || substr($1,2,4)=="Prog" || substr($1,2,4)=="Moms" || substr($1,2,4)=="Cens" || substr($1,2,2)=="KK")){i++;print ">" i substr($1,2)}else{print $1}}' ${WORKING_DIR}/phylo/${PREFIX}.fa > ${WORKING_DIR}/phylo/${PREFIX}.goodnames.fa
awk 'BEGIN{i=0}{if(substr($1,1,1)==">" && (substr($1,2,3)=="DPR" || substr($1,2,3)=="Tn5" )){i++;print ">" i substr($1,2)}else{print $1}}' ${WORKING_DIR}/phylo/${PREFIX}.fa > ${WORKING_DIR}/phylo/${PREFIX}.goodnames.fa
grep '^>' ${WORKING_DIR}/phylo/${PREFIX}.goodnames.fa | tr -d '>' | cut -c -18 | sed 's/$/\t0,0,0,0,1/' > ${WORKING_DIR}/phylo/${PREFIX}.goodnames.txt

#convert fasta to nexus
ml EMBOSS/6.6.0-intel-2021b
#seqret -sequence ${WORKING_DIR}/phylo/${PREFIX}.renamed.fa -outseq ${WORKING_DIR}/phylo/${PREFIX}.renamed.nexus -osformat nexus
#seqret -sequence ${WORKING_DIR}/phylo/${PREFIX}.unique.renamed.fa -outseq ${WORKING_DIR}/phylo/${PREFIX}.unique.renamed.nexus -osformat nexus
seqret -sequence ${WORKING_DIR}/phylo/${PREFIX}.goodnames.fa -outseq ${WORKING_DIR}/phylo/${PREFIX}.goodnames.nexus -osformat nexus


printf "\n...finished making fastas\n" | tee >(cat >&2)
