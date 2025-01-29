#!/bin/bash
#SBATCH --job-name=genotype_filter	                # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=32		                            # Number of cores per task
#SBATCH --mem=80gb			                                # Total memory for job
#SBATCH --time=24:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/organelles/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/organelles/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

###SETUP

GTAG=$1 #Mitochondria, Chloroplast
# DPMIN #Use 40

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
ml BCFtools/1.15.1-GCC-11.3.0
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

RAW_VCF=${WORKING_DIR}/VCFs/Everything_${GTAG}.allsites.vcf.gz
#PREFIX=${WORKING_DIR}/VCFs/CAC_DP40_plusOrgPanel_${GTAG}
PREFIX=${WORKING_DIR}/VCFs/DPRarea_DP40_plusOrgPanel_${GTAG}

PANEL_VCF=${WORKING_DIR}/VCFs/OrgPanel_${GTAG}.vcf.gz
#SAMPLE_LIST=${WORKING_DIR}/qualimap/${GTAG}_DP40_CAConly.txt
SAMPLE_LIST=${WORKING_DIR}/qualimap/${GTAG}_DP40_DPRonly.txt

###fix scaffold names
if [ $GTAG == $G1TAG ]; then
	#Chloroplast
	for VCF in $RAW_VCF $PANEL_VCF; do
		zcat $VCF |
			sed 's/scaffold_124:1343-85612/scaffold_124_1343_85612/g' |
			sed 's/scaffold_298:13-24197/scaffold_298_13_24197/g' |
			sed 's/scaffold_299:12-946/scaffold_299_12_946/g' |
			sed 's/scaffold_299:18885-19772/scaffold_299_18885_19772/g' |
			sed 's/scaffold_299:963-18861/scaffold_299_963_18861/g' |
			bgzip -c > ${VCF%.vcf.gz}.rn.vcf.gz
			tabix -p vcf ${VCF%.vcf.gz}.rn.vcf.gz
	done
else
	#Mitochondria
	for VCF in $RAW_VCF $PANEL_VCF; do
		zcat $VCF > ${VCF%.vcf.gz}.rn.vcf
		cat ~/Genomes/Organelles/Mito_renamer_key.txt | while read -r line; do
			OLD=$(echo -n $line | tr -s ' ' | cut -d' ' -f1)
			NEW=$(echo -n $line | tr -s ' ' | cut -d' ' -f2)
			sed -i "s/${OLD}/${NEW}/" ${VCF%.vcf.gz}.rn.vcf
		done
		bgzip ${VCF%.vcf.gz}.rn.vcf
		tabix -p vcf ${VCF%.vcf.gz}.rn.vcf.gz
	done
fi

printf "\n...removing samples with low average DP, setting DP<4 or hets to nocall, sorting\n"| tee >(cat >&2)
bcftools view -S $SAMPLE_LIST -Ou ${RAW_VCF%.vcf.gz}.rn.vcf.gz |
	bcftools sort -Ov -o ${PREFIX}.pre.vcf.gz
bcftools index -t ${PREFIX}.pre.vcf.gz

#get snplist from panel file
bcftools query -f "%CHROM\t%POS\n" ${PANEL_VCF%.vcf.gz}.rn.vcf.gz > ${GTAG}.panel.SNPlist.txt

#merge and keep snps from panel
bcftools merge -R ${GTAG}.panel.SNPlist.txt -Ou \
	${PANEL_VCF%.vcf.gz}.rn.vcf.gz ${PREFIX}.pre.vcf.gz |
	bcftools view -V indels,mnps -M 2 -Ou |
	bcftools filter -e 'GT=="het" | FORMAT/DP < 4' -Ov -o ${PREFIX}.merged.vcf.gz

bcftools filter -e 'N_MISSING > 0' -Ov -o ${PREFIX}.merged.allcalled.vcf.gz ${PREFIX}.merged.vcf.gz

printf "\n...completing all filtering steps\n" | tee >(cat >&2)

###MANUAL
printf "CHROM\tPOS\t" > ${PREFIX}.merged.table
bcftools query -l ${PREFIX}.merged.vcf.gz | tr '\n' '\t' >> ${PREFIX}.merged.table
printf "\n" >> ${PREFIX}.merged.table
bcftools query -f "%CHROM\t%POS[\t%GT]\n" ${PREFIX}.merged.vcf.gz | tr '|' '/' >> ${PREFIX}.merged.table

#
# printf "\n...making genotype tables\n" | tee >(cat >&2)
# module purge
# ml BCFtools
#
# printf 'CHROM\tPOS\t' > ${PREFIX}.SNPs.fspi.rm.table
#
# bcftools query -l ${PREFIX}.SNPs.fspi.rm.vcf.gz | tr '\n' '\t' >> ${PREFIX}.SNPs.fspi.rm.table
#
# printf '\n' >> ${PREFIX}.SNPs.fspi.rm.table
#
# bcftools query -f '%CHROM\t%POS[\t%GT]\n' ${PREFIX}.SNPs.fspi.rm.vcf.gz >> ${PREFIX}.SNPs.fspi.rm.table
#
# cat ${PREFIX}.SNPs.fspi.rm.table | sed 's#0/0#0#g' | sed 's#0|0#0#g' | sed 's#0/1#1#g' | sed 's#1/0#1#g' | sed 's#0|1#1#g' | sed 's#1|0#1#g' | sed 's#1/1#2#g' | sed 's#1|1#2#g' | sed 's#./.#NA#g' > ${PREFIX}.SNPs.fspi.rm.gt.table

###other steps
# python genocounts_groups.py gutnas_panel_bigger_popfile.txt panel38.IM62_v3.merged.fspi.rm.gt.table > panel38.IM62_v3.merged.fspi.rm.counts.txt
# awk '$4<3 && $10<2 && $13<3 && $16<3 {print $1,$2}' panel38.IM62_v3.merged.fspi.rm.counts.txt > panel38.IM62_v3.merged.fspi.rm.highqual.sites
