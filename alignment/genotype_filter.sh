#!/bin/bash
#SBATCH --job-name=genotype_filter	                # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=32		                            # Number of cores per task
#SBATCH --mem=120gb			                                # Total memory for job
#SBATCH --time=120:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/CAC2022/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/CAC2022/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

###SETUP

PROJECT=Tn5_sequencing
BATCH_NAME=CAC2022

SCRIPTS_DIR=~/Git_repos/${PROJECT}/${BATCH_NAME}
LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

GENOME=~/Genomes/IM62_v3/Mimulus_guttatus_var_IM62.mainGenome.fixed.fasta
REPEATMASK=~/Genomes/IM62_v3/Mimulus_guttatus_var_IM62.mainGenome.fixed.repeatMasked.gff

NTHREADS=32
MEM=120

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

###SETUP DIRS
mkdir -p $WORKING_DIR
cd $WORKING_DIR
mkdir -p VCFs

###MODULES
ml GATK/4.4.0.0-GCCcore-11.3.0-Java-17
module list

###create ref dict
printf "\n...creating chromosome list\n" | tee >(cat >&2)
if [ ! -f ${GENOME%.fasta}.dict ]; then
	gatk CreateSequenceDictionary -R ${GENOME}
fi

###create chr_list
printf "\n...creating chromosome list\n" | tee >(cat >&2)
if [ ! -f ${WORKING_DIR}/chr_positions.list ]; then
	head -n14 ${GENOME}.fai | awk '{print $1 ":1-"$2}' > ${WORKING_DIR}/chr_positions.list
fi


printf "\n...merging 14 chromosomes\n" | tee >(cat >&2)
if [ ! -f ${WORKING_DIR}/VCFs/Census2022.allChrs.allsites.vcf.gz ]; then
	gatk --java-options "-Xmx${MEM}G" MergeVcfs \
		I=${WORKING_DIR}/VCFs/Census2022.Chr01.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/Census2022.Chr02.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/Census2022.Chr03.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/Census2022.Chr04.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/Census2022.Chr05.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/Census2022.Chr06.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/Census2022.Chr07.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/Census2022.Chr08.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/Census2022.Chr09.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/Census2022.Chr10.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/Census2022.Chr11.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/Census2022.Chr12.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/Census2022.Chr13.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/Census2022.Chr14.allsites.vcf.gz \
		O=${WORKING_DIR}/VCFs/Census2022.allChrs.allsites.vcf.gz
fi

#inputs
RAW_VCF=${WORKING_DIR}/VCFs/Census2022.allChrs.allsites.vcf.gz
GROUP=Census2022.allChrs
REFCODE=IM62_v3
PREFIX=${WORKING_DIR}/VCFs/${GROUP}.${REFCODE}


printf "\n...extracting biallelic SNPs\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SelectVariants -V ${RAW_VCF} \
	-select-type SNP --restrict-alleles-to BIALLELIC -O ${PREFIX}.SNPs.vcf.gz

printf "\n...extracting invariant sites\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SelectVariants -V ${RAW_VCF} \
	-select-type NO_VARIATION -O ${PREFIX}.INVTs.vcf.gz

printf "\n...filtering SNPs\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" VariantFiltration -V ${PREFIX}.SNPs.vcf.gz -O ${PREFIX}.SNPs.f.vcf.gz \
	-filter "QD < 2.0" --filter-name "QD2" \
	-filter "QUAL < 40.0" --filter-name "QUAL40" \
	-filter "SOR > 3.0" --filter-name "SOR4" \
	-filter "FS > 60.0" --filter-name "FS60" \
	-filter "MQ < 40.0" --filter-name "MQ40" \
	-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
	-filter "ReadPosRankSum < -12.5" --filter-name "ReadPosRankSum-12.5" \
	-filter "ReadPosRankSum > 12.5" --filter-name "ReadPosRankSum12.5" \
	--verbosity ERROR

printf "\n...filtering INVTs\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}g" VariantFiltration -V ${PREFIX}.INVTs.vcf.gz -O ${PREFIX}.INVTs.f.vcf.gz \
	-filter "QD < 2.0" --filter-name "QD2" \
	-filter "SOR > 3.0" --filter-name "SOR4" \
	-filter "MQ < 40.0" --filter-name "MQ40" \
	--verbosity ERROR

printf "\n...sorting SNPs\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SortVcf \
	-I ${PREFIX}.SNPs.f.vcf.gz -SD ${GENOME%.fasta}.dict -O ${PREFIX}.SNPs.fs.vcf.gz

printf "\n...sorting INVTs\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SortVcf \
	-I ${PREFIX}.INVTs.f.vcf.gz -SD ${GENOME%.fasta}.dict -O ${PREFIX}.INVTs.fs.vcf.gz

printf "\n...selecting passing SNPs\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SelectVariants \
	-V ${PREFIX}.SNPs.fs.vcf.gz --exclude-filtered -O ${PREFIX}.SNPs.fsp.vcf.gz

printf "\n...selecting passing INVTs\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SelectVariants \
	-V ${PREFIX}.INVTs.fs.vcf.gz --exclude-filtered -O ${PREFIX}.INVTs.fsp.vcf.gz

module purge
ml VCFtools/0.1.16-GCC-11.2.0
ml HTSlib/1.15.1-GCC-11.3.0

printf "\n...prepare bed file of repeat-masked regions\n" | tee >(cat >&2)
printf '#chrom\tchromStart\tchromEnd\n' > ${REPEATMASK}.bed
cut -f1,4,5 ${REPEATMASK} >> ${REPEATMASK}.bed

printf "\n...filter individual SNP genotypes by depth and GQ, and exclude repeat-masked regions\n" | tee >(cat >&2)
vcftools --gzvcf ${PREFIX}.SNPs.fsp.vcf.gz -c --minGQ 20 --minDP 2 --maxDP 20 \
	--exclude-bed ${REPEATMASK}.bed \
	--recode --recode-INFO-all | bgzip -c > ${PREFIX}.SNPs.fspi.rm.vcf.gz

printf "\n...filter individual INVT genotypes by depth\n" | tee >(cat >&2)
vcftools --gzvcf ${PREFIX}.INVTs.fsp.vcf.gz -c --minDP 2 --maxDP 20 \
	--exclude-bed ${REPEATMASK}.bed \
	--recode --recode-INFO-all | bgzip -c > ${PREFIX}.INVTs.fspi.rm.vcf.gz


tabix -p vcf ${PREFIX}.SNPs.fspi.rm.vcf.gz
tabix -p vcf ${PREFIX}.INVTs.fspi.rm.vcf.gz

module purge
ml GATK/4.4.0.0-GCCcore-11.3.0-Java-17

printf "\n...merge SNP and INVT sites\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" MergeVcfs \
	-I ${PREFIX}.SNPs.fspi.rm.vcf.gz -I ${PREFIX}.INVTs.fspi.rm.vcf.gz \
	-O ${PREFIX}.merged.fpi.rm.vcf.gz

printf "\n...sort merged VCF\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SortVcf \
	-I ${PREFIX}.merged.fpi.rm.vcf.gz \
	-SD ${GENOME%.fasta}.dict -O ${PREFIX}.merged.fspi.rm.vcf.gz

printf "\n...filter SNPs by mincalled\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SelectVariants -V ${PREFIX}.SNPs.fspi.rm.vcf.gz \
	--max-nocall-number 31 --exclude-filtered -O ${PREFIX}.SNPs.fspi.rm.31called.vcf.gz

printf "\n...filter merged VCF by mincalled\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SelectVariants -V ${PREFIX}.merged.fspi.rm.vcf.gz \
	--max-nocall-number 31 --exclude-filtered -O ${PREFIX}.merged.fspi.rm.31called.vcf.gz

printf "\n...completing all filtering steps\n" | tee >(cat >&2)
