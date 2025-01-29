#!/bin/bash
#SBATCH --job-name=genotype_filter	                # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=32		                            # Number of cores per task
#SBATCH --mem=120gb			                                # Total memory for job
#SBATCH --time=120:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/gutnas_panel/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/gutnas_panel/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

###SETUP

PROJECT=Tn5_sequencing
BATCH_NAME=gutnas_panel

SCRIPTS_DIR=~/Git_repos/${PROJECT}/${BATCH_NAME}
LOG_DIR=/home/mcf96392/logs/${BATCH_NAME}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

GENOME=~/Genomes/IM62_v3/Mimulus_guttatus_var_IM62.mainGenome.fasta
REPEATMASK=~/Genomes/IM62_v3/Mimulus_guttatus_var_IM62.mainGenome.repeatMasked.gff
FASTQ_TABLE=${WORKING_DIR}/gutnas_panel_fastq_table_bigger.txt #two columns, directory and prefix

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

printf "\n...merging 14 chromosomes\n" | tee >(cat >&2)
if [ ! -f ${WORKING_DIR}/VCFs/gutnas_panel.allChrs.allsites.vcf.gz ]; then
	gatk --java-options "-Xmx${MEM}G" MergeVcfs \
		I=${WORKING_DIR}/VCFs/gutnas_panel.Chr_01.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/gutnas_panel.Chr_02.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/gutnas_panel.Chr_03.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/gutnas_panel.Chr_04.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/gutnas_panel.Chr_05.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/gutnas_panel.Chr_06.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/gutnas_panel.Chr_07.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/gutnas_panel.Chr_08.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/gutnas_panel.Chr_09.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/gutnas_panel.Chr_10.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/gutnas_panel.Chr_11.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/gutnas_panel.Chr_12.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/gutnas_panel.Chr_13.allsites.vcf.gz \
		I=${WORKING_DIR}/VCFs/gutnas_panel.Chr_14.allsites.vcf.gz \
		O=${WORKING_DIR}/VCFs/gutnas_panel.allChrs.allsites.vcf.gz
fi

#inputs
RAW_VCF=${WORKING_DIR}/VCFs/gutnas_panel.allChrs.allsites.vcf.gz
GROUP=panel38
REFCODE=IM62_v3
PREFIX=${WORKING_DIR}/VCFs/${GROUP}.${REFCODE}


###create chr_list
printf "\n...creating chromosome list\n" | tee >(cat >&2)
if [ ! -f ${WORKING_DIR}/chr_positions.list ]; then
	head -n14 ${GENOME}.fai | awk '{print $1 ":1-"$2}' > ${WORKING_DIR}/chr_positions.list
fi

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
vcftools --gzvcf ${PREFIX}.SNPs.fsp.vcf.gz -c --minGQ 15 --minDP 6 --maxDP 100 \
	--exclude-bed ${REPEATMASK}.bed \
	--recode --recode-INFO-all | bgzip -c > ${PREFIX}.SNPs.fspi.rm.vcf.gz

printf "\n...filter individual INVT genotypes by depth\n" | tee >(cat >&2)
vcftools --gzvcf ${PREFIX}.INVTs.fsp.vcf.gz -c --minDP 6 --maxDP 100 \
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
	--max-nocall-number 7 --exclude-filtered -O ${PREFIX}.SNPs.fspi.rm.31called.vcf.gz

printf "\n...filter merged VCF by mincalled\n" | tee >(cat >&2)
gatk --java-options "-Xmx${MEM}G" SelectVariants -V ${PREFIX}.merged.fspi.rm.vcf.gz \
	--max-nocall-number 7 --exclude-filtered -O ${PREFIX}.merged.fspi.rm.31called.vcf.gz

printf "\n...completing all filtering steps\n" | tee >(cat >&2)

printf "\n...making genotype tables\n" | tee >(cat >&2)
module purge
ml BCFtools/1.15.1-GCC-11.3.0

printf 'CHROM\tPOS\tREF\tALT\t' > ${PREFIX}.SNPs.fspi.rm.table
printf 'CHROM\tPOS\tREF\tALT\t' > ${PREFIX}.merged.fspi.rm.table

bcftools query -l ${PREFIX}.SNPs.fspi.rm.vcf.gz | tr '\n' '\t' >> ${PREFIX}.SNPs.fspi.rm.table
bcftools query -l ${PREFIX}.merged.fspi.rm.vcf.gz | tr '\n' '\t' >> ${PREFIX}.merged.fspi.rm.table

printf '\n' >> ${PREFIX}.SNPs.fspi.rm.table
printf '\n' >> ${PREFIX}.merged.fspi.rm.table

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' ${PREFIX}.SNPs.fspi.rm.vcf.gz >> ${PREFIX}.SNPs.fspi.rm.table
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' ${PREFIX}.merged.fspi.rm.vcf.gz >> ${PREFIX}.merged.fspi.rm.table

cat ${PREFIX}.SNPs.fspi.rm.table | sed 's#0/0#0#g' | sed 's#0|0#0#g' | sed 's#0/1#1#g' | sed 's#1/0#1#g' | sed 's#0|1#1#g' | sed 's#1|0#1#g' | sed 's#1/1#2#g' | sed 's#1|1#2#g' | sed 's#./.#NA#g' > ${PREFIX}.SNPs.fspi.rm.gt.table
cat ${PREFIX}.merged.fspi.rm.table | sed 's#0/0#0#g' | sed 's#0|0#0#g' | sed 's#0/1#1#g' | sed 's#1/0#1#g' | sed 's#0|1#1#g' | sed 's#1|0#1#g' | sed 's#1/1#2#g' | sed 's#1|1#2#g' | sed 's#./.#NA#g' > ${PREFIX}.merged.fspi.rm.gt.table

###other steps
python ${SCRIPTS_DIR}/genocounts_groups.py ${WORKING_DIR}/VCFs/gutnas_panel_bigger_popfile.txt ${PREFIX}.merged.fspi.rm.gt.table > ${PREFIX}.merged.fspi.rm.counts.txt
awk '$6<3 && $12<2 && $15<3 && $18<3 {print $1,$2,$3,$4}' ${PREFIX}.merged.fspi.rm.counts.txt > ${PREFIX}.merged.fspi.rm.highqual.sites

python ${SCRIPTS_DIR}/genocounts_groups.py ${WORKING_DIR}/VCFs/gutnas_panel_bigger_popfile.txt ${PREFIX}.SNPs.fspi.rm.gt.table > ${PREFIX}.SNPs.fspi.rm.counts.txt
awk '$6<3 && $12<2 && $15<3 && $18<3 {print}' panel38.IM62_v3.SNPs.fspi.rm.counts.txt > panel38.IM62_v3.SNPs.fspi.rm.highqual.sites.counts.txt
awk '{gutcount=$5+$14; guttotal=$7-$6+$16-$15; gutref=guttotal*2-gutcount; gutprop=gutcount/(2*guttotal); nascount=$11; nastotal=$13-$12; nasprop=nascount/(2*nastotal); nasref=nastotal*2-nascount; if (nasprop-gutprop>=0.8) {print $1,$2,$3,$4,gutref,gutcount,nasref,nascount}}' ${PREFIX}.SNPs.fspi.rm.highqual.sites.counts.txt > ${PREFIX}.SNPs.fspi.rm.AIMs_counts.txt

cut -f1,2,3,4 ${PREFIX}.SNPs.fspi.rm.AIMs_counts.txt > /scratch/mcf96392/AIMs/AIMs_panel38_final.AIMs.txt
cut -f1,2,5,6,7,8 ${PREFIX}.SNPs.fspi.rm.AIMs_counts.txt > /scratch/mcf96392/AIMs/AIMs_panel38_final.AIMs_counts.txt
cut -f1,2,3,4 ${PREFIX}.SNPs.fspi.rm.highqual.sites.counts.txt > /scratch/mcf96392/AIMs/SNPs_panel38_final.txt
