#!/bin/bash
#SBATCH --job-name=just_hmm_really	                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=32		                            # Number of cores per task
#SBATCH --mem=250gb			                                # Total memory for job
#SBATCH --time=1-00:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/Tn5_sequencing/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/Tn5_sequencing/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

PROJECT=Tn5_sequencing
BATCH_NAME=ancestry_v3

SCRIPTS_DIR=~/Git_repos/${PROJECT}/${BATCH_NAME}
LOG_DIR=/home/mcf96392/logs/${PROJECT}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

NTHREADS=32
MEM=250

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

GROUP=allnorth

cd ${WORKING_DIR}
mkdir -p ${WORKING_DIR}/ancestry_v3_${GROUP}
cd ${WORKING_DIR}/ancestry_v3_${GROUP}

#make readlist



source ~/ancestryinfer/sourceme
ml Miniconda3/4.7.10
source activate ~/ancestryinfer/ngsutil
#source activate py2
python --version
conda info
samtools --version

CHR=$1
printf "CHR is %s\n" $CHR | tee >(cat >&2)

if [ ! -f hmm_batch_${CHR}_first7.sh ]; then
	printf "Truncating hmm_batch script... %s\n" $CHR | tee >(cat >&2)
	cat hmm_batch.sh | sed "s/Chr14/${CHR}/g" | head -n7 > hmm_batch_${CHR}_first7.sh
fi

if [ ! -f ${WORKING_DIR}/allnorth_withcMs/perlscript_182.pl ]; then
	printf "Truncating perl script... %s\n" $CHR | tee >(cat >&2)
	head -n182 /home/mcf96392/ancestryinfer/ancestryinfer/combine_all_individuals_hmm_v5.pl > ${WORKING_DIR}/allnorth_withcMs/perlscript_182.pl
fi

if [ ! -f all.indivs.hmm.combined_${CHR}.filtered_focalchroms ] || [ ! -f current.samples.list_${CHR} ]; then
	bash hmm_batch_${CHR}_first7.sh

	string=$(perl ${WORKING_DIR}/allnorth_withcMs/perlscript_182.pl HMM.parental.files.list_${CHR} HMM.hybrid.files.list_${CHR} 0.5 /scratch/mcf96392/gutnas_panel/samtools_gt/AIMs_panel30_208560_counts.txt 0 /scratch/mcf96392/ancestry_v3/ancestry_v3_allnorth/focal_chrom_list_${CHR} 150 0.02 _${CHR} | grep hmm.combined.pass.formatted | grep -v formatted.posterior)
	paste $string > all.indivs.hmm.combined_${CHR}.filtered_focalchroms

	printf "Created files, checking file heads... %s\n" $CHR | tee >(cat >&2)
	head -c 20 all.indivs.hmm.combined_${CHR} | tee >(cat >&2)
	head -n2 current.samples.list_${CHR} | tee >(cat >&2)
fi


printf "Adding cMs to input matrix... %s\n" $CHR | tee >(cat >&2)
cd ${WORKING_DIR}/allnorth_withcMs
cp ${WORKING_DIR}/ancestry_v3_allnorth/current.samples.list_${CHR} current.samples.list_${CHR}
cut -f1-6 ${WORKING_DIR}/ancestry_v3_allnorth/all.indivs.hmm.combined_${CHR}.filtered_focalchroms > ${CHR}_first6

awk 'BEGIN{OFMT="%.9f";a=0}{b=$2}{print 0.000000039*(b-a)}{a=$2}' ${CHR}_first6 > ${CHR}_cMs.txt

paste ${CHR}_first6 ${CHR}_cMs.txt <(cut -f7- ../ancestry_v3_allnorth/all.indivs.hmm.combined_${CHR}.filtered_focalchroms) > all.indivs.hmm.combined_${CHR}.filtered_focalchroms_withcMs


printf "\nRunning hmm... %s\n" $CHR | tee >(cat >&2)
#perl /home/mcf96392/ancestryinfer/ancestryinfer/combine_all_individuals_hmm_v5.pl HMM.parental.files.list_${CHR} HMM.hybrid.files.list_${CHR} 0.5 /scratch/mcf96392/gutnas_panel/samtools_gt/AIMs_panel30_208560_counts.txt 0 /scratch/mcf96392/ancestry_v3/ancestry_v3_allnorth/focal_chrom_list_${CHR} 150 0.02 _${CHR}
ancestry_hmm -a 2 0.5 0.5 -p 0 -100 0.5 -p 1 -100 0.5 -e 0.02 -s current.samples.list_${CHR} -i all.indivs.hmm.combined_${CHR}.filtered_focalchroms_withcMs
printf "\nEnd of the script... %s\n" $CHR | tee >(cat >&2)

### cannot open combined HMM input data <-look for this // submitted through Chr04
perl /home/mcf96392/ancestryinfer/ancestryinfer/convert_rchmm_to_ancestry_tsv_v3.pl ${WORKING_DIR}/ancestry_v3_allnorth/current.samples.list ${WORKING_DIR}/ancestry_v3_allnorth/current.samples.read.list 1 ${WORKING_DIR}/ancestry_v3_allnorth/focal_chrom_list_${CHR}
perl /home/mcf96392/ancestryinfer/ancestryinfer/transpose_tsv.pl ancestry-probs-par1_transposed_${CHR}.tsv
perl /home/mcf96392/ancestryinfer/ancestryinfer/transpose_tsv.pl ancestry-probs-par2_transposed_${CHR}.tsv
perl /home/mcf96392/ancestryinfer/ancestryinfer/parsetsv_to_genotypes_v2.pl ancestry-probs-par1_${CHR}.tsv ancestry-probs-par2_${CHR}.tsv 0.9 ancestry-probs_${CHR}.tsv_rec.txt
Rscript /home/mcf96392/ancestryinfer/ancestryinfer/identify_intervals_ancestryinfer.R ancestry-probs_${CHR}.tsv_rec.txt /home/mcf96392/ancestryinfer/ancestryinfer
