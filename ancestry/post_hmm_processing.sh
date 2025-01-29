#!/bin/bash
#SBATCH --job-name=post_hmm_processing	                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=1		                            # Number of cores per task
#SBATCH --mem=25gb			                                # Total memory for job
#SBATCH --time=2-00:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/Tn5_sequencing/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/Tn5_sequencing/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

PROJECT=Tn5_sequencing
BATCH_NAME=ancestry_v3

SCRIPTS_DIR=~/Git_repos/${PROJECT}/${BATCH_NAME}
LOG_DIR=/home/mcf96392/logs/${PROJECT}
WORKING_DIR=/scratch/mcf96392/${BATCH_NAME}

NTHREADS=1
MEM=25

###LOG
cd $SCRIPTS_DIR
echo "${SLURM_JOB_ID},${SLURM_JOB_NAME},$(git rev-parse --short HEAD),$(date)" >> ${LOG_DIR}/joblog.txt

GROUP=allnorth

cd ${WORKING_DIR}
mkdir -p ${WORKING_DIR}/allnorth_withcMs
cd ${WORKING_DIR}/allnorth_withcMs



ml R/4.3.0-foss-2020b

for CHR in Chr{02..14}; do

	#bash ${SCRIPTS_DIR}/convert_rchmm_to_ancestry_tsv.sh ${WORKING_DIR}/ancestry_v3_allnorth/current.samples.list_${CHR} ${CHR}
	#perl /home/mcf96392/ancestryinfer/ancestryinfer/convert_rchmm_to_ancestry_tsv_v3.pl ${WORKING_DIR}/ancestry_v3_allnorth/current.samples.list ${WORKING_DIR}/ancestry_v3_allnorth/current.samples.read.list 1 ${WORKING_DIR}/ancestry_v3_allnorth/focal_chrom_list_${CHR}
	onesample=$(grep "done processing sample" ~/logs/Tn5_sequencing/post_hmm_processing.23606200.out | cut -d' ' -f4 | grep ${CHR} | head -n1)
	cut -f1,2 ${onesample}.posterior > ${CHR}_firsttwo.txt
	filestring1=$(grep "done processing sample" ~/logs/Tn5_sequencing/post_hmm_processing.23606200.out | cut -d' ' -f4 | grep ${CHR} | sed 's/\n/ /g' | sed 's/formatted/formatted.posterior.par1.results/')
	paste ${CHR}_firsttwo.txt $filestring1 > ancestry-probs-par1_transposed_${CHR}.tsv
	filestring2=$(grep "done processing sample" ~/logs/Tn5_sequencing/post_hmm_processing.23606200.out | cut -d' ' -f4 | grep ${CHR} | sed 's/\n/ /g' | sed 's/formatted/formatted.posterior.par2.results/')
	paste ${CHR}_firsttwo.txt $filestring2 > ancestry-probs-par2_transposed_${CHR}.tsv

	perl /home/mcf96392/ancestryinfer/ancestryinfer/transpose_tsv.pl ancestry-probs-par1_transposed_${CHR}.tsv
	perl /home/mcf96392/ancestryinfer/ancestryinfer/transpose_tsv.pl ancestry-probs-par2_transposed_${CHR}.tsv
	perl /home/mcf96392/ancestryinfer/ancestryinfer/parsetsv_to_genotypes_v2.pl ancestry-probs-par1_${CHR}.tsv ancestry-probs-par2_${CHR}.tsv 0.9 ancestry-probs_${CHR}.tsv_rec.txt
	Rscript /home/mcf96392/ancestryinfer/ancestryinfer/identify_intervals_ancestryinfer.R ancestry-probs_${CHR}.tsv_rec.txt /home/mcf96392/ancestryinfer/ancestryinfer
	printf "\nDone processing... %s\n" $CHR | tee >(cat >&2)
done
