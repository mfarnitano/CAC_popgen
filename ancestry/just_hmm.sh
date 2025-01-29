#!/bin/bash
#SBATCH --job-name=just_hmm	                        # Job name
#SBATCH --partition=highmem_p		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=32		                            # Number of cores per task
#SBATCH --mem=950gb			                                # Total memory for job
#SBATCH --time=168:00:00  		                            # Time limit hrs:min:sec
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
MEM=950

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

SL_ID=$(grep ${CHR} /home/mcf96392/logs/Tn5_sequencing/ancestry_hmm_jobids.txt | cut -f2)
printf "SL_ID is %s\n" $SL_ID | tee >(cat >&2)

if [ ! -f hmm_batch_${CHR}_first7.sh ]; then
	cat hmm_batch.sh | sed "s/Chr14/${CHR}/g" | head -n7 > hmm_batch_${CHR}_first7.sh
fi

bash hmm_batch_${CHR}_first7.sh

string=$(grep hmm.combined.pass.formatted slurm-${SL_ID}.out | grep -v formatted.posterior)
paste $string > all.indivs.hmm.combined_${CHR}

printf "Created file, checking file head... %s\n" $CHR | tee >(cat >&2)
head -c 20 all.indivs.hmm.combined_${CHR} | tee >(cat >&2)

printf "\nRunning perl script for hmm... %s\n" $CHR | tee >(cat >&2)
perl /home/mcf96392/ancestryinfer/ancestryinfer/combine_all_individuals_hmm_v5.pl HMM.parental.files.list_${CHR} HMM.hybrid.files.list_${CHR} 0.5 /scratch/mcf96392/gutnas_panel/samtools_gt/AIMs_panel30_208560_counts.txt 0 /scratch/mcf96392/ancestry_v3/ancestry_v3_allnorth/focal_chrom_list_${CHR} 150 0.02 _${CHR}
#ancestry_hmm -a 2 0.5 0.5 -p 0 -100 0.5 -p 1 -100 0.5 -e 0.02 -s current.samples.list_${CHR} -i all.indivs.hmm.combined_${CHR}.filtered_focalchroms
printf "\nEnd of the script... %s\n" $CHR | tee >(cat >&2)

### cannot open combined HMM input data <-look for this // submitted through Chr04
