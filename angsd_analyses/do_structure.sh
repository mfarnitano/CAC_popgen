#!/bin/bash
#SBATCH --job-name=angsd_structure	                  # Job name
#SBATCH --partition=highmem_p	                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=10		                            # Number of cores per task
#SBATCH --mem=480gb			                                # Total memory for job
#SBATCH --time=168:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/borice_hq/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/borice_hq/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

WORKING_DIR=$1
INPUT_BEAGLE=$2
OUT_PREFIX=$3

printf "Script called with parameters $1 $2 $3 \n"

ml angsd/0.940-GCC-11.2.0

cd $WORKING_DIR
NGSadmix -likes $INPUT_BEAGLE -K 2 -P 8 -outfiles $OUT_PREFIX
