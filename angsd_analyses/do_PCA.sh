#!/bin/bash
#SBATCH --job-name=angsd_structure	                  # Job name
#SBATCH --partition=batch	                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=10		                            # Number of cores per task
#SBATCH --mem=120gb			                                # Total memory for job
#SBATCH --time=24:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/borice_hq/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/borice_hq/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

WORKING_DIR=$1
INPUT_BEAGLE=$2
OUT_PREFIX=$3

printf "Script called with parameters $1 $2 $3 \n"

ml angsd/0.940-GCC-11.2.0
ml PCAngsd/1.10

cd $WORKING_DIR
pcangsd -b $INPUT_BEAGLE -t 8 -o $OUT_PREFIX --inbreedSamples
