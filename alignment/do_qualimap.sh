#!/bin/bash
#SBATCH --job-name=qualimap                     # Job name
#SBATCH --partition=batch	                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=5		                            # Number of cores per task
#SBATCH --mem=20gb			                                # Total memory for job
#SBATCH --time=4:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/borice_hq/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/borice_hq/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

###SETUP
WORKING_DIR=$1
PREFIX=$2
REFCODE=$3

ml Qualimap/2.2.1-foss-2021b-R-4.1.2
printf "\nScript called with %s %s %s \n" $1 $2 $3 | tee >(cat >&2)
printf "\n...generating coverage summary stats with qualimap\n" | tee >(cat >&2)


qualimap bamqc -bam ${WORKING_DIR}/bams/${PREFIX}.${REFCODE}.fds.bam -c -outdir ${WORKING_DIR}/qualimap/${PREFIX}

printf "\nDone" | tee >(cat >&2)
