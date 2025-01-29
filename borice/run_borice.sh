#!/bin/bash
#SBATCH --job-name=borice                     # Job name
#SBATCH --partition=batch	                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=32		                            # Number of cores per task
#SBATCH --mem=120gb			                                # Total memory for job
#SBATCH --time=7-00:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/borice_hq/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/borice_hq/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

###SETUP

WORKING_DIR=$1
BORICE_CONTROL=$2 #~/Git_repos/Tn5_sequencing/borice_hq/newinputs/Master_Control.v3.txt
BORICE_SUBPOPS=$3 #~/Git_repos/Tn5_sequencing/borice_hq/newinputs/${GROUP}_${TAG}_borice_${POPTYPE}.txt
BORICE_GENOTYPES=$4 #${WORKING_DIR}/boriceready.${GROUP}_${TAG}.IM62_v3.${SNPLIST}.${CHR}.genotypes.txt
RUNPREFIX=$5 #output files will be in ${WORKING_DIR}/${RUNPREFIX}

BORICE_PATH=~/apps/BORICE.genomic.v3/borice.v3.mod.py
LOCUSPOLICE=~/apps/BORICE.genomic.v3/locus_police/locus.police.notab.py
REMOVE_IMPOSSIBLES=~/Git_repos/Tn5_sequencing/borice_hq/remove_impossible_SNPs.py

ml scipy/1.4.1-foss-2019b-Python-3.7.4

mkdir -p $WORKING_DIR
cd $WORKING_DIR

#locus police
PREFIX=${BORICE_GENOTYPES%.genotypes.txt}
python $LOCUSPOLICE $PREFIX
python $REMOVE_IMPOSSIBLES ${PREFIX}.genotypes.txt ${PREFIX}.impossibles.bySNP.txt > ${PREFIX}.cleaned.genotypes.txt

mkdir -p ${WORKING_DIR}/${RUNPREFIX}
cd ${WORKING_DIR}/${RUNPREFIX}

mv ${PREFIX}.cleaned.genotypes.txt ${WORKING_DIR}/${RUNPREFIX}
cp -u $BORICE_CONTROL ${WORKING_DIR}/${RUNPREFIX}/Control.v3.txt
cp -u $BORICE_SUBPOPS ${WORKING_DIR}/${RUNPREFIX}

python $BORICE_PATH $RUNPREFIX
