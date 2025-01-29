#!/bin/bash
#!/bin/bash
#SBATCH --job-name=combine_ancestry	                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=1		                            # Number of cores per task
#SBATCH --mem=25gb			                                # Total memory for job
#SBATCH --time=1-00:00:00  		                            # Time limit hrs:min:sec
#SBATCH -o /home/mcf96392/logs/Tn5_sequencing/%x.%j.out		# Standard output log
#SBATCH -e /home/mcf96392/logs/Tn5_sequencing/%x.%j.err		#Standard error log
#SBATCH --mail-user=mcf96392@uga.edu                    # Where to send mail
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

sample_list=$1
outprefix=$2

for CHR in Chr{01..14}; do
	reference=$(head -n1 $sample_list)
	printf "" > temp.txt
	cat $sample_list | while read -r sample; do
		posterior=${sample}_read_1.fastq.gz_${CHR}.sam.hmm.combined.pass.formatted.posterior
		if [ ! -f $posterior ]; then
			sample=${sample}'-combined'
			posterior=${sample}_read_1.fastq.gz_${CHR}.sam.hmm.combined.pass.formatted.posterior
		fi
		printf "%s\n" $sample > ${sample}.${CHR}.calls
		tail -n +2 $posterior | grep $CHR | awk '{if ($3>=0.90) {print "G"} else if ($4>=0.90) {print "H"} else if ($5>=0.90) {print "N"} else {print "-"}}' >> ${sample}.${CHR}.calls

		filestring=$(cat temp.txt | sed 's/\n/ /g')
		filestring=${filestring}"space"${sample}.${CHR}".calls"
		printf '%s' $filestring > temp.txt
	done
	filestring=$(cat temp.txt | sed 's/\n/ /g' | sed 's/space/ /g')
	echo $filestring
	paste <(cut -f1,2 ${reference}_read_1.fastq.gz_${CHR}.sam.hmm.combined.pass.formatted.posterior | grep "chrom\|$CHR") $filestring > ${outprefix}_${CHR}.calls.txt
done
cat <(head -n1 ${outprefix}_Chr01.calls.txt) > ${outprefix}_allChrs.calls.txt
cat ${outprefix}_Chr{01..14}.calls.txt | grep -v chrom >> ${outprefix}_allChrs.calls.txt
