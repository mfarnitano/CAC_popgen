
# Usage: likes2borice.py likelihoods_file borice_ids
# likelihoods file format:
#   header line for each SNP: SNP_ID\tREFallele\tALTallele\n
#   followed by three genotype probabilities hom_ref\thet\thom_alt\n, one line per sample in order
# borice_ids file format: header line: SNP\n then two columns, family_ID and sample_type(par or off)

import sys

with open(sys.argv[2],"r") as ids:
    samples=[x.strip()+'\t' for x in ids.readlines()]

counter=0
with open(sys.argv[1],"r") as genos:
    for g in genos:
        if counter>=len(samples):
            counter=0
        print(samples[counter]+g,end="")
        counter+=1
