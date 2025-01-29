#!/usr/bin/python


### script to remove SNPs that have impossible mother-offspring configurations for Borice, allowing borice to run.

###USAGE: python ./remove_impossible_SNPs.py genotypes.txt impossibles.txt > cleaned.genotypes.txt

###to generation impossibles.txt, use the script locus.police.py provided with BORICE.genomic.v3 (https://github.com/jkkelly/BORICE.genomic.v3)
###impossibles.txt file format should be "SNP_Chr01_52894\tcount"
###by default, exclude all snps with count>0

import sys

genotypes_file=sys.argv[1]
impossibles_file=sys.argv[2]

if len(sys.argv) > 3:
    output=open(sys.argv[3],'w')
else:
    output=sys.stdout

###read impossibles file to list
removals=[]
with open(impossibles_file,'r') as impossibles:
    for line in impossibles:
        data=line.strip().split()
        if int(data[1]) > 0:
            snpid=data[0].split('_')
            removals.append(snpid[1]+'_'+snpid[2])

###read genotypes file and remove impossibles
with open(genotypes_file,'r') as genotypes:
    keep=True
    for line in genotypes:
        data=line.strip().split()
        if data[0] == "SNP":
            #new SNP, check for removal
            if data[1] in removals:
                keep=False
            else:
                keep=True
        #print line if active SNP is NOT in removals list
        if keep:
            output.write(line)
#close output
if len(sys.argv) > 3:
    output.close()
