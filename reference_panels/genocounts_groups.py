#!/usr/bin/python

import sys

popfile=sys.argv[1]
genotypes=sys.argv[2]

groupings={}
grouplist=[]

with open(popfile,'r') as pops:
    for line in pops:
        entries=line.strip().split()
        groupings[entries[0]] = entries[1]
        if entries[1] not in grouplist:
            grouplist.append(entries[1])

sys.stderr.write("Groups to count: "+str(grouplist)+"\n")
order=[]

#print output header
sys.stdout.write("CHROM\tPOS\tREF\tALT")
for g in grouplist:
    sys.stdout.write("\t"+g+"_alleles\t"+g+"_na_counts\t"+g+"_nsamples")
sys.stdout.write("\n")

with open(genotypes,'r') as genos:
    header=True
    for line in genos:
        entries=line.strip().split()
        if header:
            order=[groupings[l] for l in entries[4:]]
            header=False
            continue
        chrom=entries[0]
        pos=entries[1]
        ref=entries[2]
        alt=entries[3]
        vals=entries[4:]
        if len(vals) != len(order):
            raise Exception("Error: number of entries does not match header")
        #initialize dict
        na_counts={}
        alleles={}
        nsamples={}
        for g in grouplist:
            na_counts[g]=0
            alleles[g]=0
            nsamples[g]=0
        for i in range(len(vals)):
            group=order[i]
            try:
                alleles[group]+=int(vals[i])
                nsamples[group]+=1
            except ValueError as ve:
                na_counts[group]+=1
                nsamples[group]+=1
        sys.stdout.write(chrom+"\t"+pos+"\t"+ref+"\t"+alt)
        for g in grouplist:
            sys.stdout.write("\t"+str(alleles[g])+"\t"+str(na_counts[g])+"\t"+str(nsamples[g]))
        sys.stdout.write("\n")
#DONE
