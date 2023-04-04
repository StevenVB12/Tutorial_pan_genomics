#!/usr/bin/env python

from __future__ import division
from Bio import SeqIO
import argparse, re

parser = argparse.ArgumentParser()
parser.add_argument("-I", "--infile", help="fasta", required = True)
parser.add_argument("-o", "--outfile", help="txt", required = True)
parser.add_argument("-g", "--genomeList", help="comma seperated list of genomes in order of seq-seq-pan", required = True)
parser.add_argument("-b", "--bedfile", help="bedfile of intersect", required = True)
args = parser.parse_args()

infile  = args.infile
outfile = args.outfile
bedfile = args.bedfile
genomes = args.genomeList.split(",")

bedfile = open(bedfile, "r")

records = list(SeqIO.parse(infile, "fasta"))

seq1 = records[0].seq
seq2 = records[1].seq

#print(seq1)
#print(seq2)

print(len(seq1))
print(len(seq2))
    
bedLine = bedfile.readline()
bedLine = bedLine.strip()

f = open(outfile + '.txt',"w+")
f.write("\t".join(("start", "end", "length", "IDY", "missing1", "missing2", "\n")))
    
while bedLine:
    
    bedElem =  bedLine.split('\t')
    
    inter1 = seq1[int(bedElem[1]):int(bedElem[2])]
    inter2 = seq2[int(bedElem[1]):int(bedElem[2])]

    if(int(bedElem[1]) == int(bedElem[2])): continue
        
    score = 0
    inter1_N = 0
    inter2_N = 0
    
    
    for i in range(0,len(inter1)):

        inter1_base = inter1[i]
        inter2_base = inter2[i]

        if inter1_base == inter2_base:
            score = score + 1
        if inter1_base == '-':
            inter1_N = inter1_N + 1
        if inter2_base == '-':
            inter2_N = inter2_N + 1
   
    f.write("\t".join((bedElem[1],
                       bedElem[2],
                       str(len(inter1)),
                       str(round(score/len(inter1),4)),
                       str(round(inter1_N/len(inter1),4)),
                       str(round(inter2_N/len(inter1),4)),"\n")))        

    
    bedLine = bedfile.readline()
    bedLine = bedLine.strip()
    
    
