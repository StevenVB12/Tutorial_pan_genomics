#!/usr/bin/env python

# Use as: python seq-seq-pan_blocks_intervals.py -I SeqSeqPan_erato_melp_noNewline.xmfa -g 1,2

from __future__ import division
from Bio.Seq import Seq
import argparse, re

parser = argparse.ArgumentParser()
parser.add_argument("-I", "--infile", help="xmfa file with no newlines from seq-seq-pan", required = True)
parser.add_argument("-g", "--genomeList", help="comma separated list of genomes in order of seq-seq-pan genome", required = True)
args = parser.parse_args()

infile = args.infile
genomes = args.genomeList.split(",")

count = 0

for genome in genomes:

    identified = 0
    sequence = str()
    sequenceOr = str()
    
    f= open(genome + '_blocks_intervals.txt',"w+")
    #f.write("\t".join(("start","end","\n")))
    count +=1
    fileIn = open(infile, "r")
    fileLine = fileIn.readline()
    fileLine = fileLine.strip()
    
    print('searching genome ' + str(count) + ' ' + str(genome))
    
    while fileLine:
    
        if '> ' + str(genome) + ':' in fileLine:
                     
            identified += 1
            
            fileLine = fileIn.readline()
            fileLine = fileLine.strip()
            seq = Seq(fileLine)
                            
            sequence = sequence + seq

        if '=' not in fileLine and '>' not in fileLine:
            lenSeq = len(fileLine)

        if '=' in fileLine and identified == 0:
            sequence = sequence + '-' * lenSeq

        if '=' in fileLine and identified == 1:
            identified = 0
        
        fileLine = fileIn.readline()
        fileLine = fileLine.strip()
    
    print('total length: ' + str(len(sequence)))        
    start = 0
    for i in range(len(sequence)):
        
        if i > 0 and i < len(sequence)-1:
            if sequence[i-1] == '-':
                if sequence[i] != '-':
                    start = i
            if sequence[i+1] == '-':
                if sequence[i] != '-':
                    end = i+1
                    f.write("\t".join(("pan",str(start),str(end),"\n")))

        if sequence[i-1] != '-' and i == len(sequence)-1:
            f.write("\t".join(("pan",str(start),str(i),"\n")))
    
    sequence = str()
    sequenceOr = str()
