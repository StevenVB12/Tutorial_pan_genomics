#!/usr/bin/env python

from __future__ import division
from Bio.Seq import Seq
import argparse, re

parser = argparse.ArgumentParser()
parser.add_argument("-I", "--infile", help="extract regions from seq-seq-pan", required = True)
parser.add_argument("-g", "--genomeList", help="comma seperated list of genomes in order of seq-seq-pan", required = True)
args = parser.parse_args()

infile = args.infile
genomes = args.genomeList.split(",")


count = 0

for genome in genomes:

    identified = 0
    sequence = str()
    
    f= open('genome_' + genome + '.fasta',"a+")
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
            #seq = Seq(fileLine)
                            
            sequence = sequence + fileLine

        if '=' not in fileLine and '>' not in fileLine:
            lenSeq = len(fileLine)

        if '=' in fileLine and identified == 0:
            sequence = sequence + '-' * lenSeq

        if '=' in fileLine and identified == 1:
            identified = 0
        
        fileLine = fileIn.readline()
        fileLine = fileLine.strip()
    
    print('total length: ' + str(len(sequence)))  

    f.write("".join((">genome_"+str(genome),"\n")))
    f.write("".join((sequence,"\n")))
      
    sequence = str()
                    

