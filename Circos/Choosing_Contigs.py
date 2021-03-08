#!/usr/bin/env python
import os
import sys
import inspect
import re
import argparse
import random
from Bio import SeqIO
from heapq import nlargest 

parser = argparse.ArgumentParser(description="""           
                                 Description
                                 -----------
                                 This script write selected contigs to a text file
                                 """,formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog="""
                                 Authors
                                 -------
                                 Vincent Merel
                                 """)

#Input files
parser.add_argument("--Q_Species", type=str, required=True, dest="Q_Species", default=None, help="the query species fasta")
parser.add_argument("--S_Species", type=str, required=True, dest="S_Species", default=None, help="the subject species fasta")
#Output files
parser.add_argument("--output", type=str, required=True, dest="output", default=None, help="the output file")
#Additional arugement
parser.add_argument("--Q_Contigs", type=str, required=True, dest="Q_Contigs", default=None, help="either an integer corresponding to the X longest contigs to use or a space separated list of contigs to use")
parser.add_argument("--S_Contigs", type=str, required=True, dest="S_Contigs", default=None, help="either an integer corresponding to the X longest contigs to use or a space separated list of contigs to use")

args = parser.parse_args()

#Opening Output files
output = open(args.output,"w")

######################Let's start with the Query######################

QDic={} #A dictionnary with Contigs as Keys and length as value

for record in SeqIO.parse(args.Q_Species, "fasta"):
  
  QDic[record.id]=len(record)

if str.isdigit(args.Q_Contigs)==True :
  Highest = nlargest(int(args.Q_Contigs), QDic, key = QDic.get)
  for key in Highest :
    output.write(key+" 0 "+str(QDic[key])+"\n")
else :
  QContigs=args.Q_Contigs.split(" ")
  for key in QDic :
    if key in QContigs :
      output.write(key+" 0 "+str(QDic[key])+"\n")

######################################################################
########################And now the Subject###########################

SDic={} #A dictionnary with Contigs as Keys and length as value

for record in SeqIO.parse(args.S_Species, "fasta"):
  
  SDic[record.id]=len(record)

if str.isdigit(args.S_Contigs)==True :
  Highest = nlargest(int(args.S_Contigs), SDic, key = SDic.get)
  for key in Highest :
    output.write(key+" 0 "+str(SDic[key])+"\n")
else :
  SContigs=args.S_Contigs.split(" ")
  for key in SDic :
    if key in SContigs :
      output.write(key+" 0 "+str(SDic[key])+"\n")

output.close()

'''python ~/REDSKIN/Choosing_Contigs.py \
--Q_Species /media/lerat/SWING/Review/Dsuz/Drosophila_suzukii.fasta.masked \
--S_Species /media/lerat/SWING/Review/Dmel/dmel-all-chromosome-r6.28.fasta.masked \
--Q_Contigs "15" \
--S_Contigs "2L 2R 3L 3R X 4" \
--output /media/lerat/SWING/REDSKIN/Chosen_Contigs.txt
'''
