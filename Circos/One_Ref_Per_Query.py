#!/usr/bin/env python
import os
import sys
import inspect
import re
import argparse
import random

parser = argparse.ArgumentParser(description="""           
Description
-----------
This script keep only alignements for the best Ref-Query pair.
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Vincent Merel
""")

#Input file
parser.add_argument("--input", type=str, required=True, dest="input", default=None, help="the input file")
#Output file
parser.add_argument("--output", type=str, required=True, dest="output", default=None, help="the output file")

args = parser.parse_args()

input = open(args.input,"r")
output = open(args.output,"w")

##############################################################################################################################
########################################Creating a List with selected couples of Contig########################################
##############################################################################################################################
myDict = {}
List = []

for line in input:
  
  Q_length=int(line.split()[5])
  Q_Contig=line.split()[7]
  R_Contig=line.split()[6]

  if Q_Contig not in myDict :
    
    myDict[Q_Contig]={}
    myDict[Q_Contig][R_Contig]=Q_length
    
  else :
    
    if R_Contig not in myDict[Q_Contig] :
      
      myDict[Q_Contig][R_Contig]=Q_length
    
    else :
      
      myDict[Q_Contig][R_Contig]=myDict[Q_Contig][R_Contig]+Q_length

for key in myDict :

	List.append(max(myDict[key], key=myDict[key].get)+" "+key)
  
input.close()
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
##############################################Writting selected alignement####################################################
##############################################################################################################################
input = open(args.input,"r")

for line in input:
  
  Q_Contig=line.split()[7]
  R_Contig=line.split()[6]
  Couple=R_Contig+" "+Q_Contig
  
  if Couple in List:
    output.write(line)
    
output.close()
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

'''
python /home/lerat/REDSKIN/one_ref_per_query.py  \
--input /home/lerat/Software/circos-0.69-9/data/Links.txt  \
--output /home/lerat/Software/circos-0.69-9/data/one_ref_per_query.txt
'''
