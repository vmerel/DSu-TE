#!/usr/bin/env python
import os
import sys
import inspect
import re
import argparse
import random
import itertools

parser = argparse.ArgumentParser(description="""           
Description
-----------
This script create a Circos Density Track from a GFF/GTF and the corresponding fasta
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Vincent Merel
""")

#Input files
parser.add_argument("--fasta", type=str, required=True, dest="fasta", default=None, help="the fasta file")
parser.add_argument("--gff", type=str, required=True, dest="gff", default=None, help="the gff file")
parser.add_argument("--karyotype", type=str, required=True, dest="karyotype", default=None, help="the karyotype file")
#Additional argument
parser.add_argument("--window", type=int, required=True, dest="window", default=None, help="the window size")
parser.add_argument("--feature", type=str, required=True, dest="feature", default=None, help="the feature")
parser.add_argument("--prefix", type=str, required=True, dest="prefix", default=None, help="the prefix")
#Output files
parser.add_argument("--output", type=str, required=True, dest="output", default=None, help="the output file")

#parser.add_argument("--links", type=str, required=True, dest="links", default=None, help="the links file")

args = parser.parse_args()

output=open(args.output,"w")




#########################Getting size of the Assembly on the plot#########################

karyotype = open(args.karyotype,"r")

Assembly_Size=0

for line in karyotype :
  
  line=line[:-1]
  
  if args.prefix in line:
    
    Assembly_Size=Assembly_Size+int(line.split()[5])

print("#########################################################")
print("Ready to process "+str(Assembly_Size)+"pb of Assembly")

karyotype.close()
##########################################################################################

#Get number of chunks
NChunks=round(Assembly_Size/args.window) #/!\ Round 
print(str(NChunks)+" Chunks of "+str(args.window)+"pb")
print("#########################################################")

################################Parsing assembly by chunk#################################


################
#####Input######
################
'''
chr - dq_Contig_87 87 0 639529 chr2
chr - dq_Contig_25 25 0 1895205 chr2
chr - dq_Contig_62 62 0 894058 chr2
chr - dq_Contig_143 143 0 340096 chr2
chr - dq_Contig_363 363 0 53226 chr2
chr - dq_Contig_51 51 0 1047734 chr2
chr - dq_Contig_109 109 0 471589 chr2
chr - dq_Contig_282 282 0 87555 chr2
chr - dq_Contig_174 174 0 235515 chr2
chr - dq_Contig_194 194 0 183327 chr2
'''
################
################
################

################
#Desired output#
#for 1mb Chunks#
################
'''
[Chunk] Contig      Start       End                           (Remaining)
1       Contig_87   0           639.529                       1.000.000-639.529=360.471
        Contig_25   0           360.471                        360.471-360.471=0 -> 1.000.00
2      Contig_25   360.471     1.360.47                       1.000.000-1.000.000=0 -> 1.000.00
3      Contig_25   1.360.472   1.895.205                      1.000.000-(1.895.205-1.360.472)=465267
       Contig_62   0           465267                         0 -> 1.000.000
4      Contig_62   465267      894058                         1.000.000-(894.058-465.267)=571209
       Contig_143  0           340096                         571.209-340.096=231.113
       Contig_363  0           53226                          231.113-53.226=177.887
       Contig_51   0           177887                         0 -> 1.000.000
5      Contig_51   177887      1047734                        1000000-(1047734-177887)=130153
       Contig_109  0           130153 
'''
################
################
################

with open(args.karyotype) as karyotype:
  Lines = karyotype.readlines()

from itertools import cycle
Lines=cycle(Lines)

Line=next(Lines) #Get first line and associated info
Contig=Line.split()[2]
Remaining_Contig_Length=int(Line.split()[5])
Start=0

Chunks={}
Chunk=0 #Initiation of the chunk cpt

while Chunk!=NChunks: #4

  Chunk=Chunk+1
  Chunks[Chunk]=[]
  Remaining=args.window #What you need to complete the chunk
  
  while Remaining > Remaining_Contig_Length: # Let's put full contigs if we can
  
    Chunks[Chunk].append(Contig+" "+str(Start)+" "+str(Start+Remaining_Contig_Length))
    Remaining=Remaining-Remaining_Contig_Length 
    
    Line=next(Lines)
    
    Contig=Line.split()[2] #7
    Start=0
    Remaining_Contig_Length=int(Line.split()[5])

  #Splitted Contig
  Chunks[Chunk].append(Contig+" "+str(Start)+" "+str(Start+Remaining)) #Not full contig

  Start=Start+Remaining
  Remaining_Contig_Length=Remaining_Contig_Length-Remaining

'''  
for Chunk in Chunks:
  print(str(item)+" "+"\n".join(Chunks[item]))
'''  

##########################################################################################

############################Calculating number of Genes by chunk##########################
import gffutils

gff = gffutils.create_db(args.gff, ":memory:", disable_infer_transcripts=True, disable_infer_genes=True)

for Chunk in Chunks:
  
  N=0
  
  for item in Chunks[Chunk]:
    
    
    
    #print(str(Chunk)+" "+item)
    Contig=item.split()[0][3:]
    Start=item.split()[1]
    End=item.split()[2]

    features =list(gff.region(seqid=Contig, start=Start, end=End, completely_within=True))
    
    for feature in features:
      
      if feature.featuretype==args.feature :
      
          N=N+1  
          
  for item in Chunks[Chunk]:
    
    output.write(item+" "+str(N)+"\n")      
  	        
'''

'''
