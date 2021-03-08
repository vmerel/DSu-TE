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
This script color links according to a track threshold
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Vincent Merel
""")

#Input files
parser.add_argument("--input", type=str, required=True, dest="input", default=None, help="the input file")
parser.add_argument("--genes", type=str, required=True, dest="genes", default=None, help="the genes file")

#Output file
parser.add_argument("--output", type=str, required=True, dest="output", default=None, help="the output file")

#Additional
parser.add_argument("--colors", type=str, required=True, dest="colors", default=None, help="the colors")
parser.add_argument("--alternate_colors", type=str, required=True, dest="alternate_colors", default=None, help="the alternate colors")

args = parser.parse_args()

output = open(args.output,"w")


#Colors
colors=args.colors
colors=colors.split(";")

alternate_colors=args.alternate_colors
alternate_colors=alternate_colors.split(";")

ColorDic={}
cpt=-1

for color in colors:
  cpt=cpt+1
  ColorDic[color]=alternate_colors[cpt]

print(ColorDic)
###################Dic with gene density####################
genes = open(args.genes,"r")

Dic={}

for line in genes:
  
  Contig=line.split()[0]
  Start=line.split()[1]
  End=line.split()[2]
  Density=line.split()[3]
  
  if Contig in Dic:
    
    Dic[Contig].append([Start,End,Density])
  
  else:
    
    Dic[Contig]=[]
    Dic[Contig].append([Start,End,Density])
###########################################################

links = open(args.input,"r")
lines = links.readlines()

for i in range(0, len(lines)):
	
	if (i%2 == 0): #Pair line
	
	  R_line=lines[i][:-1]
	  Q_line=lines[i+1][:-1]
	  
	  #print(Q_line)
	  
	  R_splitted=R_line.split()
	  Q_splitted=Q_line.split()
	  
	  Q_Contig=Q_splitted[1]
	  Q_Start=int(Q_splitted[2])
	  Q_End=int(Q_splitted[3])
	  #Looking for the corresponding density
	  
	  if Q_Contig in Dic:
	    
	    #print(Chunk)
	    Union_A=0
	    Union_B=0
	    
	    Found=False
	    Fully_Inside=False
	    TheGoodChunk=[]
	     
	    Density=0
	    
	    for Chunk in Dic[Q_Contig]:
	      
	      Chunk_Start=int(Chunk[0])
	      Chunk_End=int(Chunk[1])

	      
	      if (Chunk_Start<=Q_Start<=Chunk_End and Chunk_Start<=Q_End<=Chunk_End):
	        
	        TheGoodChunk=Chunk
	        Found=True
	        Fully_Inside=True
	        break
	        
	      elif Chunk_Start<=Q_Start<=Chunk_End:
	        
	        Union_A=Chunk_End-Q_Start
	        Chunk_A=Chunk
	        Found=True
	     
	      elif Chunk_Start<=Q_End<=Chunk_End:
	        
	        Union_B=Q_End-Chunk_Start
	        Chunk_B=Chunk
	        Found=True
          
	    if Found==True and Fully_Inside==False:
	      
	      if Union_A>Union_B:
	        TheGoodChunk=Chunk_A
	      else:
	        TheGoodChunk=Chunk_B
      
	    if len(TheGoodChunk)==3:
	      
	      Density=int(Chunk[2])
	      
	      if Density<6.5:
	        
	        R_splitted[4]="color="+ColorDic[R_splitted[4][6:]]
	        Q_splitted[4]="color="+ColorDic[Q_splitted[4][6:]]
	  
	  output.write(" ".join(R_splitted)+"\n")
	  output.write(" ".join(Q_splitted)+"\n")

'''
links = open(args.input,"r")
lines = links.readlines()

for i in range(0, len(lines)):
	
	if (i%2 == 0): #Pair line

		R_line=lines[i][:-1]
		Q_line=lines[i+1][:-1]
		
		R_splitted=R_line.split()
		Q_splitted=Q_line.split()

		Q_Contig=Q_splitted[1]
		R_Contig=R_splitted[1]
		
		if Q_Contig in Dic :
		  
		  Q_Color=Dic[Q_Contig]
		  
		if R_Contig in Dic :
		  
		  R_Color=Dic[R_Contig]
		
		if Q_Color==R_Color :  
		
		  output.write(R_line+" color="+Q_Color+"\n")
		  output.write(Q_line+" color="+R_Color+"\n")
		
		#else :
		  
		  #output.write(R_line+" color="+"lgrey"+"\n")
		  #output.write(Q_line+" color="+"lgrey"+"\n")

output.close()
karyotype.close()

		if R_Contig in myDict :

			if Q_Contig in myDict[R_Contig] :

				myDict[R_Contig][Q_Contig]=myDict[R_Contig][Q_Contig]+R_Start


for key in myDict :
	
	plop=sorted(myDict[key], key=myDict[key].get)
	
	for item in plop :
	  
	  output.write(item+" chr"+key.split("_")[1]+"\n")
	#output.write(max(myDict[key], key=myDict[key].get)+" "+key+"\n")

for key in myDict :
	
	output.write(key+" chr"+key.split("_")[1]+"\n")
	


python /home/lerat/REDSKIN/coloring_links.py  \
--input /home/lerat/Software/circos-0.69-9/data/Links.txt  \
--karyotype /media/lerat/SWING/REDSKIN/Circos/Karyotype.txt \
--output /media/lerat/SWING/REDSKIN/Circos/Links.txt
'''

