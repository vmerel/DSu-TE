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

""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Vincent Merel
""")

#Input files
parser.add_argument("--input", type=str, required=True, dest="input", default=None, help="the input file")
parser.add_argument("--Contig_size", type=str, required=True, dest="Contig_size", default=None, help="the Contig size file")
#Output file
parser.add_argument("--output", type=str, required=True, dest="output", default=None, help="the output file")
parser.add_argument("--spacing", type=str, required=True, dest="spacing", default=None, help="the spacing")
parser.add_argument("--reversing", type=str, required=True, dest="reversing", default=None, help="the reversing")
parser.add_argument("--karyotype", type=str, required=True, dest="karyotype", default=None, help="the karyotype")
parser.add_argument("--links", type=str, required=True, dest="links", default=None, help="links")
parser.add_argument("--radius", type=str, required=True, dest="radius", default=None, help="radius")
#Additional arg
parser.add_argument("--colors", type=str, required=True, dest="colors", default=None, help="the colors")

args = parser.parse_args()

output = open(args.output,"w")

#Importing colors
colors=args.colors
colors=colors.split(";")


#Creating a Dic of contig Size
Sizes={}

Contig_size=open(args.Contig_size,'r')
for line in Contig_size:
  Sizes[line.split()[0]]=line.split()[1]
##############################################################################################################################
#########################################################To Sort##############################################################
##########################################Creating a Dictionnary of Dictionnary (Pos)#########################################
##################################Pos[R_Contig][Q_Contig]=[mean pos of aln 1, mean pos of aln 2, ...]#########################
##############################################################################################################################
######################################################To Reverse##############################################################
##############################################################################################################################
##########################################Creating a Dictionnary of Dictionnary (Sense)#######################################
##################################Sense[R_Contig][Q_Contig]=(R_End-Rstart)_aln1+(R_End-Rstart)_aln2+...#######################
##############################################################################################################################

#Creating an alignment Dictionnary		
input = open(args.input,"r")

Pos = {}
Sense = {} #To choose contig to reverse

for line in input :
	
  line=line[:-1]
	
  R_Start=int(line.split()[0])
  R_End=int(line.split()[1])
  Q_Start=int(line.split()[2])
  Q_End=int(line.split()[3])
  
  R_Contig=line.split()[6]
  Q_Contig=line.split()[7]

  if (Q_End-Q_Start)==0:
    print(line)
  if R_Contig not in Pos : #New Reference Contig and so Query Contig
	  
    Pos[R_Contig]={}
    Sense[R_Contig]={}
	  
    Pos[R_Contig][Q_Contig]=[(R_Start+R_End)/2]
    Sense[R_Contig][Q_Contig]=(Q_End-Q_Start)/abs(Q_End-Q_Start)*abs(R_End-R_Start) #if Q_Start > Q_End (R_End < R_Start doesn't exist ?), (Q_End-Q_Start)/(Q_End-Q_Start)=-1, else (Q_End-Q_Start)/(Q_End-Q_Start)+1
    
  elif Q_Contig not in Pos[R_Contig] : #Known Reference Contig but new Query Contig
  
    Pos[R_Contig][Q_Contig]=[(R_Start+R_End)/2]
    Sense[R_Contig][Q_Contig]=(Q_End-Q_Start)/abs(Q_End-Q_Start)*abs(R_End-R_Start)
    
  else :
    
    Pos[R_Contig][Q_Contig].append((R_Start+R_End)/2)
    Sense[R_Contig][Q_Contig]=Sense[R_Contig][Q_Contig]+(Q_End-Q_Start)/abs(Q_End-Q_Start)*abs(R_End-R_Start)
		  
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################


# Calculate average Ref pos for each Q contig
# Python program to get average of a list 
def Average(lst): 
    return sum(lst) / len(lst) 
  
for key in Pos :
  
  for i in Pos[key] :
    
    Pos[key][i]=Average(Pos[key][i])



##############################################################################################################################
################################################Ideogram section <spacing>####################################################
##############################################################################################################################
spacing=open(args.spacing,"w")

spacing.write(" <spacing>"+'\n')
spacing.write("   default = 0.01r"+'\n')

for key in sorted(Pos.keys()):
	
	plop=sorted(Pos[key], key=Pos[key].get, reverse=False)

	for i in range(0, len(plop)):
	  
	  if (i > 0): 
	
	    spacing.write("  <pairwise dq_"+plop[i-1]+" dq_"+plop[i]+">"+'\n')
	    spacing.write("   spacing = 0r"+'\n')
	    spacing.write("  </pairwise>"+'\n')

spacing.write(" </spacing>"+'\n')
spacing.close()
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
################################################Radius####################################################
##############################################################################################################################
radius=open(args.radius,"w")


beginning="chromosomes_radius="
Text=[]

for key in Pos:
  
  for key2 in Pos[key]:
    
    Text.append("dq_"+key2+":1.05r")

radius.write(beginning+",".join(Text))
	
radius.close()
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
#########################################################To reverse###########################################################
##############################################################################################################################
reversing=open(args.reversing,"w")

ToReverse=""

beginning="chromosomes_reverse="
Text=[]

#All reference chromosomes are reversed
for key in sorted(Pos.keys(), reverse=True):
  
  Text.append("dr_"+key)

for key in sorted(Pos.keys()):
	
	plop=sorted(Pos[key], key=Pos[key].get, reverse=False)

	for i in range(0, len(plop)):
	  
	  if Sense[key][plop[i]]<0:
	    
	    Text.append("dq_"+plop[i])
	    
reversing.write(beginning+",".join(Text))

reversing.close()
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################


##############################################################################################################################
#########################################################Karyotype############################################################
##############################################################################################################################	    
karyotype=open(args.karyotype,"w")

color_cpt=-1

Colors={}

for key in sorted(Pos.keys()):
	
	color_cpt=color_cpt+1

	plop=sorted(Pos[key], key=Pos[key].get, reverse=False)
	
	for i in range(0, len(plop)):
	  
	  splitted=plop[i].split("_")
	  while("" in splitted) : 
	    splitted.remove("")
	  
	  Colors[plop[i]]=colors[color_cpt]
	  karyotype.write("chr - dq_"+plop[i]+" "+splitted[len(splitted)-1]+" 0 "+Sizes[plop[i]]+" "+colors[color_cpt]+"\n")

chr=["4","X","3R","3L","2R","2L"]

print(colors)

for key in sorted(Pos.keys(), reverse=True):
  
  print(key)
  print(Sizes[key])
  Colors[key]=colors[color_cpt]
  karyotype.write("chr - dr_"+key+" "+key+" 0 "+Sizes[key]+" "+"255,255,255"+"\n") #185,199,205
  color_cpt=color_cpt-1
  
'''	    
for key in sorted(Pos.keys(), reverse=True):
  
  Colors[key]=color_cpt
  karyotype.write("chr - dr_"+key+" "+key+" 0 "+Sizes[key]+" chr"+str(color_cpt)+"\n")
  color_cpt=color_cpt-1
'''
karyotype.close()
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################


##############################################################################################################################
#########################################################Links############################################################
##############################################################################################################################	    
links=open(args.links,"w")

input = open(args.input,"r")

cpt=0

for line in input :
  
  cpt=cpt+1
  
  Q_Contig=line.split()[7]
  R_Contig=line.split()[6]
  Q_Start=line.split()[2]
  Q_End=line.split()[3]
  R_Start=line.split()[0]
  R_End=line.split()[1]
  
  links.write(str(cpt)+" dr_"+R_Contig+" "+str(R_Start)+" "+str(R_End)+" color="+str(Colors[R_Contig])+"\n")
  links.write(str(cpt)+" dq_"+Q_Contig+" "+str(Q_Start)+" "+str(Q_End)+" color="+str(Colors[Q_Contig])+"\n")

links.close()
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

