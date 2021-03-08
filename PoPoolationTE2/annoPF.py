#!/usr/bin/env python
import sys
import os
import argparse
import subprocess
import shutil
import Bio
import subprocess
import glob

parser = argparse.ArgumentParser(description="""           
Description
-----------
    This script annotate pseudoF""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Author
-------
    Vincent Merel    
""")


parser.add_argument("--ToAnnotate", type=str, required=True, dest="ToAnnotate", default=None, help="Fasta to annotate")
parser.add_argument("--CM", type=str, required=True, dest="CM", default=None, help="Cross Mapping output")
parser.add_argument("--fasta", type=str, required=True, dest="fasta", default=None, help="fasta")

parser.add_argument("--Output", type=str, required=True, dest="Output", default=None, help="Output folder")

args = parser.parse_args()

################################################################
def checkEqual1(iterator):
    iterator = iter(iterator)
    try:
        first = next(iterator)
    except StopIteration:
        return True
    return all(first == rest for rest in iterator)
################################################################    
################################################################
#########A dictionnary of information about each Cluster########
#########ClusterInfo[Cluster_X]=[SF, Order, [Sequences]]########
##############A sequence to cluster dictionnary#################
################################################################
ClusterInfo={} 
Seq_to_Clust={}

ToAnnotate=open(args.ToAnnotate,'r')

for line in ToAnnotate:
  
  line=line[:-1]
  
  Cluster=line.split("\t")[4]
  SF=line.split("\t")[5]
  Order=line.split("\t")[6]
  Seq=line.split("\t")[0]
  Seq=Seq.split("#")[0]
  
  Seq_to_Clust[Seq]=Cluster
  
  #If the cluster is not already in ClusterInfo
  #Create an entry with SF, Order and the first sequence
  if Cluster not in ClusterInfo:
    ClusterInfo[Cluster]=[]
    ClusterInfo[Cluster].append(SF)
    ClusterInfo[Cluster].append(Order)
    ClusterInfo[Cluster].append([Seq])
  
  else:
    ClusterInfo[Cluster][2].append(Seq)

ToAnnotate.close()
'''
for Cluster in ClusterInfo:
  print(Cluster)
  print(ClusterInfo[Cluster])
'''
print("################################")
print("Ready to analyse "+str(len(ClusterInfo))+" Clusters")
print("################################")
################################################################
################################################################
################################################################


################################################################
######################For the id column#########################
###############ListOfSequencesNames[Short]=Full#################
################################################################
fasta=open(args.fasta,'r')

from Bio import SeqIO

ListOfSequencesNames={}

for seq_record in SeqIO.parse(fasta, "fasta"):
    Full=seq_record.id
    Short=Full.split("#")[0]
    ListOfSequencesNames[Short]=Full

#print(ListOfSequencesNames)

################################################################
######A dic with number of mapped Reads for each sequence#######
###################SeqReads[Sequence]=NReads####################
################################################################
SeqReads={}

CM=open(args.CM,'r')

for line in CM:

  Seq=line.split("\t")[4]
  Seq=Seq.split("#")[0]
  NReads=int(line.split("\t")[2])
  
  if Seq not in SeqReads:
    SeqReads[Seq]=NReads
  else:
    SeqReads[Seq]=SeqReads[Seq]+NReads
    
'''
for Sequence in SeqReads:
  print(Sequence)
  print(SeqReads[Sequence])
'''
################################################################
################################################################
################################################################

################################################################
###########Creating a list with CrossMapping Clusters###########
###################SeqReads[Sequence]=NReads####################
################################################################
CMClust=[]

CM=open(args.CM,'r')

nSeqWithoutMapping=0

for line in CM:

  Seq=line.split("\t")[4]
  Seq=Seq.split("#")[0]
  nReads=int(line.split("\t")[2])
  NReads=SeqReads[Seq]
  
  if NReads==0:
    nSeqWithoutMapping=nSeqWithoutMapping+1
  
  else:
    if nReads/NReads>0.01: #Here is a mapping
    
      Cluster=Seq_to_Clust[Seq]
      CMSeq=line.split("#")[0]
      CMCluster=Seq_to_Clust[CMSeq]
      
      if Cluster!=CMCluster: #Here is a new pair
        
        Pair=[Cluster,CMCluster]
        #print("#")
        #print("Pair: "+(" ").join(Pair))
        
        NewCluster=True #
        NewCMClust=True #Tant qu'on a pas la preuve du contraire
        
        cpt=0
        ClusterIndex=""
        CMClusterIndex=""
        
        for item in CMClust:
          
          if Cluster in item:
            ClusterIndex=cpt
          if CMCluster in item:
            CMClusterIndex=cpt
        
          cpt=cpt+1
          
        if ClusterIndex=="" and CMClusterIndex=="":  #New group
          CMClust.append([Cluster,CMCluster])
          #print("New group")
          
        if ClusterIndex=="" and CMClusterIndex!="":  #Add Cluster to a group
          CMClust[CMClusterIndex].append(Cluster)
          #print("Add Cluster to a group")
  
        if ClusterIndex!="" and CMClusterIndex=="":  #Add CMCluster to a group
          CMClust[ClusterIndex].append(CMCluster)
          #print("Add CMCluster to a group")
        
        if ClusterIndex!="" and CMClusterIndex!="":  #Merging two groups
          
          if ClusterIndex==CMClusterIndex:
            
            DoNothing=True
            #print("Already seen pair")
          
          else:
  
            #print("Merging two groups")
          
            for itemou in CMClust[CMClusterIndex]:
              CMClust[ClusterIndex].append(itemou)
          
            del CMClust[CMClusterIndex]


#for item in CMClust:
#  print(item)

################################################################
###########################Wrtting Output#######################
################################################################

ToAnnotate=open(args.ToAnnotate,'r')
output=open(args.Output,'w')
output.write('Sequence'+'\t'+'BH/Seq'+'\t'+'SF'+'\t'+'O'+'\t'+'f'+'\t'+'fSF'+'\t'+'fO'+'\t'+'family'+'\t'+'superfamily'+'\t'+'order'+'\t'+'id'+'\n')

NF=0 #to count F
NPF=0 #to count PF

for line in ToAnnotate:
  
  line=line[:-1]
  
  Short=line.split("\t")[0]
  
  Cluster=line.split("\t")[4]
  SF=line.split("\t")[5]
  isPF=False
  
  for PF in CMClust:
    if Cluster in PF:
        isPF=True
        break

  if isPF==False:
    
    NF=NF+1
    output.write(line+"\t"+Cluster+"\t"+ClusterInfo[Cluster][0]+"\t"+ClusterInfo[Cluster][1]+"\t"+ListOfSequencesNames[Short]+"\n")
  
  else:
    PFSFs=[]
    PFOs=[]
    for Cluster in PF:
      PFSFs.append(ClusterInfo[Cluster][0])
      PFOs.append(ClusterInfo[Cluster][1])
    #PF SF
    if checkEqual1(PFSFs)==True:
      PFSF=PFSFs[0]
    else:
      PFSF="Unknown"
      
    #PF O
    if checkEqual1(PFOs)==True:
      PFO=PFOs[0]
    else:
      PFO="Unknown"
    
    NPF=NPF+1
    output.write(line+"\t"+("_").join(PF)+"\t"+PFSF+"\t"+PFO+"\t"+ListOfSequencesNames[Short]+"\n")

print("################################")
print("NPseudoFamilies: "+str(len(CMClust)))
print("################################")

output.close()

'''
python ~/redskin/TE_db/CrossMapping/annoPF.py \
--ToAnnotate /home/lerat/REDSKIN/TE_db/all_TEs.txt \
--CM /home/lerat/REDSKIN/TE_db/CrossMapping/CrossMapping.txt \
--Output ~/REDSKIN/TE_db/CrossMapping/all_TEs.txt
  '''
           
