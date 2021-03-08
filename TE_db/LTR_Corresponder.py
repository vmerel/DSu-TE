#!/usr/bin/env python
import sys
import os
import argparse
import Bio

parser = argparse.ArgumentParser(description="""           
Description
-----------
    This script establish a correspondancy table LTR - Internal part""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Author
-------
    Vincent Merel    
""")



parser.add_argument("--input",
 type=str,
 required=True,
 dest="input",
 default=None,
 help="Fasta to merge")

parser.add_argument("--output",
 type=str,
 required=True,
 dest="output",
 default=None,
 help="Correspondancy table")

args = parser.parse_args()

from Bio import SeqIO

OK = {} #Dictionnary with solved families

cpt_Ambiguous=0
cpt_Paired=0 
cpt_NonPaired=0

with open(args.output,"w") as output:

	for seq_record in SeqIO.parse(args.input, "fasta"):

		BeforeHash = seq_record.id.split("#")[0]
		
		###e.g Gypsy-2_DBp-I#LTR/Gypsy and Gypsy-2_DBp-LTR#LTR/Gypsy###
		if BeforeHash.endswith('R'): #LTR
			family=BeforeHash[:-4]
			OK.setdefault(family,[]).append(seq_record.id)
		elif BeforeHash.endswith('I'): #Internal part
			family=BeforeHash[:-2]
			OK.setdefault(family,[]).append(seq_record.id)
		###############################################################

		###e.g BEL1-I_DV#LTR/Pao RepbaseID: BEL1-I_DVXX and BEL1-LTR_DV#LTR/Pao RepbaseID: BEL1-LTR_DVXX###
		elif "-I_" in BeforeHash :
			family=BeforeHash.replace("-I","")
			OK.setdefault(family,[]).append(seq_record.id)
		elif "-LTR_" in BeforeHash :
			family=BeforeHash.replace("-LTR","")
			OK.setdefault(family,[]).append(seq_record.id)
		###############################################################

		else :
			output.write(seq_record.id+"\n")			
			print("You will have to solve manually : " + seq_record.id)
			cpt_Ambiguous=cpt_Ambiguous+1

	for key in OK:
		output.write(key+"\t"+'\t'.join(OK[key])+"\n")

		if len(OK[key])==2:
			cpt_Paired=cpt_Paired+1
		elif len(OK[key])==1:
			cpt_NonPaired=cpt_NonPaired+1

print("Carefull ! Each line must have two or three columns for LTR_Merger.py to work correctly") 
print(str(cpt_Paired) + " sequences paired") 
print(str(cpt_NonPaired) + " sequences unpaired") 
print(str(cpt_Ambiguous) + " sequences ambiguous") 

output.close()


