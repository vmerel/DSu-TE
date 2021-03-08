#!/usr/bin/env python
import sys
import os
import argparse
import Bio

parser = argparse.ArgumentParser(description="""           
Description
-----------
    This script merged LTR and Internal Part using a correspondancy table""",formatter_class=argparse.RawDescriptionHelpFormatter,
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

parser.add_argument("--cor",
 type=str,
 required=True,
 dest="cor",
 default=None,
 help="Correspondancy table")

parser.add_argument("--output",
 type=str,
 required=True,
 dest="output",
 default=None,
 help="Output")

args = parser.parse_args()

from Bio import SeqIO

TE_dict = SeqIO.index(args.input, "fasta")	

with open(args.output,"w") as output:

	for line in open(args.cor, "r"):

		line=line[:-1]
		splittedline=line.split("\t")
		family=splittedline[0]


		if len(splittedline)==2 : #Unpaired

			SeqIO.write(TE_dict[splittedline[1]], output, "fasta")

		elif len(splittedline)==3 :
			if "-LTR" in splittedline[1] or "_LTR_" in splittedline[1] :
				#print("Family"+"\t"+"LTR"+"\t"+"I")
				LTR=TE_dict[splittedline[1]]
				I=TE_dict[splittedline[2]]
				TE=LTR+I+LTR
				TE.id=TE_dict[splittedline[1]].id
				TE.description=""
				SeqIO.write(TE, output, "fasta")
			else :
				#print("Family"+"\t"+"I"+"\t"+"LTR")
				LTR=TE_dict[splittedline[2]]
				I=TE_dict[splittedline[1]]
				TE=LTR+I+LTR
				TE.id=TE_dict[splittedline[1]].id
				TE.description=TE_dict[splittedline[1]].id
				SeqIO.write(TE, output, "fasta")
		else :
			print("You have a problem")
			print(splittedline)


		
output.close()
