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
This script filters alns in a coords file according to a user provided list of Contigs.
Are kept alns where both query and subject are included in this list.
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Vincent Merel
""")

#Input files
parser.add_argument("--coords", type=str, required=True, dest="coords", default=None, help="the coords File")
parser.add_argument("--Chosen_Contigs", type=str, required=True, dest="Chosen_Contigs", default=None, help="the Chosen_Contigs file")
#Output files
parser.add_argument("--output", type=str, required=True, dest="output", default=None, help="the output file")

args = parser.parse_args()

#########################Creating a list of choosen contigs#########################
Chosen_Contig_list=[]

Chosen_Contigs = open(args.Chosen_Contigs,"r")

for line in Chosen_Contigs :

	Chosen_Contig_list.append(line.split(" ")[0])

Chosen_Contigs.close()
#print(Chosen_Contig_list)
#####################################################################################

###########################Parsing input coords file#################################
#############Writting desired alns to output input coords file#######################
coords = open(args.coords,"r")
output = open(args.output,"w")

for line in coords :

	S_Contig=line.split()[6]
	Q_Contig=line.split()[7]

	if (S_Contig in Chosen_Contig_list and Q_Contig in Chosen_Contig_list) :

		output.write(line)

coords.close()
output.close()
#####################################################################################

