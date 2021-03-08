# TE-DSu

In this folder you will find scripts and TE database associated with the publication "The worldwide invasion of *Drosophila suzukii* is accompanied by a large increase of transposable element load and a small number of putatively adaptive insertions" (10.1101/2020.11.06.370932)

## TE database

In the TE_db folder:

- all_TEs.fa: all sequences in fasta format;
- all_TEs.txt: the corresponding annotation.
Colnames: Sequence, Best Hit (BH), BH Superfamily, BH Order, Family, Superfamily, Order.

## Scripts

In the SLiM folder:

- Run.SLiM: the SLiM script used to study watterson theta evolution.

In the TE_db folder:

- LTR_Corresponder.py: the script used to establish a correspondancy table LTR - Internal part
- LTR_Merger.py: the script used to merge LTR and Internal Part using a correspondancy table

In the PoPoolationTE2 folder:

- ItMap.sh: the script used to perfom an iterative masking of the genome assembly;
- CrossMapping.sh: the script used to generate info. about crossmapping TEs;
- annoPF.py: the script used to create a TE hierarchy (PoPoolationTE2 input file) using info. about CrossMapping TEs to create pseudofamilies.

In the circos folder:

- Circos.Rmd: the pipeline used to produce fig 1.C
- *.py: necessary python scripts


