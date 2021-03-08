#!/bin/bash
date

for CPT in $(seq 1 1 16)
do
  ####################Reads simulation####################
  
  # Multithreaded simulation
  
  for cpt in $(seq 1 1 32)
  do
    python2 /home/ubuntu/create-reads-for-te-sequences.py \
  --te-sequences /home/ubuntu/all_TEs.fa \
  --boost 2 \
  --read-length 125 \
  --output /mnt/mydatalocal/Reads_"$cpt".fastq &
  done
  wait
  
  # Merging
  
  cat /mnt/mydatalocal/Reads_*fastq > /mnt/mydatalocal/Reads.fastq &
  wait
  
  # Removing unecessary files
  for cpt in $(seq 1 1 32)
  do
   rm /mnt/mydatalocal/Reads_"$cpt".fastq
  done
  
  ###############################################
  ####################Mapping####################
  cd /mnt/mydatalocal/
  
  # Index
  bwa index /mnt/mydatalocal/Drosophila_suzukii.fasta.masked
  
  # Mapping
  bwa bwasw /mnt/mydatalocal/Drosophila_suzukii.fasta.masked \
  /mnt/mydatalocal/Reads.fastq -t 32|\
  samtools view -Sb --threads 32 |\
  samtools sort --threads 32 > /mnt/mydatalocal/Aligned.bam
  
  ###############################################
  ####################Masking####################
  bedtools bamtobed -i /mnt/mydatalocal/Aligned.bam > /mnt/mydatalocal/bed.bed
  
  bedtools maskfasta -fi /mnt/mydatalocal/Drosophila_suzukii.fasta.masked \
    -fo /mnt/mydatalocal/Drosophila_suzukii.fasta.masked.tmp -bed /mnt/mydatalocal/bed.bed
  
    cat /mnt/mydatalocal/Drosophila_suzukii.fasta.masked.tmp|grep -v "^>"|grep "N"|perl -pe "s/[^N]//g"|wc -c >> /home/ubuntu/IterativeMasking.txt
    
    cp /mnt/mydatalocal/Drosophila_suzukii.fasta.masked.tmp /mnt/mydatalocal/Drosophila_suzukii.fasta.masked
  done
date
