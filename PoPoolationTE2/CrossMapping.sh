IDs=`grep "^>" /media/vincent/0E9A457E451F81F1/REDSKIN_2/PoPoolationTE2/all_TEs.500.fa | tr -d ">"`

for ID in $IDs
do

#Extracting Sequence
  echo $ID
  awk -v P="$ID" '!c++{printf RS}$0~P' RS='>' /media/vincent/0E9A457E451F81F1/REDSKIN_2/PoPoolationTE2/all_TEs.500.fa > /media/vincent/0E9A457E451F81F1/REDSKIN_2/TE_db/CrossMapping/Sequence.fa

#Simulating Reads
  for cpt in 1 2 3 4 5 
  do
    python2 ~/Software/simulaTE/create-reads-for-te-sequences.py \
  --te-sequences /media/vincent/0E9A457E451F81F1/REDSKIN_2/TE_db/CrossMapping/Sequence.fa \
  --boost 20 \
  --read-length 125 \
  --output /media/vincent/0E9A457E451F81F1/REDSKIN_2/TE_db/CrossMapping/Reads_$cpt.fastq &
  done
  wait
  cat /media/vincent/0E9A457E451F81F1/REDSKIN_2/TE_db/CrossMapping/Reads_*fastq > /media/vincent/0E9A457E451F81F1/REDSKIN_2/TE_db/CrossMapping/Reads.fastq &
  wait
  rm /media/vincent/0E9A457E451F81F1/REDSKIN_2/TE_db/CrossMapping/Reads_*fastq

#Mapping
  bwa bwasw \
  -t 5 \
  /media/vincent/0E9A457E451F81F1/REDSKIN_2/PoPoolationTE2/all_TEs.500.fa \
  /media/vincent/0E9A457E451F81F1/REDSKIN_2/TE_db/CrossMapping/Reads.fastq |\
  samtools view -Sb --threads 5 |\
  samtools sort --threads 5 > /media/vincent/0E9A457E451F81F1/REDSKIN_2/TE_db/CrossMapping/Aligned.bam
  
  samtools index /media/vincent/0E9A457E451F81F1/REDSKIN_2/TE_db/CrossMapping/Aligned.bam
  samtools idxstats /media/vincent/0E9A457E451F81F1/REDSKIN_2/TE_db/CrossMapping/Aligned.bam > /media/vincent/0E9A457E451F81F1/REDSKIN_2/TE_db/CrossMapping/Tmp.txt

  awk -v id="$ID" 'BEGIN{FS=OFS="\t"}{$5=id}1' /media/vincent/0E9A457E451F81F1/REDSKIN_2/TE_db/CrossMapping/Tmp.txt >> /media/vincent/0E9A457E451F81F1/REDSKIN_2/TE_db/CrossMapping/CrossMapping.txt
done
