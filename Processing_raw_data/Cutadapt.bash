#!/bin/bash
#SBATCH -J Cutadapt_PF
#SBATCH -A xxx
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10GB
#SBATCH --cpus-per-task=1

module load cutadapt/3.3-gimkl-2020a-Python-3.8.2

####NOTES####
#Please note, ALL primers have been given in this script and the will need to be adjusted for the different raw data files. This could be written into a loop.

#raw data directories and primers used:
#rawdata_16S  rawdata_18S_G1  rawdata_18S_G2
#Herlemann  FD-Bra;Zim  Stoeck;Zhan;Sherwood

#This can also be run in the terminal, as little computing power is needed.
###########

for i in ./rawdata/*_R1_001.fastq.gz;
do
  SAMPLE=$(echo ${i} | sed "s/_R1_\001\.fastq.gz//")
  echo ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_001.fastq.gz
  cutadapt \
        --discard-untrimmed \
        --action=trim \
        -e 0 \
        --no-indels \
        -g FD=^TTAAAAAGCKCGTAGTTG;min_overlap=18 -G ^ACTTTCGTTCTTGAT;min_overlap=15 \
        -g Zim=^ATTCCAGCTCCAATAGCG;min_overlap=18 -G ^GACTACGATGGTATCTAATC;min_overlap=20 \
        -g Stoeck=^CCAGCASCYGCGGTAATTCC;min_overlap=20 -G ^ACTTTCGTTCTTGATYRA;min_overlap=18 \
        -g Zhan=^AGGGCAAKYCTGGTGCCAGC; min_overlap=20 -G ^GRCGGTATCTRATCGYCTT;min_overlap=19 \
        -g Sherwood=^GGACAGAAAGACCCTATGAA;min_overlap=20 -G ^TCAGCCTGTTATCCCTAGAG;min_overlap=20 \
        -g Herlemann=^CCTACGGGNGGCWGCAG;min_overlap=17 -G ^GACTACHVGGGTATCTAATCC;min_overlap=21 \
        -o ./"${SAMPLE}_{name}_R1_001.fastq.gz" \
        -p ./"${SAMPLE}_{name}_R2_001.fastq.gz" \
        ${SAMPLE}_R1_001.fastq.gz  ${SAMPLE}_R2_001.fastq.gz

done
