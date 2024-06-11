#!/bin/bash
#SBATCH -J Q2_P1_PF
#SBATCH -A xxx
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --cpus-per-task=2

####NOTES####
#Create paths for each primers raw data sets where primer has been found and trimmed.
#This code was created post analysis for another project, but could easily be implemented here to make things easier.
#Create a loop to run primer datasets independently.
#primers.txt is a list of primers targeted with each name on a newline.

cd raw_data

for i in $(cat ../primers.txt); \
do \
find $(pwd) -maxdepth 1 -type f -not -path '*/\.*' | grep ${i} | grep "_R1_" | awk '{print $0"\t"$1}' | awk '{print $0"\t"$2}' | awk '{sub(/.*\//,"",$1)}1' | awk '{sub(/_.*/,"",$1)}1' | sed 's/_R1_/_R2_/2' | sed '1s/^/sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n/' > path_${i}.txt; \
done

mv path*.txt ../
cd ../

sed -i 's/ /\t/g' path*.txt
############

module load QIIME2/2021.2

for i in $(cat primers.txt); \
do \

#Importing our files into qiime format
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ./path_${i}.txt \
  --output-path paired_end_demux_${i}.qza \
  --input-format PairedEndFastqManifestPhred33V2


#Visualization of sequence quality - sequence quality control - Viewing a summary of joined data with read quality
qiime demux summarize \
  --i-data paired_end_demux_${i}.qza \
  --o-visualization paired_end_demux_${i}.qzv \
  --verbose

done
