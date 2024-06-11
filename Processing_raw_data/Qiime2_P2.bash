#!/bin/bash
#SBATCH -J Q2_PF_P2
#SBATCH -A xxx
#SBATCH --time=10:00:00
#SBATCH --mem=40GB
#SBATCH --cpus-per-task=8

module load QIIME2/2021.2

export TMPDIR=xxx/tmp_${SLURM_JOB_ID}
mkdir -p $TMPDIR
export TMPDIR

####NOTES####
#After assessing the sequence quality and where it should be truncated, create a file that will loop through the primer datasets and truncate accordingly.
#trunc.tsv
#Herlemann#270#230
#FD#230#230
#Zim#230#230
#Stoeck#230#230
#Zhan#230#230
#Sherwood#230#230
############

for i in $(cat trunc.tsv); \
do \

Primer=$(echo ${i} | cut -f 1 -d '#'); \
echo $Primer; \
Fwd=$(echo ${i} | cut -f 2 -d '#'); \
echo $Fwd; \
Rev=$(echo ${i} | cut -f 3 -d '#'); \
echo $Rev; \

#Quality and chimera filtering using DADA2
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./paired_end_demux_${Primer}.qza \
  --p-trim-left-f 0 \
  --p-trunc-len-f ${Fwd} \
  --p-trim-left-r 0 \
  --p-trunc-len-r ${Rev} \
  --p-chimera-method consensus \
  --p-min-fold-parent-over-abundance 1 \
  --p-n-threads 0 \
  --o-representative-sequences rep_seqs_dada2_${Primer}.qza \
  --o-table table_dada2_${Primer}.qza \
  --o-denoising-stats Sriv-stats-dada2_${Primer}.qza \
  --verbose

#Looking at the feature table and feature data summary
qiime feature-table summarize \
  --i-table table_dada2_${Primer}.qza \
  --o-visualization table_dada2_${Primer}.qzv \

qiime feature-table tabulate-seqs \
  --i-data rep_seqs_dada2_${Primer}.qza \
  --o-visualization rep_seqs_dada2_${Primer}.qzv

#if no features/OTUs are present in the table this creates a problem trying to read in the tbale in R
qiime feature-table filter-samples \
 --i-table table_dada2_${Primer}.qza \
 --p-min-features 1 \
 --o-filtered-table rm_0_features_table_${Primer}.qza

# Exporting table
qiime tools export \
  --input-path ./rm_0_features_table_${Primer}.qza\
  --output-path ./

mv feature-table.biom ./feature-table_${Primer}.biom

#Exporting rep seqs
qiime tools export \
--input-path ./rep_seqs_dada2_${Primer}.qza \
--output-path ./

mv dna-sequences.fasta ./dna-sequences.fasta_${Primer}

#looking at denoising stats from dada2
  qiime metadata tabulate \
  --m-input-file Sriv-stats-dada2_${Primer}.qza \
  --o-visualization Sriv-stats-dada2_${Primer}.qzv

done
