#!/bin/bash
#SBATCH -J Q2_PF_P3
#SBATCH -A xxx
#SBATCH --time=24:00:00
#SBATCH --mem=60GB
#SBATCH --cpus-per-task=8

module load QIIME2/2021.2

export TMPDIR=xxx/tmp_${SLURM_JOB_ID}
mkdir -p $TMPDIR
export TMPDIR

####NOTES####
#Again, this could be made into a loop - calling the primer and the respective database.
#This is the generic code used, and at the time was used independently and not in a loop.

#Primer  Forward  Reverse
#FD  TTAAAAAGCKCGTAGTTG  ACTTTCGTTCTTGAT
#Zim  ATTCCAGCTCCAATAGCG  GACTACGATGGTATCTAATC
#Zhan  AGGGCAAKYCTGGTGCCAGC  GRCGGTATCTRATCGYCTT
#Stoeck  CCAGCASCYGCGGTAATTCC  ACTTTCGTTCTTGATYRA
#Sherwood  GGACAGAAAGACCCTATGAA  TCAGCCTGTTATCCCTAGAG
#Herlemann  CCTACGGGNGGCWGCAG  GACTACHVGGGTATCTAATCC

#SILVA databases used SSURef_NR99
####

#Training a Naive Bayes classifier - assigning taxonomy from different databases
# Importing the fasta and taxonomy files
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /full_path/SILVA_138.1.fasta \
  --output-path rdp_seq.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path /full_path/SILVA_138.1.tax \
  --output-path rdp_ref_taxonomy.qza

# Extracting the sequence portion we have targeted
qiime feature-classifier extract-reads \
  --i-sequences rdp_seq.qza \
  --p-f-primer FORWARD \
  --p-r-primer REVERSE \
  --o-reads rdp_seq.qza

# Trainning the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads rdp_seq.qza \
  --i-reference-taxonomy rdp_ref_taxonomy.qza \
  --o-classifier rdp_classifier.qza \
  --verbose

#Assigning taxonomy with the trainned classifier and looking at the results
qiime feature-classifier classify-sklearn \
  --i-classifier rdp_classifier.qza \
  --i-reads rep_seqs_dada2.qza \
  --p-confidence 0.97 \
  --p-n-jobs $SLURM_CPUS_PER_TASK \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# Exporting table
qiime tools export \
   --input-path taxonomy.qza \
   --output-path ./
