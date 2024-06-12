#make a databases of interest using cutadapt
module load cutadapt/2.10-gimkl-2020a-Python-3.8.2

#algae_phylum.txt - which phyla are to be targeted - each phylum on a newline
#xxx.tax - ID*tab*taxonomyDetails - as used in DADA2
#yyy.fasta - >ID*newline*sequence - as used in DADA2

for algae_phylum in $(cat algae_phylum.txt); do grep $algae_phylum /path/to/Databases/xxx.tax | sed 's/\t.*//' > ${algae_phylum}_tax.txt; grep -A 1 -Fwf ${algae_phylum}_tax.txt /path/to/Databases/yyy.fasta | sed '/--/d' > ${algae_phylum}.fasta; done

#You could append the algae databases here or change the grep to search multiple phyla

#to test primers individually or in pairs against a database

#primer_pairs.txt - PairName_Forward_ReverseComplement
#Here the database sequence entries are in the forward direction

for pp in $(cat primer_pairs.txt); do Pair=$(echo $pp | sed 's/_.*//'); Forward=$(echo "$pp" | cut -d '_' -f 2); Reverse=$(echo $pp | sed 's/.*_//g'); cutadapt \
--discard-untrimmed \
--action=trim \
-e 0 \ #maximum error rate – no errors allowed between primer and sequence
--no-indels \ #no insertions or deletions allowed
-g "${Pair}_linked=$Forward;min_overlap=${#Forward};required...$Reverse;min_overlap=${#Reverse};required" \ #primer set to test – reverse primer input as reverse complement
-g "${Pair}_F_only=$Forward;min_overlap=${#Forward}"
-g "${Pair}_RevComp_only=$Reverse;min_overlap=${#Reverse}"
-o "${Pair}.fasta" \ #this provides a curated database based on the primer set provided
/path/to/algae/database.fasta; done > ${Pair}_cutadapt.slurm

#the slurm file provides stats on how many sequences matched the various primer variations
#To note, the database is investigated sequentially and thus if matched for the linked set, it cannot be matched for any other variation

#If wanting to use the curated database in software like DADA2, then using something like:
grep '^>' ${Pair}.fasta > temp.tax
sed -i 's/^>//' Paired_DB.tax
grep -Fwf Paired_DB.tax tax_file.txt > ${Pair}.tax
#can pull out the tax IDs for the curated database
