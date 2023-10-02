#make a databases of interest

#algae_phylum.txt - which phyla are to be targeted - each phylum on a newline
#xxx.tax - ID*tab*taxonomyDetails - as used in DADA2
#yyy.fasta - >ID*newline*sequence - as used in DADA2

for algae_phylum in $(cat algae_phylum); do grep $algae_phylum /path/to/Databases/xxx.tax | sed 's/\t.*//' > ${algae_phylum}_tax.txt; grep -A 1 -Fwf ${algae_phylum}_tax.txt /path/to/Databases/yyy.fasta | sed '/--/d' > ${algae_phylum}_fasta.txt; done

#to test primers individually or in pairs against a database

#

for primer in $(cat primers2test.txt); do cutadapt --discard-untrimmed -e 0 --no-indels -g "$primer;min_overlap=${#primer}" -o xxx.fasta Ochrophyta.fasta >> heads.txt; done
-e 0 \ #maximum error rate – no errors allowed between primer and sequence
--no-indels \ #no insertions or deletions allowed
-g primer_sequence;min_overlap=no_primer_bp \ #primer to test –reverse primer input as reverse complement
-o output_file_name.fasta 
path_to_database

for pp in $(cat primer_pairs); do Pair=$(echo $pp | sed 's/_.*//'); Forward=$(echo "$pp" | cut -d '_' -f 2); Reverse=$(echo $pp | sed 's/.*_//g'); cutadapt --discard-untrimmed --action=trim -e 0 --no-indels -g "$Pair=$Forward;min_overlap=${#Forward};required...$Reverse;min_overlap=${#Reverse};required" -o "${Pair}.txt" Database.fasta; done
do Pair=$(echo $pp | sed 's/_.*//'); Forward=$(echo "$pp" | cut -d '_' -f 2); Reverse=$(echo $pp | sed 's/.*_//g'); DB_array=("Chl" "Rho" "Och"); i=0; len=${#DB_array[@]}; while [ $i -lt $len ]; do cutadapt --discard-untrimmed --action=trim -e 0 --no-indels -g "$Pair=$Forward;min_overlap=${#Forward};required...$Reverse;min_overlap=${#Reverse};required" -o "${Pair}.txt" "${DB_array[$i]}.fasta" > "${Pair}_${DB_array[$i]}_slurm.txt"; let i++; done; done
