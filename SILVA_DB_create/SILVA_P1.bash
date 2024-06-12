#The SILVA database used was 138.1
#https://www.arb-silva.de/documentation/release-1381/
#https://www.arb-silva.de/no_cache/download/archive/release_138.1/Exports/
#SSURef_NR99 or LSU for 23S
#eg - SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz

#Need to convert the database to get the fasta and taxonomy files as QIIME2 requires for the RDP classifier
gunzip SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
#number of entries
grep '^>' SILVA_138.1_SSURef_NR99_tax_silva.fasta | wc -l #510508
#some 'U' bases and need to be changed to T
sed -i '/^[A-Z]/s/U/T/g' SILVA_138.1_SSURef_NR99_tax_silva.fasta
#transforming to single line sequence - this makes using head etc easier
#If using '#' instead of #NEWLINE#, some tax lineages have #'s in them and it causes there to be line breaks when there should not be
sed -i '/^>/s/$/#NEWLINE#/' SILVA_138.1_SSURef_NR99_tax_silva.fasta
sed -i '/^>/s/^/START/' SILVA_138.1_SSURef_NR99_tax_silva.fasta
tr -d '\n' < SILVA_138.1_SSURef_NR99_tax_silva.fasta > SILVA_138.1_SSURef_NR99.fasta
sed -i 's/#NEWLINE#/\n/g' SILVA_138.1_SSURef_NR99.fasta
sed -i '/[A-Z]/s/START>/\n>/g' SILVA_138.1_SSURef_NR99.fasta

#remove any empty lines - like line 1
sed -i '/^$/d' SILVA_138.1_SSURef_NR99.fasta

#must match number of entries as before
grep '^>' SILVA_138.1_SSURef_NR99.fasta | wc -l
grep '^[A-Z]' SILVA_138.1_SSURef_NR99.fasta | wc -l

#If used '#' instead of newline - join it back to the preceding line
grep -v '^>\|^[A-Z]' SILVA_138.1_SSURef_NR99.fasta
grep -A 1 -B 1 -v '^>\|^[A-Z]' SILVA_138.1_SSURef_NR99.fasta
sed -i '/^[1-9]/s/^/###/g' SILVA_138.1_SSURef_NR99.fasta
sed -i -z 's/\n###//g' SILVA_138.1_SSURef_NR99.fasta

#Create tax file
grep '^>' SILVA_138.1_SSURef_NR99.fasta > SILVA_138.1_SSURef_NR99.tax
sed -i 's/ /\t/' SILVA_138.1_SSURef_NR99.tax
sed -i 's/^>//' SILVA_138.1_SSURef_NR99.tax
mv SILVA_138.1_SSURef_NR99.tax SILVA_138.1.tax

sed 's/ /_/g' SILVA_138.1.tax > SILVA_138.1_mod.tax
sed -i 's/;/\t/g' SILVA_138.1_mod.tax

wc -l SILVA_138.1_SSURef_NR99.tax

#rm tax from fasta file
sed -i 's/ .*//g' SILVA_138.1_SSURef_NR99.fasta
mv SILVA_138.1_SSURef_NR99.fasta SILVA_138.1.fasta

#taxo.tsv is the tax levels from SILVA
#https://www.arb-silva.de/no_cache/download/archive/release_138.1/Exports/taxonomy/
#tax_slv_ssu_138.1.txt.gz
gunzip tax_slv_ssu_138.1.txt.gz
mv tax_slv_ssu_138.1.txt taxo.tsv
