#FD v STOECK
#ensure that there is a complete overlap for FD sequences, they need to be trimmed to remove the last three characters

sed '/^[A-Z]/s/...$//' FD.fasta > FD_rm3char.fasta   

while IFS='' read -r line; do fastaID=$(echo $line | sed 's/^>.*/YES/' ); if [[ "$fastaID" == "YES" ]]; then ID=$(echo $line); else grep -B 1 $line Stoeck.fasta | grep '^>' | sed 's/$/, /' > match.txt; MATCH=$(cat match.txt | tr -d '\n'); printf "%s\t%s\n" "$ID" "$MATCH" >> FDvSTOECK.txt; fi; done < FD_rm3char.fasta

#FD v ZHAN
#FD lies within ZHAN and so there is no need to trim to assess the overlap

while IFS='' read -r line; do fastaID=$(echo $line | sed 's/^>.*/YES/' ); if [[ "$fastaID" == "YES" ]]; then ID=$(echo $line); else grep -B 1 $line Zhan.fasta | grep '^>' | sed 's/$/, /' > match.txt; MATCH=$(cat match.txt | tr -d '\n'); printf "%s\t%s\n" "$ID" "$MATCH" >> FDvZHAN.txt; fi; done < FD.fasta

#FD v ZIM
#FD lies within ZIM and so there is no need to trim to assess the overlap

while IFS='' read -r line; do fastaID=$(echo $line | sed 's/^>.*/YES/' ); if [[ "$fastaID" == "YES" ]]; then ID=$(echo $line); else grep -B 1 $line Zim.fasta | grep '^>' | sed 's/$/, /' > match.txt; MATCH=$(cat match.txt | tr -d '\n'); printf "%s\t%s\n" "$ID" "$MATCH" >> FDvZIM.txt; fi; done < FD.fasta
