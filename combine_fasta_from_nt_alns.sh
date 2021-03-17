mkdir -p aga-alignments

for FILE in *_aln.fasta ;
do
	cp $FILE aga-alignments
	seqkit grep -v -r -p "NC_045512.2" $FILE > $FILE%1
	awk '/^>/ {gsub(/.fa(sta)?$/,"",FILENAME);printf(">%s\n",FILENAME);next;} {print}' $FILE%1 > renamed_$FILE
	awk -F- '/^>/{print $1;next}{print}' renamed_$FILE > renamed2_$FILE
	rm $FILE%1
	mv renamed2_$FILE $FILE
	#rm renamed2_$FILE
	rm renamed_$FILE
done
cat *_aln.fasta > run_consensus_aln.fasta
sed -i '' "s/\./N/g" run_consensus_aln.fasta
sed -i '' "s/_alnNfasta%1//g" run_consensus_aln.fasta 



