#Step 1: run the onelining.R script on your original genome fasta file
 
#Step 2: line by line blast to itself. Make sure that linebyline.R is in the directory with your fasta file and tempout file.
#Change db to whatever yours is called in line 19
noseqs=`wc -l tempout | awk '{print $1}'`
echo "seqid" "pident" "match_length" "mismatch" "gapopen" "qstart" "qend" "sstart" "send" "evalue" "bitscore" "scaffold_length" > blastoutput.txt
for i in `seq 1 2 $noseqs`;
do tail -n+$i tempout | head -n1 > tempinseq;
j=`expr $i + 1`;
tail -n+$j tempout | head -n1 >> tempinseq;
~/bin/ncbi-blast-2.2.30+/bin/dustmasker -in tempinseq -out tempoutseq -outfmt fasta;
Rscript onelining_tempseq.R;
rm tempoutseq;
tail -n+2 tempseq | head -n1 | awk '{print length}' > length.txt;
echo "qseqid" "sqeqid" "pident" "match_length" "mismatch" "gapopen" "qstart" "qend" "sstart" "send" "evalue" "bitscore" > rawblast.txt;
~/bin/ncbi-blast-2.2.30+/bin/makeblastdb -in tempseq -dbtype nucl;
~/bin/ncbi-blast-2.2.30+/bin/blastn -db tempseq -query tempseq -perc_identity 75 -outfmt 6 -lcase_masking >> rawblast.txt;
Rscript linebyline.R >> linebylineblastlog.txt;
done;
