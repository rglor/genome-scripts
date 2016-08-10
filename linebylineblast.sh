#Step 1: run the onelining.R script on your original genome fasta file
 
#Step 2: make the blast database. Change the file name following -in to whatever your original fasta file was called.
~/bin/ncbi-blast-2.2.30+/bin/makeblastdb -in a.lines.fasta -dbtype nucl

#Step 3: line by line blast. Make sure that linebyline.R is in the directory with your fasta file and tempout file.
for i in `seq 1 $nosamples`;
do name=`tail -n+$i samplenames.txt | head -n1`;
