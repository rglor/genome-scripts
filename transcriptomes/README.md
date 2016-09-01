Loosely following the guidelines at:
http://sfg.stanford.edu/quality.html
https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running%20Trinity
http://pasapipeline.github.io/#A_ComprehensiveTranscriptome
http://informatics.fas.harvard.edu/best-practices-for-de-novo-transcritome-assembly-with-trinity.html


Step 1: concatenate all the files for R1 together, and do the same thing for R2

Step 2: run fastqc

Step 3: make a log directory
```
mkdir logs
```

Step 4: run cutadapt: discard 0 length reads and perform very gentle trimming (following http://journal.frontiersin.org/article/10.3389/fgene.2014.00013/full)
```
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 5 -m 1 -o Brain_trimmed_R1.fastq.gz -p Brain_trimmed_R2.fastq.gz Brain_ATCACG_R1.fastq.gz Brain_ATCACG_R1.fastq.gz > logs/cutadapt.log
```

Step 5: run trinity with strand specificity 

Trinity -seqType fq --max_memory 50G --left brain_pear.unassembled.forward.fastq --right brain_pear.unassembled.reverse.fastq --CPU 8 --full_cleanup --normalize_reads --min_kmer_cov 2 >> trinity.log

