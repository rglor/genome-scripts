Loosely following the guidelines at:
http://sfg.stanford.edu/quality.html
https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running%20Trinity

Step 1: concatenate all the files for R1 together, and do the same thing for R2

Step 2: run fastqc

Step 3: make a log directory
```
mkdir logs
```

Step 4: run cutadapt
```
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o Brain_trimmed_R1.fastq.gz -p Brain_trimmed_R2.fastq.gz Brain_ATCACG_R1.fastq.gz Brain_ATCACG_R1.fastq.gz > logs/cutadapt.log
```

Step 5: run pear with quality filtering
