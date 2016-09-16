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

Step 4: run cutadapt: discard reads <25 bp in length (otherwise Trinity will fail because its kmer = 25) and perform very gentle trimming (following http://journal.frontiersin.org/article/10.3389/fgene.2014.00013/full)
```
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 5 -m 25 -o Brain_trimmed_R1.fastq.gz -p Brain_trimmed_R2.fastq.gz Brain_ATCACG_R1.fastq.gz Brain_ATCACG_R2.fastq.gz > logs/cutadapt.log
```

Step 5: run trinity with strand specificity 
```
Trinity -seqType fq --max_memory 50G --left Brain_trimmed_R1.fastq.gz --right Brain_trimmed_R2.fastq.gz -SS_lib_type RF --CPU 8 --full_cleanup --normalize_reads --min_kmer_cov 2 >> trinity.log
```

Step 6A: QC of transcriptome assembly: Mapping reads back to transcriptome
```
bowtie2-build trinity_out_dir.Trinity.fasta Trinity.fasta

/public/bowtie2-2.2.9/bowtie2 --local --no-unal -x trinity.fasta -q -1 Brain_trimmed_R1.fastq.gz -2 Brain_trimmed_R2.fastq.gz | /public/samtools-1.3.1/samtools view -Sb | samtools sort -no - - > bowtied2.nameSorted.bam

 /public/trinityrnaseq-2.2.0/util/SAM_nameSorted_to_uniq_count_stats.pl bowtied2.nameSorted.bam 
```
From the Trinity instructions (link above): "A typical Trinity transcriptome assembly will have the vast majority of all reads mapping back to the assembly, and ~70-80% of the mapped fragments found mapped as proper pairs."

Step 6B: QC of transcriptome assembly:
```
/public/ncbi-blast-2.2.30+/bin/makeblastdb -in uniprot_sprot.fasta -dbtype prot

#This step takes forever btw
/public/ncbi-blast-2.2.30+/bin/blastx -query trinity_out_dir.Trinity.fasta -db uniprot_sprot.fasta -out blastx.outfmt6 -evalue 1e-20 -num_threads 6 -max_target_seqs 1 -outfmt 6
```

Step 6C: QC of transcriptome assembly: Using BUSCO to look at presence of conserved orthologs
```
 python3 /public/BUSCO_v1.22/BUSCO_v1.22.py -o busco_brain -in trinity_out_dir.Trinity.fasta -l /public/BUSCO_v1.22/buscolib/vertebrata/ -m trans
```

Step 6D: QC of transcriptome assembly: E90N50 transcript contig length 
```
#First need to estimate transcript abundance. Doing this using RSEM based on http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4842274/
/public/trinityrnaseq-2.2.0/util/align_and_estimate_abundance.pl --transcripts trinity_out_dir.Trinity.fasta --seqType fq --left Brain_trimmed_R1.fastq.gz --right Brain_trimmed_R2.fastq.gz --SS_lib_type RF --thread_count 4 --est_method RSEM --output_dir trin_rsem --aln_method bowtie --trinity_mode --prep_reference
```
