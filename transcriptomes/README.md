Loosely following the guidelines at:
http://sfg.stanford.edu/quality.html
https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running%20Trinity
http://pasapipeline.github.io/#A_ComprehensiveTranscriptome
http://informatics.fas.harvard.edu/best-practices-for-de-novo-transcritome-assembly-with-trinity.html


Step 1 Concatenation: If your project was run on different lanes or at different times, you must fist use the cat command to concatenate  fastq formatted sequences obtained each individual tissue sample. During this concatenation phase, you will need to create separate concatenated files for each end of your paired end reads (e.g., one concatenated file for R1 and another for R2). In the example, below, we are concatenating from reads from Illumina runs that were done at two separate times (7June2016 and 19Aug2016) using the same library.
```
cat ../anole_RNAseq_7June2016_run1/Project_Glor_Alexander/Sample_Digestiv/*R1* ../anole_RNAseq_19Aug2016_run2/Project_Glor_Alexander/Sample_Digestiv/*R1* >> Digestive_CTTGTA_R1.fastq.gz

cat ../anole_RNAseq_7June2016_run1/Project_Glor_Alexander/Sample_Digestiv/*R2* ../anole_RNAseq_19Aug2016_run2/Project_Glor_Alexander/Sample_Digestiv/*R2* >> Digestive_CTTGTA_R2.fastq.gz

```

Step 2 Preliminary QC: Preliminary QC of your sequences can be completed by applying the fastqc function (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to your fastq sequence files.
```
fastqc -k 6 *
```

Step 3 Start Logging: Make a log directory to store subsequent output.
```
mkdir logs
```

Step 4A Trim Sequences: In this step, you will use the function cutadapt to (1) trim Illumina adapter sequences, (2) discard reads <25 bp in length (otherwise de novo assembly in Trinity will fail because its kmer = 25) and (3) perform gentle trimming of low quality basecalls. Recent studies suggest that trimming to a relatively low PHRED score of 5 results in transcriptomes that are considerably more complete that those that result from more aggressive quality trimming, without a commensurate increase in errors (http://journal.frontiersin.org/article/10.3389/fgene.2014.00013/full).
```
#PBS -N cutadapt.sh
#PBS -l nodes=1:ppn=1:avx,mem=16000m,walltime=5:00:00
#PBS -M alana.alexander@ku.edu
#PBS -m abe
#PBS -d /scratch/a499a400/anolis/transcriptome/Sample_Digestiv
#PBS -j oe
#PBS -o cutadapterror

cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 5 -m 25 -o Brain_trimmed_R1.fastq.gz -p Brain_trimmed_R2.fastq.gz Brain_ATCACG_R1.fastq.gz Brain_ATCACG_R2.fastq.gz > logs/cutadapt.log
```
Step 4B QC of Trimmed Reads: Run fastqc again on your trimmed sequences to make sure the seqeunces still look good and adaptors have been trimmed. 

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

Step 6B: QC of transcriptome assembly: Full-length transcript analysis for model and non-model organisms using BLAST+
```
/public/ncbi-blast-2.2.30+/bin/makeblastdb -in uniprot_sprot.fasta -dbtype prot

#This step takes forever btw
/public/ncbi-blast-2.2.30+/bin/blastx -query trinity_out_dir.Trinity.fasta -db uniprot_sprot.fasta -out blastx.outfmt6 -evalue 1e-20 -num_threads 6 -max_target_seqs 1 -outfmt 6

/public/trinityrnaseq-2.2.0/util/analyze_blastPlus_topHit_coverage.pl blastx.outfmt6 trinity_out_dir.Trinity.fasta uniprot_sprot.fast > analyze_blastPlus_topHit_coverage.log

/public/trinityrnaseq-2.2.0/util/misc/blast_outfmt6_group_segments.pl blastx.outfmt6 trinity_out_dir.Trinity.fasta uniprot_sprot.fasta > blast.outfmt6.grouped

/public/trinityrnaseq-2.2.0/util/misc/blast_outfmt6_group_segments.tophit_coverage.pl blast.outfmt6.grouped > analyze_groupsegments_topHit_coverage.log
```

Step 6C: QC of transcriptome assembly: Using BUSCO to look at presence of conserved orthologs
```
 python3 /public/BUSCO_v1.22/BUSCO_v1.22.py -o busco_brain -in trinity_out_dir.Trinity.fasta -l /public/BUSCO_v1.22/buscolib/vertebrata/ -m trans
```

Step 6D: QC of transcriptome assembly: The Transcriptome Contig Nx Statistic
```
/public/trinityrnaseq-2.2.0/util/TrinityStats.pl trinity_out_dir.Trinity.fasta > transcriptome_contig_nx_stat.log
```

Step 6E: QC of transcriptome assembly: Contig Ex90N50 Statistic and Ex90 Transcript Count
```
#First need to estimate transcript abundance. Doing this using RSEM based on http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4842274/
/public/trinityrnaseq-2.2.0/util/align_and_estimate_abundance.pl --transcripts trinity_out_dir.Trinity.fasta --seqType fq --left Brain_trimmed_R1.fastq.gz --right Brain_trimmed_R2.fastq.gz --SS_lib_type RF --thread_count 4 --est_method RSEM --output_dir trin_rsem --aln_method bowtie --trinity_mode --prep_reference

#Constructing a matrix of counts and a matrix of normalized expression values for isoforms and genes
/public/trinityrnaseq-2.2.0/util/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix trans_counts --name_sample_by_basedir trin_rsem/RSEM.isoforms.results

/public/trinityrnaseq-2.2.0/util/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix gene_counts --name_sample_by_basedir trin_rsem/RSEM.genes.results

#Counting numbers of expressed transcripts or genes
/public/trinityrnaseq-2.2.0/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl trans_counts.TPM.not_cross_norm > trans_matrix.TPM.not_cross_norm.counts_by_min_TPM

/public/trinityrnaseq-2.2.0/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl gene_counts.TPM.not_cross_norm > genes_matrix.TPM.not_cross_norm.counts_by_min_TPM

#Follow the guidance and R-script at: https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification

#Inside the trin_rsem folder, run the following commands to extract the columns of interest (because we are only working on one sample)
cat RSEM.isoforms.results  | perl -lane 'print "$F[0]\t$F[5]";' >  RSEM.isoforms.results.mini_matrix
cat RSEM.genes.results  | perl -lane 'print "$F[0]\t$F[5]";' >  RSEM.genes.results.mini_matrix

#cd up to where the trinity.fasta file is
/public/trinityrnaseq-2.2.0/util/misc/contig_ExN50_statistic.pl trin_rsem/RSEM.isoforms.results.mini_matrix trinity_out_dir.Trinity.fasta > ExN50_trans.stats

/public/trinityrnaseq-2.2.0/util/misc/contig_ExN50_statistic.pl trin_rsem/RSEM.genes.results.mini_matrix trinity_out_dir.Trinity.fasta > ExN50_genes.stats
```

Step 6F: QC of transcriptome assembly: Compute DETONATE scores
```
#RSEM (reference-free mode)
/public/detonate-1.11-precompiled/rsem-eval/rsem-eval-estimate-transcript-length-distribution trinity_out_dir.Trinity.fasta /public/detonate-1.11-precompiled/rsem-eval/true_transcript_length_distribution/anolis_distichus.txt

/public/detonate-1.11-precompiled/rsem-eval/rsem-eval-calculate-score --paired-end --bam trin_rsem/bowtie.bam trinity_out_dir.Trinity.fasta sample_brain 200 --transcript-length-parameters /public/detonate-1.11-precompiled/rsem-eval/true_transcript_length_distribution/anolis_distichus.txt --strand-specific -p 4  >& rsem_eval.log
```

Step 6G: QC of transcriptome assembly: TransRate (as of the 15-Oct-2016, transrate has some issues as mentioned in https://github.com/blahah/transrate/issues/201, so running both v1.0.1 and v.1.0.3 and splicing together their output.
```transrate --assembly trinity_out_dir.Trinity.fasta --left Brain_trimmed_R1.fastq.gz --right Brain_trimmed_R2.fastq.gz --threads 8 >> transrate103.log

mv transrate_results transrate_results_103

transrate _1.0.1_ --assembly trinity_out_dir.Trinity.fasta --left Brain_trimmed_R1.fastq.gz --right Brain_trimmed_R2.fastq.gz --threads 8 >> transrate101.log
```

