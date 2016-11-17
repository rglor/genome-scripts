Loosely following the guidelines at:
http://sfg.stanford.edu/quality.html
https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running%20Trinity
http://pasapipeline.github.io/#A_ComprehensiveTranscriptome
http://informatics.fas.harvard.edu/best-practices-for-de-novo-transcritome-assembly-with-trinity.html


Step 1: Concatenate Fastq Files
======
If your project was run on different lanes or at different times, you must first use the `cat` command to concatenate fastq formatted sequences obtained each individual tissue sample. During this concatenation phase, you will need to create separate concatenated files for each end of your paired end reads (e.g., one concatenated file for R1 and another for R2). In the example, below, we are concatenating from reads from Illumina runs that were done at two separate times (7June2016 and 19Aug2016) using the same library.
```
cat ../anole_RNAseq_7June2016_run1/Project_Glor_Alexander/Sample_Digestiv/*R1* ../anole_RNAseq_19Aug2016_run2/Project_Glor_Alexander/Sample_Digestiv/*R1* >> Digestive_CTTGTA_R1.fastq.gz

cat ../anole_RNAseq_7June2016_run1/Project_Glor_Alexander/Sample_Digestiv/*R2* ../anole_RNAseq_19Aug2016_run2/Project_Glor_Alexander/Sample_Digestiv/*R2* >> Digestive_CTTGTA_R2.fastq.gz

```

Step 2: Preliminary QC & Quality Trimming
======
Preliminary QC of your sequences can be completed by applying the [`fastqc`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to your fastq sequence files. [`fastqc`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a popular tool that "aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines." Carefully inspect output from [`fastqc`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) because this is the main way that you are going to make sure that nothing is seriously wrong with your data before you delve into the time consuming series of analyses discussed below.
```
fastqc -k 6 *
```

If your sequences look OK after preliminary QC, its time to get your sequences ready for downstream analyses by trimming adaptors and eliminating low quality sequences and basecalls. Before going any further, you should start a log for subsequent output.
```
mkdir logs
```

Next use the function `cutadapt` to (1) trim Illumina adapter sequences, (2) discard reads <25 bp in length (otherwise *de novo* assembly in Trinity will fail because its kmer = 25) and (3) perform gentle trimming of low quality basecalls. Recent studies suggest that trimming to a relatively low PHRED score of 5 results in transcriptomes that are considerably more complete that those that result from more aggressive quality trimming, without a commensurate increase in errors ([MacManes 2014](http://journal.frontiersin.org/article/10.3389/fgene.2014.00013/full)).
```
#PBS -N cutadapt.sh
#PBS -l nodes=1:ppn=1:avx,mem=16000m,walltime=5:00:00
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome_RNAseq/Testes
#PBS -j oe
#PBS -o cutadapterror

cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 5 -m 25 -o Testes_trimmed_R1.fastq.gz -p Testes_trimmed_R2.fastq.gz Testes_CTTGTA_R1.fastq.gz Testes_CTTGTA_R2.fastq.gz > logs/cutadapt.log
```
After trimming is complete, use `fastqc` on the resulting files to check that adapters have been trimmed and that the newly generated `fastq` files look good. 

Step 3: *de novo* Assembly with Trinity
======
After you've conducted the basic QC steps described above, you're ready to do your fist *de novo* assembly. We use the `Trinity` package for *de novo* transcriptome assembly. Because *de novo* assembly with `Trinity` requires a large amount of computer memory (~1GB RAM/1 million sequence reads) and generates a large number of files (~900,000) you should be sure that your quotas are able to accommodate these requirements prior to initiating assembly. In the script below, we will run *de novo* assembly using a high memory node via the `bigm` queue; moreover, hundreds of thousands of files that are temporarily required during the assembly process are stored on the node running the analyses (`file=200gb') rather than in the scratch space, which has stricter constraints on file size and number. Importantly, all of the temporary files generated during the course of the `Trinity` analyses are deleted at the end of the script to avoid gumming up the cluster's hard drives.

```
#PBS -N trinity_heart
#PBS -q bigm -l nodes=1:ppn=24:avx,mem=512000m,walltime=72:00:00,file=200gb
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome_RNAseq/Heart
#PBS -j oe
#PBS -o trinity_heart_error

work_dir=$(mktemp -d) #Generates name for working directory where temporary files will be stored.
mkdir $work_dir #Generates working directory for temporary file storage.
cp /scratch/glor_lab/rich/distichus_genome_RNAseq/Heart/Heart_trimmed* $work_dir #Copies sequence files for assembly to new working directory on node doing analyses.
Trinity -seqType fq --max_memory 512G --left $work_dir/Heart_trimmed_R1.fastq.gz --right $work_dir/Heart_trimmed_R2.fastq.gz -SS_lib_type RF --CPU 24 --full_cleanup --normalize_reads --min_kmer_cov 2 --output $work_dir/trinity_out_dir >> trinity.log
mv $work_dir/trinity_out_dir.Trinity.fasta /scratch/glor_lab/rich/distichus_genome_RNAseq/Heart/ #Moves files from node doing analyses to scratch drive.
rm -rf $work_dir #Deletes working directory. This step is necessary because otherwise you will really gum up the node.

```

Step 4: QC of Transcriptome
======
We conduct several separate QC steps with the Trinity assembly, most of which are directly recommended by the authors of Trinity (the commands below are mostly copied directly from the Trinity github).

Step 4a: Mapping Reads to Assembled Transcriptome
------
Our first assembly QC step will involve mapping our individual sequencing reads back to the assembly. This step is necessary to ensure that most reads (70-80%) can be mapped to the assembly as "proper pairs;" those that do not are likely the result of extremely low frequency, erroneous sequences, or instances where one of the two paired end sequences failed. In order to complete this operation, we must have access to `bowtie` and `samtools`, both of which are available via the `env-selector-menu` option on the cluster. Once all the appropriate software is available, we need to build a bowtie2 index. This operation will produce six files with the extension `bt2` and should take 10-15 minutes run on a single processor in interactive mode.
```
bowtie2-build trinity_out_dir.Trinity.fasta Trinity_testes.fasta
```
Next we do the actual mapping, which is best submitted to the default queue because this function may take a day or more to complete. I don't know how parallelized this process is but assign it to 24 processors.
```
#PBS -N mapping_testes
#PBS -q default -l nodes=1:ppn=24:avx,mem=50000m,walltime=24:00:00
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome_RNAseq/Testes
#PBS -j oe
#PBS -o mapping_testes_error
bowtie2 --local --no-unal -x trinity.fasta -q -1 Testes_trimmed_R1.fastq.gz -2 Testes_trimmed_R2.fastq.gz | samtools view -Sb - | samtools sort -no - - > bowtie2.nameSorted.bam
```
To get a summary of the mapping results, we can use a perl script from Trinity.
```
SAM_nameSorted_to_uniq_count_stats.pl bowtie2.nameSorted.bam

#read_type      count   pct
proper_pairs    30446172        94.39
improper_pairs  928085  2.88
left_only       699683  2.17
right_only      181046  0.56

Total aligned rnaseq fragments: 32254986
```

Step 4b: Transcriptome Contig Statistics (Nx  and ExNy)
------
In this step, we calculate basic statistics to assess contig assembly. Nx is the transcriptome equivalent of N50. 

```
TrinityStats.pl trinity_out_dir.Trinity.fasta > transcriptome_contig_nx_stat.log
```
The authors of Trinity encourage assessment of Nx based on the single longest transcript per isoform given that the tendency of assemblers to produce too many isoforms, particularly for larger transcripts (as a result, we paradoxically expect that our N50 for the dataset with only the longest isoform will actually be lower than that for all transcripts).

```
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':  456614
Total trinity transcripts:      590779
Percent GC: 42.48

########################################
Stats based on ALL transcript contigs:
########################################

        Contig N10: 3701
        Contig N20: 2409
        Contig N30: 1623
        Contig N40: 1072
        Contig N50: 703

        Median contig length: 304
        Average contig: 535.45
        Total assembled bases: 316335067


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

        Contig N10: 3065
        Contig N20: 1523
        Contig N30: 837
        Contig N40: 561
        Contig N50: 426

        Median contig length: 283
        Average contig: 425.57
        Total assembled bases: 194318977
```

Trinity's authors encourage use of the Ex statistic, a variant of Nx that incoporporations information on transcript frequency. In order to do this, we must first estimate transcript abundance using 'RSEM' ([Teng et al. 2016](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4842274/). We will be using `RSEM` via a script in the `Trinity` toolkit, but still need to be sure to place `RSEM` in our path before proceeding. You can run the following command to prepare the reference and run alignment and abundance estimation in interactive mode; it should only take 10 or 15 minutes to run, but will generate at least 1GB of files with `trinity_out_dir.Trinity.fasta.bowtie` or `trinity_out_dir.Trinity.fasta.RSEM` prefixes.

```
align_and_estimate_abundance.pl --transcripts trinity_out_dir.Trinity.fasta --seqType fq --left Testes_trimmed_R1.fastq.gz --right Testes_trimmed_R2.fastq.gz --SS_lib_type RF --thread_count 24 --est_method RSEM --output_dir trin_rsem --aln_method bowtie --trinity_mode --prep_reference
```
We then construct matrices containing counts and a matrix of normalized expression values for isoforms (with `trans_counts` prefix)  and genes (with `gene_counts` prefix). Normalized values will be calculated as transcripts per million transcripts (TPM). Running these perl scripts included in the `Trinity` toolkit should only take a few seconds in interactive mode.
```
abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix trans_counts --name_sample_by_basedir trin_rsem/RSEM.isoforms.results
abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix gene_counts --name_sample_by_basedir trin_rsem/RSEM.genes.results
```
Next we are going to count numbers of expressed transcripts or genes. Note that the required function is part of the `Trinity` toolkit, but the folder its in does not automatically install in your path when you load `Trinity` explaining the explicit directory information below (if you don't do this or ad the `misc` folder in the `Trinity` `util` folder to your path you will get an error).
```
/tools/cluster/6.2/trinityrnaseq/2.2.0/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl trans_counts.TPM.not_cross_norm > trans_matrix.TPM.not_cross_norm.counts_by_min_TPM

/tools/cluster/6.2/trinityrnaseq/2.2.0/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl gene_counts.TPM.not_cross_norm > genes_matrix.TPM.not_cross_norm.counts_by_min_TPM
```
Follow the guidance and R-script at: https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification. Inside the trin_rsem folder, run the following commands to extract the columns of interest (because we are only working on one sample)
```
cat RSEM.isoforms.results  | perl -lane 'print "$F[0]\t$F[5]";' >  RSEM.isoforms.results.mini_matrix
cat RSEM.genes.results  | perl -lane 'print "$F[0]\t$F[5]";' >  RSEM.genes.results.mini_matrix
```
Then move back to the folder containing the trinity.fasta file and run the following commands.
```
/tools/cluster/6.2/trinityrnaseq/2.2.0/util/misc/contig_ExN50_statistic.pl trin_rsem/RSEM.isoforms.results.mini_matrix trinity_out_dir.Trinity.fasta > ExN50_trans.stats

/tools/cluster/6.2/trinityrnaseq/2.2.0/util/misc/contig_ExN50_statistic.pl trin_rsem/RSEM.genes.results.mini_matrix trinity_out_dir.Trinity.fasta > ExN50_genes.stats
```

Step 4c: Assess Full-length Transcripts Relative to Reference Via BLAST+
------
```
/public/ncbi-blast-2.2.30+/bin/makeblastdb -in uniprot_sprot.fasta -dbtype prot

/public/ncbi-blast-2.2.30+/bin/blastx -query trinity_out_dir.Trinity.fasta -db uniprot_sprot.fasta -out blastx.outfmt6 -evalue 1e-20 -num_threads 6 -max_target_seqs 1 -outfmt 6

/public/trinityrnaseq-2.2.0/util/analyze_blastPlus_topHit_coverage.pl blastx.outfmt6 trinity_out_dir.Trinity.fasta uniprot_sprot.fast > analyze_blastPlus_topHit_coverage.log

/public/trinityrnaseq-2.2.0/util/misc/blast_outfmt6_group_segments.pl blastx.outfmt6 trinity_out_dir.Trinity.fasta uniprot_sprot.fasta > blast.outfmt6.grouped

/public/trinityrnaseq-2.2.0/util/misc/blast_outfmt6_group_segments.tophit_coverage.pl blast.outfmt6.grouped > analyze_groupsegments_topHit_coverage.log
```

Step 4d: Assess Presence of Conserved Orthologs Via BUSCO
------
```
 python3 /public/BUSCO_v1.22/BUSCO_v1.22.py -o busco_brain -in trinity_out_dir.Trinity.fasta -l /public/BUSCO_v1.22/buscolib/vertebrata/ -m trans
```

Step 4e: Compute DETONATE scores
------
```
#RSEM (reference-free mode)
/public/detonate-1.11-precompiled/rsem-eval/rsem-eval-estimate-transcript-length-distribution trinity_out_dir.Trinity.fasta /public/detonate-1.11-precompiled/rsem-eval/true_transcript_length_distribution/anolis_distichus.txt

/public/detonate-1.11-precompiled/rsem-eval/rsem-eval-calculate-score --paired-end --bam trin_rsem/bowtie.bam trinity_out_dir.Trinity.fasta sample_brain 200 --transcript-length-parameters /public/detonate-1.11-precompiled/rsem-eval/true_transcript_length_distribution/anolis_distichus.txt --strand-specific -p 4  >& rsem_eval.log
```

Step 4f: QC of transcriptome assembly: TransRate (as of the 15-Oct-2016, transrate has some issues as mentioned in https://github.com/blahah/transrate/issues/201, so running both v1.0.1 and v.1.0.3 and splicing together their output.
------
```transrate --assembly trinity_out_dir.Trinity.fasta --left Brain_trimmed_R1.fastq.gz --right Brain_trimmed_R2.fastq.gz --threads 8 >> transrate103.log

mv transrate_results transrate_results_103

transrate _1.0.1_ --assembly trinity_out_dir.Trinity.fasta --left Brain_trimmed_R1.fastq.gz --right Brain_trimmed_R2.fastq.gz --threads 8 >> transrate101.log
```

