library(plyr)
intable <- read.table("blastoutput.txt",header=TRUE)
counts <- ddply(intable, .(intable$seqid, intable$scaffold_length), nrow)
names(counts) <- c("seqid", "scaffold_length", "Freq")
write.table(counts, "selfhit_frequency_summary.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)


qseqcounts <- unique(intable[,c("seqid","qstart","qend","scaffold_length")])
names(qseqcounts) <- c("seqid","start","end","scaffold_length")
sseqcounts <- unique(intable[,c("seqid","sstart","send","scaffold_length")])
names(sseqcounts) <- c("seqid","start","end","scaffold_length")
seqcounts <- rbind(qseqcounts,sseqcounts)
sumseqcounts <- unique(seqcounts[,c("seqid","start","end","scaffold_length")])

write.table(sumseqcounts, "list_of_selfhit_regions.txt",quote=FALSE, col.names=TRUE,row.names=FALSE)

towrite <- ddply(sumseqcounts, .(sumseqcounts$seqid, sumseqcounts$scaffold_length), nrow)
names(towrite) <- c("seqid","scaffold_length","no_selfhit_regions")

write.table(towrite,"frequency_of_selfhit_regions_by_scaffold.txt",quote=FALSE, col.names=TRUE,row.names=FALSE)


