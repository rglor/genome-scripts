length <- as.matrix(read.table("length.txt"))
blast <- as.matrix(read.table("rawblast.txt"))
blast <- blast[(which(blast[,1]==blast[,2])),]
blast <- blast[(which(blast[,7]!=blast[,9])),]
length <- matrix(length,ncol=1,nrow=(dim(blast)[1]))
blast <- cbind(blast,length)
blast <- blast[,-1]
write.table(blast, "blastoutput.txt",quote=FALSE, col.names=FALSE,row.names=FALSE, append=TRUE)
