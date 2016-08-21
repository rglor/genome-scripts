length <- as.matrix(read.table("length.txt"))
blast <- as.matrix(read.table("rawblast.txt"))
print(blast[2,1])
if(dim(blast)[1]>1) {
blast <- blast[-1,]
blast <- blast[(which(blast[,1]==blast[,2])),]
if(is.matrix(blast)) {
blast <- blast[(which(as.numeric(blast[,7])!=as.numeric(blast[,9]))),]
if(is.matrix(blast)) {
length <- matrix(length,ncol=1,nrow=(dim(blast)[1]))
blast <- cbind(blast,length)
blast <- blast[,-1]
write.table(blast, "blastoutput.txt",quote=FALSE, col.names=FALSE,row.names=FALSE, append=TRUE)
} else {
if(!(blast[7]==blast[9])) {
blast <- t(as.matrix(blast))
length <- matrix(length,ncol=1,nrow=(dim(blast)[1]))
blast <- cbind(blast,length)
blast <- blast[,-1]
blast <- t(as.matrix(blast))
write.table(blast, "blastoutput.txt",quote=FALSE, col.names=FALSE,row.names=FALSE, append=TRUE)
}
}
} else {
if(!(blast[7]==blast[9])) {
blast <- t(as.matrix(blast))
length <- matrix(length,ncol=1,nrow=(dim(blast)[1]))
blast <- cbind(blast,length)
blast <- blast[,-1]
blast <- t(as.matrix(blast))
write.table(blast, "blastoutput.txt",quote=FALSE, col.names=FALSE,row.names=FALSE, append=TRUE)
}
}
}
