#Pairwise Sequence Alignment
#www.pavannichani.com


library(Biostrings)
library(seqinr)

prokaryotes<-read.fasta(file = "prok.fasta", seqtype = "DNA")

seq1<-as.character(prokaryotes[[1]])
seq1<-paste(seq1,collapse = "")

seq2<-as.character(prokaryotes[[2]])
seq2=paste(seq2, collapse = "")

#Align seq1 and seq2 using the default settings of Biostrings
pairalign<-pairwiseAlignment(pattern = seq2, subject = seq1)

#Convert Alignment to FASTA and save
pairalignString = BStringSet( c( toString( subject(pairalign) ), toString(pattern(pairalign))))
writeXStringSet(pairalignString, "aligned.txt", format="FASTA")