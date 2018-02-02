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

#Create a Dotplot using seqinr 
coxgenes<-read.fasta(file = "cox1multi.fasta", seqtype="AA")
cox1<-as.character(coxgenes[[1]])
cox2<-as.character(coxgenes[[2]])

#Simple Dotplot
dotPlot(cox1, cox2, main = "Human vs Mouse Cox1 Dotplot")
#Windowed Dotplot
dotPlot(cox1, cox2, wsize = 3, wstep = 3, nmatch = 3, main = "Human vs Mouse Cox1 Dotplot\nwsize = 3, wstep = 3, nmatch = 3")
#Dotplot of only the first 100 residues of the sequence
dotPlot(cox1[1:100], cox2[1:100], wsize = 3, wstep = 3, nmatch = 3, main = "Human vs Mouse Cox1 first 100 AA Dotplot\nwsize = 3, wstep = 3, nmatch = 3")

#Store the sequence as a string set
coxAlignStr = as(coxAligned, "AAStringSet")
#Store the sequence as a FASTa file
writeXStringSet(coxAlignStr, file="coxAligned.fasta")
write.phylip(coxAligned, "coxAligned.phylip")