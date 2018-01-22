#Script to work with FASTA files
#www.pavannichani.com

library(seqinr)
library(ape)
library(rentrez)

#Read a local FASTA file
cox1<-read.fasta(file = "cox1.fasta",seqtype = "AA")

#Retrieve a GENBANK sequence as a binary object
AB003468<-read.GenBank("AB003468", as.character = "TRUE")

#Save GENBANK Sequence as FASTA format
write.dna(AB003468, file ="AB003468.fasta", format = "fasta", append = FALSE, nbcol = 6, colsep = " ", colw = 10)

#Retrieve sequences
entrez_search(db="nucleotide", term="human superoxide dismutase")

#Word Count Analysis
#Covert sequence into simple string
CloningVector<-AB003468[[1]]

#Nucleotide count
count<-count(CloningVector,1) 

#We can just as easily see all dinucleotide combinations by using a 2, or trinucleotide combinations by using a 3
#count(CloningVector,2)
#count(CloningVector,3)

#GC content tells us what the ratio of G (guanines) and C (cystosines) residues compared to A(adenine) and T (thymidine) residues in a sequence is.
#This is important because coding regions tend to be higher in GC.
#GC content also affects the "melting" temperature of DNA 

GC<-GC(CloningVector)

#We’ll use a variable called GCwindow to break the sequence into chunks of 200
GCwindow<-seq(1, length(CloningVector)-200, by = 200)

#For Loop to compute GC per chunk
n<-length(GCwindow)
Chunks <- numeric(n)
for(i in 1:n){
  chunk<-CloningVector[GCwindow[i]:(GCwindow[i]+199)]
  chunkGC<-GC(chunk)
  print(chunkGC)
  Chunks[i]<-chunkGC
}

plot(GCwindow,Chunks,type="b",xlab="Nucleotide start position",ylab="GC content")

#Custom Function for GC Window Plotting
slidingwindowGCplot<-function(windowsize,inputseq)
{
  GCwindow<-seq(1, length(inputseq)-windowsize, by = windowsize)
  n<-length(GCwindow)
  Chunks <- numeric(n)
  for(i in 1:n){
    chunk<-inputseq[GCwindow[i]:(GCwindow[i]+windowsize-1)]
    chunkGC<-GC(chunk)
    print(chunkGC)
    Chunks[i]<-chunkGC
  }
  plot(GCwindow,Chunks,type="b",xlab="Nucleotide start position",ylab="GC content",main=paste("GC Plot with windowsize ", windowsize))
}
slidingwindowGCplot(200,CloningVector)
slidingwindowGCplot(175,CloningVector)
slidingwindowGCplot(100,CloningVector)
