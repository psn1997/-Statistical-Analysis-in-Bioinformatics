#Multiple Sequence Alignment
#www.pavannichani.com

library(msa)
library(seqinr)
library(ape)
library(phangorn)

coxAA<-readAAStringSet("cox1multi.fasta")
prokDNA<-readDNAStringSet("prok.fasta")

#Aligns the cox1 protein sequences using CLUSTALW
coxAligned<-msa(coxAA)
#Aligns the cox1 DNA sequences using CLUSTALW
prokAligned<-msa(prokDNA)

print(prokAligned, show="complete")

msa(prokDNA, "ClustalW")
msa(prokDNA, "ClustalOmega")
msa(prokDNA, "Muscle")

#Store the sequence as a string set
prokAlignStr = as(prokAligned, "DNAStringSet")
#Store the sequence as a FASTa file
writeXStringSet(prokAlignStr, file="prokAligned.fasta")
write.phylip(prokAligned, "prokAligned.phylip")

#Phylogenetic Reconstruction
#Convert prokaryotic alignment to seqinr format
prokAligned2<-msaConvert(prokAligned, type="seqinr::alignment")

#Generate distance matrix using seqinr
prokdist<-dist.alignment(prokAligned2, "identity")

#Generate neighbour joining distance tree
prokTree<-nj(prokdist)
plot(prokTree)

#Maximum Parsimony Tree
prokAligned3<-msaConvert(prokAligned, type="phangorn::phyDat")
ParsTree<-pratchet(prokAligned3)
plot(ParsTree)

#Maximum Likelihood Tree
fit<-pml(prokTree, prokAligned3)
#Attempt to improve the tree using Jukes-Cantor Model
fitJC<-optim.pml(fit, model = "JC", rearrangement = "stochastic")
plot(fitJC)

#Bootstrap the optimized ML tree
#bootstrapped<-bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
#Plot the Bootstrapped tree
#plotBS(midpoint(fitJC$tree), bootstrapped, p = 50, type="p")


