#Basic R
#www.pavannichani.com

#The purpose of this script is to learn how to play around with variables and data in R 

count<-0
primes<-c(1,3,5,7,11)
names<-c("Bob","Ted","Joe","Alisha")
Truth<-c(TRUE,FALSE)

organism<-c("Human","Chimpanzee","Yeast")
chromosomes<-c(23,24,16)
multicellular<-c(TRUE,TRUE,FALSE)

OrganismTable<- data.frame(organism,chromosomes,multicellular)

write.table(OrganismTable,file = "Basic.csv",row.names = FALSE, na="", col.names = FALSE, sep = ",")

NewDataFrame<-read.csv("Basic.csv", header=FALSE, sep=",")

count<-count*count
count<-count-1
count<-count/5

barplot(OrganismTable$chromosomes)
