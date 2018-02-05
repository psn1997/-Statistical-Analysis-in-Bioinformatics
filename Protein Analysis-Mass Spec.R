#Protein Analysis-Mass Spec
#This Script will demonstrate several ways to analyse proteins
#www.pavannichani.com


library(gplots)
library(timeSeries)
library(MASS)
library(rgl)
library(ggplot2)

#Mass Spec Fragment Histogram
peptides.txt<-read.table("peptidefrags.txt", header=FALSE)
peptides<-as.vector(peptides.txt$V1)
hist(peptides,breaks=400) 

#Load Data & Convert it into Vector Form
mascot.txt<-read.table("mascot.txt", header=FALSE)
xtandem.txt<-read.table("xtandem.txt", header=FALSE)
protpro.txt<-read.table("protpro.txt", header=FALSE)
mascot<-as.vector(mascot.txt$V1)
xtandem<-as.vector(xtandem.txt$V1)
protpro<-as.vector(protpro.txt$V1)

#Merge data vectors in a single dataframe
combinedMSdata <- list(Mascot=mascot, XTandem=xtandem, ProtPro=protpro)
#Create Venn Diagram comparing protein ID's from both Mascot and xTandem
venn(combinedMSdata)

#LDA Analysis

#Load Data
Dataset<-read.csv("ms.csv", header=TRUE, na.strings="NA", dec=".", strip.white=TRUE)

#Limit to actual numeric data columns
RawData <- Dataset[,2:14]

filledcols = colSds(RawData) != 0.0
RawData <- RawData[,filledcols]

#The "Group" column contains an integer that sets the group for the test subject
test1.lda <- lda(Dataset$X1 ~ . , data=Dataset)
test1.lda.values <- predict(test1.lda,Dataset)

#Plotting the results
x <- test1.lda.values$x[,1]
y <- test1.lda.values$x[,2]

#Define the conditions as a seperate variable for use in plot and legend
class <- Dataset$X1

#Generate initial 2D plot data with centroids
plotdata <- data.frame(class, x, y)
centroids <- aggregate(cbind(x,y)~class,plotdata,mean)

#Calculating the distance between centroids
CentroidDistances <- dist(centroids, method = "euclidean", diag = TRUE, upper = FALSE, p = 2)
attr(CentroidDistances, "Labels") <- centroids$class

#Generate the plot
ggplot(plotdata,aes(x,y,color=factor(class))) + geom_point(size=3) + geom_point(data=centroids,size=7) + geom_text(data=centroids, size=7, label=centroids$class, colour="black") + ggtitle("LDA of Conditions 1-3") + geom_text(aes(label=Dataset$X1),hjust=0, vjust=0, colour="black")

#Ouput the Distance Matrix of Centroids as a file
write.csv(as.matrix(CentroidDistances, file = "centroiddistances.csv"))