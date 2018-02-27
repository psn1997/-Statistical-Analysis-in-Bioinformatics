#Linkage Analysis-GWAS
#www.pavannichani.com

library(GenABEL)

#Load the data
data<-load.gwaa.data(phenofile="phenotype.dat", genofile="genotype.raw")

#Generate a histogram for distribution of the continuous trait
ct<-phdata(data)$ct
hist(ct, col="slateblue", main="Distribution of ct")

#Generate a boxplot to see if ct depends on gender
boxplot(ct ~ phdata(data)$sex, col=c("blue", "red"), xlab="sex", names=c("F", "M"), main="Gender Distribution of ct")

#Do Quality Control Check on the data
qc <- check.marker(data, call = 0.95, perid.call = 0.95, maf = 1e-08, p.lev = 1e-08) 

#Create a new dataset containing only individuals and SNPs that passed QC
data.qc <- data[qc$idok, qc$snpok]

#Test to see whether or not ct depends on any of the markers
an <- qtscore(ct~1, data=data.qc, trait="gaussian")

#Plot the result in a "Manhatten" Plot
plot(an, col=c("olivedrab","slateblue"), cex=.5, main="Manhattan plot")

#Plot the genomic inflation corrected plot
plot(an, col=c("olivedrab","slateblue"), cex=.5, main="Manhattan plot", df="Pc1df")

#Add Bonferroni correction
pval.threshold <- 0.05
bonferroni <- -log10(pval.threshold / nids(data.qc))

#Replot with a line showing the Bonferroni score as a red line
plot(an, col=c("olivedrab","slateblue"), cex=.5, main="Manhattan plot", df="Pc1df")
abline(h=bonferroni, lty=3, color="red")