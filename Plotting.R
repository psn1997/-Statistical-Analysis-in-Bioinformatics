#Plotting graph using ggplot2
#www.pavannichani.com

library(ggplot2)
library(reshape2)

rawdata<-read.csv("Plotdata.csv", header = TRUE)
melted=melt(rawdata,id.vars ="Subject", measure.vars = c("a","c","d","e","f","g","j","k"))

ggplot(rawdata, aes(x=Subject, y=a)) + geom_point()
myPlot<-ggplot(melted, aes(x=variable, y=value, col=Subject, group=Subject)) + geom_point() + geom_line() + xlab("Sample") + ylab("# Observed") + ggtitle("Few Observations")

myPlot