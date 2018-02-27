#Protein Protein Interaction
#www.pavannichani.com

library(graph)
library(Rgraphviz)
library(RBGL)
library(yeastExpData)

#From the yeastDataExp library let's load the litG dataset
data("litG")

##Extract all the nodes from litG
litGnodes <- nodes(litG)

#Subgraphs, or individual connected components
connectedComponents <- connectedComp(litG)

#Pullout subgraph 3
component3 <- connectedComponents[[3]]
subgraph3 <- subGraph(component3, litG)

#plot subgraph 3
subgraph3plot <- layoutGraph(subgraph3, layoutType="neato")
renderGraph(subgraph3plot)

#Degree Plot
numdegrees <- graph::degree(litG)
numdegrees <- sort(numdegrees)
meandeg <- mean(numdegrees)
hist(numdegrees, col="red", main = paste("Degree Distribution - Protein Interactions in litG with a mean of " , meandeg))