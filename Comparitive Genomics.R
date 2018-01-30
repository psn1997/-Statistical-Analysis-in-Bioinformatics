#Comparative Genomics
#www.pavannichani.com

library(biomaRt)
ensembl=useMart("ensembl")

listDatasets(ensembl)
# Get Datasets for Chimp and Gorilla
Chimpanzee <- useDataset("ptroglodytes_gene_ensembl", mart = ensembl)
Gorilla <- useDataset("ggorilla_gene_ensembl", mart = ensembl)

# Retrieve Attribute Lists
ChimpanzeeAttributes <- listAttributes(Chimpanzee)
GorillaAttributes <- listAttributes(Gorilla)

# Retrieve Ensemble GENE ID's
ChimpanzeeGenes <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"), mart=Chimpanzee)
GorillaGenes <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"), mart=Gorilla)

# Parse Data into names, types vectors
ChimpanzeeGeneNames <- ChimpanzeeGenes[[1]]
ChimpanzeeGeneTypes <- ChimpanzeeGenes[[2]]
GorillaGeneNames <- GorillaGenes[[1]]
GorillaGeneTypes <- GorillaGenes[[2]]

# Create Table of Gene Types
ChimpanzeeTypesTable <- table(ChimpanzeeGeneTypes)
GorillaTypesTable <- table(GorillaGeneTypes)

# Count Protein Coding Genes
ChimpanzeeProteins <- ChimpanzeeTypesTable["protein_coding"]
GorillaProteins <- GorillaTypesTable["protein_coding"]

ChimpGGNum <- getBM(attributes = "ggorilla_homolog_ensembl_gene", mart = Chimpanzee)
ChimpanzeeGG <- getBM(attributes = c("ensembl_gene_id","ens_hs_gene"), mart = Chimpanzee)

GorillaCNum <- getBM(attributes = "ptroglodytes_homolog_ensembl_gene", mart = Gorilla)
GorillaC <- getBM(attributes = c("ensembl_gene_id","ptroglodytes_homolog_ensembl_gene"), mart = Gorilla)
