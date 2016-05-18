library(ggbio)
library(BSgenome.Athaliana.TAIR.TAIR9) #arabidopsis genome from 2009
library(org.At.tair.db) #arabidopsis annotation 2015
library("AnnotationHub")
library("biomaRt")
library("TxDb.Athaliana.BioMart.plantsmart22")
#from https://www.bioconductor.org/help/workflows/annotation/annotation/
ah <- AnnotationHub() #create local annotation hub object
unique(ah$species)#see all species
arabidopsis <- subset(ah, ah$species == "Arabidopsis thaliana")
#there are two At records, choose one?
grs <- query(arabidopsis, "GRanges") #look for one with GRanges - neither have this grrrrr
myat <- arabidopsis[["AH49574"]] #piked the one from ncbi
myat #see what we've got
columns(myat)#all types of data that can be retrieved from this object
keytypes(myat) #the types of data that can be used as keys
head(keys(myat, keytype="TAIR")) #AGIs
head(keys(myat, keytype="SYMBOL")) #abbrev names of genes
keys(myat, keytype="SYMBOL", pattern="COX") # a way to pick out a type of class of genes using keys
keys(myat, keytype="ENTREZID", pattern="COX", column="SYMBOL") #gives back the entrezid instead of symbol


genome <- BSgenome.Athaliana.TAIR.TAIR9
genome$Chr2[16672493:16675748] #coi1
autoplot(Athaliana, "COI1")
autoplot(genome, which = "genome$Chr2[16672493:16675748]", geom = NULL)
