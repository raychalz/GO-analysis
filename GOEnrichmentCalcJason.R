####GO category enrichment
library("org.At.tair.db")
library("topGO")
library("goProfiles")


GeneLists <- read.csv("GeneList.csv")
GeneUniverse <- vector(levels(GeneLists$AllGenes))
CamDiffGenes <- levels(GeneLists$CamTot.HEM)[-1]
LesDiffGenes <- levels(GeneLists$LesTot.HEM)[-1]

camGeneList <- factor(as.integer(GeneUniverse %in% CamDiffGenes))
names(camGeneList) <- GeneUniverse
lesGeneList <- factor(as.integer(GeneUniverse %in% LesDiffGenes))
names(lesGeneList) <- GeneUniverse

camGOdata <- new("topGOdata", ontology = "BP", allGenes = camGeneList, geneSel = function(p) p == 1, 
              description = "Test", annot = annFUN.org, mapping="org.At.tair.db")
lesGOdata <- new("topGOdata", ontology = "BP", allGenes = lesGeneList, geneSel = function(p) p == 1, 
                 description = "Test", annot = annFUN.org, mapping="org.At.tair.db")

camResultFisher <- runTest(camGOdata, algorithm = "classic", statistic = "fisher")
camGenTable <- GenTable(camGOdata,classicFisher = camResultFisher, orderBy="classicFisher", ranksOf="classicFisher", topNodes = 500)
lesResultFisher <- runTest(lesGOdata, algorithm = "classic", statistic = "fisher")
lesGenTable <- GenTable(lesGOdata,classicFisher = lesResultFisher, orderBy="classicFisher", ranksOf="classicFisher", topNodes = 500)

write.csv(camGenTable, "CamBPEnrichment.csv")
write.csv(lesGenTable, "LesBPEnrichment.csv")




######################
###Networks specific GO categories

data.net <- read.csv("GeneNetworkAdherence.csv")

networkResults <- list()
counter <- 1

for(i in unique(data.net$Phenotype)){
  for(j in unique(data.net$Network)){
        
    network.vec <- data.net$AGI[data.net$Phenotype == i & data.net$Network == j]
    networkGeneList <- factor(as.integer(GeneUniverse %in% network.vec))
    names(networkGeneList) <- GeneUniverse
    
    netGOdata <- new("topGOdata", ontology = "BP", allGenes = networkGeneList, geneSel = function(p) p == 1, 
                 description = "Test", annot = annFUN.org, mapping="org.At.tair.db")
    
    netResultFisher <- runTest(netGOdata, algorithm = "classic", statistic = "fisher")
    netGenTable <- GenTable(netGOdata,classicFisher = netResultFisher, orderBy="classicFisher", ranksOf="classicFisher", topNodes = 500)
    
    networkResults[[counter]] <- netGenTable
    names(networkResults)[counter] <- paste(i, j, sep = ".")
    
    counter <- counter + 1
    
  }
}


TotalNetworkResults <- do.call("rbind", networkResults)
write.csv(TotalNetworkResults, "NetworkGOEnrichment.csv")
