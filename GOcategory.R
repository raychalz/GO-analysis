library(goseq)
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/GWAS/Significant SNPs lists")


library(topGO); library(org.At.tair.db); library(goProfiles)
str(org.At.tair.db)
data <- read.csv("YlwLesionProportion.HEM.95Sig.genelist.csv")
genes <- as.character(data[,2]) #create character vector of AGIs of interest
genestest <- head(genes)

GOtest <- as.list(org.At.tairGENENAME[genestest]) #subset AGI list from Bimap interface
#GOtest is now a list of gene descriptions!

#=========================EXAMPLE FROM org.At.tair.pdf=========================================

# ## Bimap interface:
# x <- org.At.tairGENENAME
# # Get the TAIR identifiers that are mapped to a gene name
# mapped_tairs <- mappedkeys(x) #mapped_tairs is a list of AGIs (character vector?)
# # Convert to a list
# xx <- as.list(x[mapped_tairs])
# if(length(xx) > 0) {
#     # Get the GENENAME for the first five tairs
#     xx[1:5]
#     # Get the first one
#     xx[[1]]
# }
#=============================================================================================
#                                            {goProfiles}
#=============================================================================================

# basicProfile(genelist, idType = "Entrez", onto = "ANY", level = 2,orgPackage=NULL, anotPackage=NULL,
#              ord = TRUE, multilevels = NULL, empty.cats = TRUE, cat.names = TRUE, na.rm = TRUE)
aprof <- basicProfile(genestest, orgPackage= "org.At.tair.db", onto = "ANY")


# compareGeneLists(genelist1, genelist2, idType = "Entrez", onto = "ANY",
#                  level = 2, orgPackage,
#                  method = "lcombChisq", ab.approx = "asymptotic", confidence = 0.95, compareFunction="compareGOProfile
l <- compareGeneLists(genes[1:10], genes[11:20], onto = "BP", orgPackage= "org.At.tair.db")

compSummary(l, decs = 6) #returns brief summary of the comparison of two gene lists

expandedtest <- expandedProfile(genestest, orgPackage= "org.At.tair.db", onto = "BP")

fisherGOProfiles(expandedtest)#must enter object pn, object of class ExpandedGOProfile, GO Class-by-class Fisher tests in lists of genes characterized by their functional profiles

GOTermsList(genestest, orgPkg= "org.At.tair.db") #returns GOTerms for each AGI

# If several profiles have to be plot together they must be first merged using the ’mergeProfiles’ function.
#use mergeProfilesLists()
# plotProfiles(aProf, aTitle = "Functional Profile", anOnto = NULL, percentage = FALSE,
#              HORIZVERT = TRUE, legendText = NULL, colores = c("white", "red"), multiplePlots = F, multipleWindows 
plotProfiles(aprof)

prof1 <- basicProfile(genes[1:100], orgPackage= "org.At.tair.db", onto = "BP")
prof2 <- basicProfile(genes[101:200], orgPackage= "org.At.tair.db", onto = "BP")
twoprofiles.BP <-mergeProfilesLists(prof1, prof2, profNames=c("1:100", "101:200")) 
printProfiles(twoprofiles.BP) #prints table
plotProfiles(twoprofiles.BP, percentage = T, legend = T, aTitle = "comparison")

#=============================================================================================
#                                            {topGO}
#=============================================================================================

#retuns all GO terms associated with all AGIs in feasibleGenes list, see .pdf for more options...
# annFUN.org(whichOnto = "BP", feasibleGenes = genestest, mapping = "org.At.tair.db", ID = "entrez")
# 
# xx <- annFUN.org("BP", mapping = "org.At.tair.db", ID = "symbol")
# head(xx)
# allGenes <- unique(unlist(xx)) # allGenes is what i'm comparing against
#  myInterestedGenes <- sample(allGenes, 500) 
# mygenes <- rep(1,length(myInterestedGenes)) #make a vector of ones same length as genes of interest
# names(mygenes) <- myInterestedGenes #name this vector of ones with yer gene names
# alldagenes <- rep(1, length(allGenes))#do the same to make a named vector of allgenes
# names(alldagenes) <- allGenes
# GOdata <- new("topGOdata",
#                 ontology = "BP",
#                 allGenes = alldagenes,
#                 geneSel = mygenes,
#                 nodeSize = 5,
#                 annot = annFUN.org,
#                 mapping = "org.At.tair.db",
#                 ID = "symbol")
# 
# GenTable(object, ...)
# showGroupDensity(object, whichGO, ranks = FALSE, rm.one = TRUE)
# printGenes(object, whichTerms, file, ...)

# ------------------------------CREATE TOPGODATA OBJECT ----------------------------------------
xx <- annFUN.org("BP", mapping = "org.At.tair.db", ID = "entrez")
allGenes <- unique(unlist(xx)) 
allGenes <- allGenes[1:1000] #allGenes is yer list of AGIs

geneID2GO <- GOTermsList(allGenes, orgPkg= "org.At.tair.db", na.rm=TRUE) #get a named list of AGI where list names are AGI and list members are GOTerms
geneNames <- names(geneID2GO)
head(geneNames)

myInterestingGenes <- sample(geneNames, length(geneNames) / 10) #take random sample 1/10 length of big list
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
str(geneList)

GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOdata

description(GOdata)
graph(GOdata)

weight01.fisher <- runTest(GOdata, statistic = "fisher")
