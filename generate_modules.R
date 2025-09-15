#Generate modules for WGCNA
#Danielle Jordan
#3.11.23

#load libraries
library(WGCNA)
library(flashClust)
allowWGCNAThreads()

#Read in dataset
log_bros <- read.csv("/uoa/home/r01dj21/sharedscratch/longread_analyses/WGCNA2/log_bros.csv", check.names = FALSE, stringsAsFactors = FALSE)

#Remove columns to ensure dataset is numeric
log_bros_num <- log_bros[ , unlist(lapply(log_bros, is.numeric))]

#build an adjacency "correlation" matrix
enableWGCNAThreads()
softPower = 10  #this will change depending on soft_threshold.R script
adjacency = adjacency(log_bros_num, power = softPower, type = "signed") #specify network type
head(adjacency)

# Construct Networks- USE A SUPERCOMPUTER IRL -----------------------------
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type #this might change
dissTOM = 1-TOM

# Generate a clustered gene tree
geneTree = flashClust(as.dist(dissTOM), method="average")
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
#This sets the minimum number of genes to cluster into a module
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)
dynamicColors= labels2colors(dynamicMods)
#Calculate eigengenes
MEList= moduleEigengenes(log_bros_num, colors= dynamicColors,softPower = softPower)
MEs= MEList$eigengenes
#Calculate dissimilarity of module eigengenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
# Create hierarchical clustering result
hc_result = flashClust(as.dist(MEDiss), method = "average")

# Convert the result to a dendrogram
METree = as.dendrogram(hc_result)
save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_allSamples_signed_RLDfiltered.RData")

#plots tree showing how the eigengenes cluster together
#INCLUE THE NEXT LINE TO SAVE TO FILE
pdf(file="clusterwithoutmodulecolors.pdf")
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
#set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
MEDissThres = 0.25 #change this to alter how many modules there are
merge = mergeCloseModules(log_bros_num, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
#INCLUE THE NEXT LINE TO SAVE TO FILE
dev.off()

#plot dendrogram with module colors below it
#INCLUE THE NEXT LINE TO SAVE TO FILE
pdf(file="cluster.pdf")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
#INCLUE THE NEXT LINE TO SAVE TO FILE
dev.off()

save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_allSamples_signed_nomerge_RLDfiltered.RData")

#######################################
#Module-trait relationships script ####
#######################################
#Read in trait dataset
datTraits_all <- read.csv("/uoa/home/r01dj21/sharedscratch/longread_analyses/WGCNA2/trait_data_dpf.csv", check.names = FALSE, stringsAsFactors = FALSE)
#Remove columns to ensure dataset is numeric
datTraits <- datTraits_all[ , unlist(lapply(datTraits_all, is.numeric))]

save(log_bros_num, datTraits, file = "SamplesAndTraits.RData")
#load("SamplesAndTraits.RData")

##Correlate traits to gene expression modules
#Define number of genes and samples
nGenes = ncol(log_bros_num)
nSamples = nrow(log_bros_num)
#Recalculate MEs with color labels
MEs0 = moduleEigengenes(log_bros_num, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
par(mar= c(6, 8.5, 3, 3))

#display the corelation values with a heatmap plot
#INCLUE THE NEXT LINE TO SAVE TO FILE
pdf(file="heatmap.pdf")
labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(datTraits),
               yLabels= names(MEs),
               ySymbols= names(MEs),
               colorLabels= FALSE,
               colors= blueWhiteRed(50),
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 0.5,
               zlim= c(-1,1),
               main= paste("Module-trait relationships"))
#INCLUE THE NEXT LINE TO SAVE TO FILE
dev.off()
