#This script is for filtering and formatting the C.fornicata counts matrix 
#for WGCNA, which will be run on the UoA HPC
#Danielle Jordnan
#3.11.23

#Load libaries
library(edgeR)
library(tidyverse)
library(WGCNA)
library(flashClust)

#Read in data 
#NORMALIZED counts matrix 
cforn_data <- read.csv("cforn_hybrid_isoform_counts_matrix.csv", check.names = FALSE,
                       stringsAsFactors = FALSE)

#Make DGE list 
bros <- DGEList(counts = cforn_data[,2:27], genes = cforn_data[,1])

#Filtering the data to remove genes with very low counts (>0.1cpm in 6 or more libs)
#Want between 20k-40k genes
keep <- rowSums(cpm(bros)>0.1) >=6
bros <- bros[keep, , keep.lib.sizes=FALSE]

#How many genes do I have after filtering?
summary(bros$genes)
#74867, too many 

#Filter again (>5cpm in 3 or more libraries)
keep <- rowSums(cpm(bros)>3) >=3
bros <- bros[keep, , keep.lib.sizes=FALSE]

#How many genes do I have after filtering?
summary(bros$genes)
#29032--this is a good number

#Log transformation 
bros$logCPM <- log2(bros$counts[,-1] +1)

#Write out as csv and load back in to change from DGEList to csv
write.csv(bros, file = "bros_filtered.csv", row.names = FALSE, col.names = TRUE)

log_bros <- read.csv("bros_filtered.csv", check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)

