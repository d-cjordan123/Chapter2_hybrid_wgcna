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
cforn_data <- read.csv("cforn_hybrid_isoform_counts_matrix_raw.csv", check.names = FALSE,
                       stringsAsFactors = FALSE)

#Make DGE list 
bros <- DGEList(counts = cforn_data[,2:27], genes = cforn_data[,1])

#Filtering the data to remove genes with very low counts (>0.1cpm in 6 or more libs)
#Want between 20k-40k genes
keep <- rowSums(cpm(bros)>0.1) >=6
bros <- bros[keep, , keep.lib.sizes=FALSE]

#How many genes do I have after filtering?
summary(bros$genes)
#75872, too many 

#Filter again (>5cpm in 3 or more libraries)
keep <- rowSums(cpm(bros)>5) >=3
bros <- bros[keep, , keep.lib.sizes=FALSE]

#How many genes do I have after filtering?
summary(bros$genes)
#26724--this is a good number

#Log transformation 
bros$logCPM <- log2(bros$counts[,-1] +1)

#Write out as csv and load back in to change from DGEList to csv
write.csv(bros, file = "bros_filtered.csv", row.names = FALSE, col.names = TRUE)

log_bros <- read.csv("bros_filtered.csv", check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)

#Manipulate file to match WGCNA requested format 
#Basically need to flip the columns and rows 
row.names(log_bros) = log_bros$X
log_bros$X = NULL
log_bros = as.data.frame(t(log_bros))
dim(log_bros)

#Columns and rows are now flipped, but column names have been replaced. Need to set them back to Trinity_ID
#This will make the first row (the trinity IDs) the column names 
names(log_bros) <- log_bros[1,]
log_bros <- log_bros[-1,]

#For some reason the dataset is not read as a numeric file, and the easiest way to fix this 
#is to export it as a csv and then read it back in and it will read the data as a numeric dataset
#Also need to specify that the row names are indeed row names otherwise that messes up the numeric read as well
log_bros$row.names <- rownames(log_bros)

write.csv(log_bros, file = "log_bros.csv", row.names = FALSE, col.names = TRUE)

log_bros_num <- read.csv("log_bros.csv", row.names = "row.names", stringsAsFactors = FALSE, header = TRUE)

#This should now be a numeric dataset
#Now moving to see if there are any gene outliers 

gsg = goodSamplesGenes(log_bros_num, verbose = 3)

gsg$allOK

#True, so there are no gene outliers that need to be filtered 

###Now need to move to supercomputer for the remaining steps###