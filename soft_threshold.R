#Choosing a soft threshold for WGCNA
#Danielle Jordan
#3.11.23

#load libraries
library(WGCNA)
library(flashClust)

#Read in dataset
log_bros_num <- read.csv("/uoa/home/r01dj21/sharedscratch/longread_analyses/WGCNA2/log_bros.csv", check.names = FALSE, stringsAsFactors = FALSE)

#Check the dimensions of the loaded dataset
cat("Dimensions of log_bros_num:", dim(log_bros_num), "\n")

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(log_bros_num, powerVector = powers, verbose = 5)

#Check the results of pickSoftThreshold
print(sft)

# Plot the results:
pdf(file = "softthresholdpower.pdf")
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.60,col="red", lty = 2)
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
