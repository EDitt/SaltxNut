
#module load R/3.6.1-foss-2018a-X11-20180131-GACRC
#R

#########################
######## SETUP ##########
#########################

library(DESeq2)
library(tximport)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

setwd("/scratch/eld72413/Salty_Nut/CultivatedOnly/DE_Analyses_Inbred")


#############################
### LOAD SAMPLE METADATA ####
#############################

### Experimental Data
#open up grouping file
experiment <- read.csv("StudyDesign.csv", header=T)

### Subset to include only cultivated HA + RHA samples inbred samples
experiment_inbred <- subset(experiment, Group == "HA" | Group == "RHA") #113
experiment_inbred <- droplevels(experiment_inbred)
dim(experiment_inbred)

#############################
# IMPORT DATA WITH TXIMPORT #
#############################

# create vector pointing to quantification files - vector of filenames that contains sample IDs and combine with
files_cult <- file.path("/scratch/eld72413/Salty_Nut/CultivatedOnly/DE_Analyses/RSEMOut_Files", paste0("RSEMOut_", experiment_inbred$Plant, ".merged_.genes.results"))
names(files_cult) <- experiment_inbred$Plant #assign sampleIDs

all(file.exists(files_cult)) #make sure files actually exist (=TRUE)

txi.rsem.cult <- tximport(files_cult, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem.cult$counts)
dim(txi.rsem.cult$counts) #78260 x 113

###############################
######## LIBRARY SIZES ########
###############################

## check library sizes

sums <- as.data.frame(colSums(txi.rsem.cult$counts), colnames(txi.rsem.cult$counts))
mean(sums[,1]) #13,003,905
sd(sums[,1]) #4,240,742
#standard error
sd(sums[,1])/sqrt(length(sums[,1])) #398,935.5
write.csv(sums, file.path("QC_Files/LibSums.csv"), row.names=TRUE)

png("QC_Files/LibSizeDist.png")
hist(sums)
dev.off()

low_lib <- rownames(sums)[which(sums < 5000000)] #45, 63, 254, 261, 336
MeanMinusSD <- mean(sums[,1]) - sd(sums[,1]) #8,763,163
rownames(sums)[which(sums < MeanMinusSD)] #45, 56, 63, 186, 244*, 254, 261, 298*, 336, 339*, 340*
rownames(sums)[which(sums < 8000000)] #45, 63, 186, 254, 261, 336

png("LibSizes.png")
barplot(sums[,1])
dev.off()

###############################
# CONSTRUCT DESEQ DATA OBJECT #
###############################

str(experiment_inbred)
### Variables
experiment_inbred$SampleDay<-factor(experiment_inbred$SampleDay) #day sampled 1-2
experiment_inbred$Cross<-factor(experiment_inbred$Cross) #1,2,3

### Variables of interest:
levels(experiment_inbred$Treatment)
experiment_inbred$Treatment <- relevel(experiment_inbred$Treatment, ref="Control") #relevel to make Control reference
experiment_inbred$Accession<-as.factor(make.names(experiment_inbred$Accession))  #make syntactically valid names (no characters R doesn't like)
levels(experiment_inbred$Accession)
levels(experiment_inbred$Group) #HA, RHA

### Other variables to account for:
levels(experiment_inbred$Reproductive) #whether plant was reproductive when sampled


## make DESeq2 data object:
ddsTxi_cult <- DESeqDataSetFromTximport(txi.rsem.cult,
                                        colData = experiment_inbred,
                                        design = ~ 0 + Group + Group:Cross + SampleDay + Reproductive + Treatment)
ddsTxi_cult #dim is 78260 x 113
### Pre-filtering to reduce size of data object:
keep <- rowSums(counts(ddsTxi_cult) >= 1) >= 3 ### at least 3 samples have a count of 1 or higher
length(which(keep==1)) #56,076
ddsTxi_cult_red <- ddsTxi_cult[keep,]
dim(ddsTxi_cult_red) #56076 x 113


###############################
###### TRANSFORMATION #########
###############################

## Variance stabilizing transformation
vsd <- vst(ddsTxi_cult_red, blind = FALSE) #return a DESeqTransform object. Transformed values are no longer counts. The colData attaached to dds is still accessible
# using blind=FALSE means that variables in design will not contribute to the expected variance-mean trend of the exp. 
# The design is not used directly in the transformation, only in estimating the global amount of variability in the counts
# Use blind = TRUE for a fully unsupervised transformation
# sequencing depth correction is done automatically with this function

head(assay(vsd), 3)

###############################
###########  PCA  #############
###############################

pcaData <- plotPCA(vsd, intgroup = c("Accession", "Treatment"), ntop = 500, returnData = TRUE) 

p <- ggplot(pcaData, aes(PC1, PC2, colour=Accession))

png("PCA.png")
p + geom_point(aes(shape=Treatment))
dev.off()

png("PCA_wlabels.png")
p + geom_text(aes(label=name))
dev.off()

# outliers: 63, 335, 37, 261, 336


###############################
###### SAMPLE DISTANCE ########
###############################

### Extract sample distances (using all genes)
sampleDists <- dist(t(assay(vsd)))

# Re-cluster samples
sampleTree = hclust(sampleDists, method = "average")

pdf(file = "SampleClustering_wline.pdf", width = 12, height = 12);
plot(sampleTree, main = "Sample clustering on all genes", 
     xlab="", sub="", cex = 0.7);
abline(h=180, col = "red");
dev.off()


###############################
###### REMOVE OUTLIERS ########
###############################

aggregate(experiment_inbred$Plant, by=list(experiment_inbred$Accession, experiment_inbred$Treatment), length)
# 3 reps: HA 445 Combo, HA 445 High Salt
# 4 reps: HA 412 HO Combo, RHA 373 Control, HA 445 Low Nutrient
# everything else 5

outliers <- c(37, 45, 63, 186, 254, 261, 323, 335, 336)
experiment_inbred_noOut1 <- experiment_inbred[-which(experiment_inbred$Plant %in% outliers),] 
aggregate(experiment_inbred_noOut1$Plant, by=list(experiment_inbred_noOut1$Accession, 
                                                  experiment_inbred_noOut1$Treatment), length)
# 3 reps: HA 445 combo, HA 445 High Salt

# previously removed:
#39 - RHA 373 Combo (currently 5 reps)
#67 - RHA 274 High Salt (currently 5 reps)

# also outlier
#56 - RHA 373 Salt (currently 4 reps)

head(experiment_inbred_noOut1)
dim(experiment_inbred_noOut1) #104 x 48
#write.csv(experiment_inbred_noOut1, file.path("DataFiles/StudyDesign_Inbred_noOut.csv"))
write.csv(experiment_inbred_noOut1, "StudyDesign_Inbred_noOut.csv", row.names=FALSE)
