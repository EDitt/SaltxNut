#qsub -I -q s_interq -l walltime=4:00:00 -l nodes=1:ppn=8 -l mem=50gb
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

### To parallelize:
library("BiocParallel")
register(MulticoreParam(8)) #register cores so can specify parallel=TRUE when need to parallelize. Need to increase memory if increasing cores.
#here I used 8 cores and 50 gb memory/core

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

outliers <- c(37, 45, 63, 186, 254, 261, 323, 335, 336)

experiment_inbred_noOut1 <- experiment_inbred[-which(experiment_inbred$Plant %in% outliers),] #104
rownames(experiment_inbred_noOut1) <- experiment_inbred_noOut1$Plant

metadata_red <- experiment_inbred_noOut1[,c(1,8,14,17)]
write.csv(metadata_red, "Metadata_forNiki.csv")

#############################
# IMPORT DATA WITH TXIMPORT #
#############################

# Directory containing the .genes.results files,
#setwd("/scratch/eld72413/Salty_Nut/CultivatedOnly/DE_Analyses/RSEMOut_Files")

# create vector pointing to quantification files - vector of filenames that contains sample IDs and combine with
files_cult <- file.path("/scratch/eld72413/Salty_Nut/CultivatedOnly/DE_Analyses/RSEMOut_Files", paste0("RSEMOut_", experiment_inbred_noOut1$Plant, ".merged_.genes.results"))
names(files_cult) <- experiment_inbred_noOut1$Plant #assign sampleIDs

all(file.exists(files_cult)) #make sure files actually exist (=TRUE)

txi.rsem.cult <- tximport(files_cult, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem.cult$counts)
dim(txi.rsem.cult$counts) #78260 x 104

#Save this R object
#save(txi.rsem.cult, file = "txi.rsem.cult.RData")


###############################
# CONSTRUCT DESEQ DATA OBJECT #
###############################

### Variables
experiment_inbred_noOut1$SampleDay<-factor(experiment_inbred_noOut1$SampleDay) #day sampled 1-2
experiment_inbred_noOut1$Cross<-factor(experiment_inbred_noOut1$Cross) #1,2,3

### Variables of interest:
levels(experiment_inbred_noOut1$Treatment)
experiment_inbred_noOut1$Treatment <- relevel(experiment_inbred_noOut1$Treatment, ref="Control") #relevel to make Control reference
experiment_inbred_noOut1$Accession<-as.factor(make.names(experiment_inbred_noOut1$Accession))  #make syntactically valid names (no characters R doesn't like)
levels(experiment_inbred_noOut1$Accession)
levels(experiment_inbred_noOut1$Group) #HA, RHA

### Other variables to account for:
levels(experiment_inbred_noOut1$Reproductive) #whether plant was reproductive when sampled


## make DESeq2 data object:
ddsTxi_cult <- DESeqDataSetFromTximport(txi.rsem.cult,
                                        colData = experiment_inbred_noOut1,
                                        design = ~ 0 + Group + Group:Cross + SampleDay + Reproductive + Treatment)
ddsTxi_cult #dim is 78260 x 104
### Pre-filtering to reduce size of data object:
keep <- rowSums(counts(ddsTxi_cult) >= 1) >= 3 ### at least 3 samples have a count of 1 or higher
length(which(keep==1)) #55,728
ddsTxi_cult_red <- ddsTxi_cult[keep,]
dim(ddsTxi_cult_red) #55728 x 104



### For Niki
head(assay(ddsTxi_cult_red), 3)


transform <- t(assay(ddsTxi_cult_red))
write.csv(transform, "GeneCounts_forNiki.csv")

## Variance stabilizing transformation
vsd <- vst(txi.rsem.cult, blind = FALSE) #return a DESeqTransform object. Transformed values are no longer counts. The colData attaached to dds is still accessible
# using blind=FALSE means that variables in design will not contribute to the expected variance-mean trend of the exp. 
# The design is not used directly in the transformation, only in estimating the global amount of variability in the counts
# Use blind = TRUE for a fully unsupervised transformation
# sequencing depth correction is done automatically with this function

head(assay(vsd), 3)

vsd <- vst(ddsTxi_cult, blind = TRUE)

png("PCA_wlabels_Unsupervised.png")
p + geom_text(aes(label=name))
dev.off()