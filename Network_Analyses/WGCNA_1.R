#### WGCNA ANALYSES ####

#srun --pty  -p inter_p  --mem=50G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
#module load R/4.0.0-foss-2019b
#R

setwd("/scratch/eld72413/Salty_Nut/CultivatedOnly/DE_Analyses_Inbred")

# needed to install WGCNA- created a ~.Renviron file to specify location
install.packages('WGCNA')

### after errors in calculating eigengenes (even with tutorial), tried re-installing WGCNA with a diff version of R
# module load R/3.6.2-foss-2019b
# module load R/3.5.0-foss-2019b <- installed WGCNA in this version
# when I first tried to re-install WGCNA, I was warned of a non-zero exit status. I then used:
install.packages("BiocManager")
BiocManager::install("WGCNA")

# also realized that when I do library(WGCNA) it says the following objects are masked from 'package:stats': hclust, cor
# I think the first time it only listed cor here?
# This did not fix the problem

#########################
######## SETUP ##########
#########################

library(DESeq2)
library(tximport)
library(WGCNA)
source("Functions.R")

### To parallelize:
library("BiocParallel")
register(MulticoreParam(8)) 


#############################
##### LOAD DATA OBJECTS #####
#############################

options(stringsAsFactors = FALSE)

load("dds_cult.RData") #DESeq data object saved in DESeq.R script
dim(dds_cult) # 55728 x 104

# associated metadata
metadata_cult <- read.csv("StudyDesign_Inbred_noOut.csv", header=T)
head(metadata_cult)
dim(metadata_cult) # 104 x 48


#####################################
######## FILTER & TRANSFORM ########
#####################################

### Pre-filtering to reduce size of data object:
#keep <- rowSums(counts(dds_cult) >= 10) >= 80 ### at least 80 samples have a count of 10 or higher
keep <- rowSums(counts(dds_cult) >= 10) >= 10 ### at least 10 samples have a count of 10 or higher (based on WGCNA FAQ's)
length(which(keep==1)) #34748 (was 23236 when requiring 80 samples)
dds_cult_red1 <- dds_cult[keep,]
dim(dds_cult_red1) #23236 x 104

vsd_cult <- vst(dds_cult_red1 , blind = FALSE)

head(assay(vsd_cult), 3)

### To extract matrix of normalized values
mat_norm <- assay(vsd_cult)
dim(mat_norm) #34748 x 104

colnames(mat_norm) #sample names
head(rownames(mat_norm)) #gene names

options(stringsAsFactors = FALSE)

#make vector for sample names:
ArrayName = colnames(mat_norm) #N=104
#make vector for gene names:
GeneName = rownames(mat_norm) #N=23236

#transpose so that rows correspond to samples and columns correspond to genes
datExpr=data.frame(t(mat_norm)) #row and column names automatically transposed

#check
datExpr[1:5, 1:5]

#array names in the metadata file need to line up with those in the datExpr table
table( dimnames(datExpr)[[1]] == metadata_cult$Plant) #TRUE, 104

#define sample traits
Treatment=metadata_cult$Treatment
metadata_cult$Accession <- make.names(metadata_cult$Accession) # make syntactically valid names
Accession=metadata_cult$Accession
Accession.num = as.numeric(as.factor(metadata_cult$Accession)) #needs to be numeric for plotClusterTreeSamples function

#remove genes with no variance
#first, calculate the variances of the probes (genes) and the number of present entries
variancedatExpr=as.vector(apply(as.matrix(datExpr),2,var, na.rm=T))
no.presentdatExpr=as.vector(apply(!is.na(as.matrix(datExpr)),2, sum) ) # number of entries=number of samples/genes

#keep only genes whose variance is non-zero and have at least 4 present entries
KeepGenes= variancedatExpr>0 & no.presentdatExpr>=4

table(KeepGenes) # true for all

#############################
##### SEPARATE HA & RHA #####
#############################
metadata_cult$Group <- factor(metadata_cult$Group)
levels(metadata_cult$Group)

HA_metadata <- subset(metadata_cult, Group=="HA")
RHA_metadata <- subset(metadata_cult, Group=="RHA")

# to subset DESEQ OBJECT (ended up subsetting after transformation instead)
#dds_HA <- dds_cult[,which(colnames(dds_cult) %in% HA_metadata$Plant)]
#dim(dds_HA) # 55728 x 50
#dds_HA <- dds_cult[,dds_cult@colData$Plant %in% HA_metadata$Plant] # same thing
#dds_RHA <- dds_cult[,which(colnames(dds_cult) %in% RHA_metadata$Plant)]
#dim(dds_RHA) # 55728 x 54

datExpr_HA <- datExpr[which(rownames(datExpr) %in% HA_metadata$Plant),] 
dim(datExpr_HA) # 50 x 34748

datExpr_RHA <- datExpr[which(rownames(datExpr) %in% RHA_metadata$Plant),] 
dim(datExpr_RHA) # 54 x 34748

#############################
###### QUALITY CONTROL ######
#############################

#determine mean expression per array
meanAll <- mean(apply( datExpr,1,mean, na.rm=T))

# HA
meanExpressionByArrayHA=apply( datExpr_HA,1,mean, na.rm=T)
pdf(file = "Consensus_Network/mean_expression_array_HA.pdf", width = 8, height = 8)
barplot(meanExpressionByArrayHA,
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across HA samples",
        names.arg = c(HA_metadata$Plant), cex.names = 0.7,
        ylim=c(0,7))
abline(h=meanAll, lty=2)
dev.off()

# RHA
meanExpressionByArrayRHA=apply( datExpr_RHA,1,mean, na.rm=T)
pdf(file = "Consensus_Network/mean_expression_array_RHA.pdf", width = 8, height = 8)
barplot(meanExpressionByArrayRHA,
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across RHA samples",
        names.arg = c(RHA_metadata$Plant), cex.names = 0.7,
        ylim=c(0,7))
abline(h=meanAll, lty=2)
dev.off()


#############################
##### CLUSTERING/PLOTS ######
#############################

#y should be a matrix in which each column corresponds to a different trait and each row to a sample
y = metadata_cult[,c(2:8,14,17,23,33,42,43)]
y$Treatment <- as.numeric(y$Treatment)
y$Accession <- as.numeric(as.factor(y$Accession))
y$Group <- as.numeric(y$Group)
y$Reproductive <- as.numeric(y$Reproductive)
#(can automatically transform the other numbers to numeric)

# detection of outliers
y_model <- model.matrix(~ 0 + Group + Accession + Treatment + Reproductive, data=metadata_cult)
y_model_plot <- as.data.frame(y_model)

# Cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")

pdf(file = "Consensus_Network/cultivated_cluster.pdf", width = 10, height = 8)
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(y_model_plot, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(y_model_plot),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# sample 56 is an outlier


#############################
######## FORMAT DATA ########
#############################

nSets = 2
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("HA", "RHA")
shortLabels = c("HA", "RHA")

# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)

multiExpr[[1]] = list(data = datExpr_HA)
#names(multiExpr[[1]]$data) corresponds to gene names
#rownames(multiExpr[[1]]$data) corresponds to sample #'s'
multiExpr[[2]] = list(data = datExpr_RHA)
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)

# check
gsg = goodSamplesGenesMS(multiExpr, verbose = 3)
# output:
# Flagging genes and samples with too many missing values...
#  ..step 1
#  ..Excluding 340 genes from the calculation due to too many missing samples or zero variance.
#   ..bad gene count: 340, bad sample counts: 0, 0
#  ..step 2
#   ..bad gene count: 340, bad sample counts: 0, 0
gsg$allOK # FALSE


# remove offending samples & genes:
if (!gsg$allOK)
{
# Print information about the removed genes:
if (sum(!gsg$goodGenes) > 0)
printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
collapse = ", ")))
for (set in 1:exprSize$nSets)
{
2
if (sum(!gsg$goodSamples[[set]]))
printFlush(paste("In set", setLabels[set], "removing samples",
paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
# Remove the offending genes and samples
multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
}
# Update exprSize
exprSize = checkSets(multiExpr)
}

# now 34,408 genes (50 & 54 samples)

#############################
### CLUSTERING SEPARATELY ###
#############################

sampleTrees = list()
for (set in 1:nSets)
{
sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

pdf(file = "Consensus_Network/SampleClustering_groups.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
xlab="", sub="", cex = 0.7);
dev.off();

# one outlier in RHA dataset

# Choose the "base" cut height for the RHA data set
baseHeight = 200
cutHeights = c(200, 200)
# Re-plot the dendrograms including the cut lines
pdf(file = "Consensus_Network/SampleClustering_groups_wline.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
{
plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
xlab="", sub="", cex = 0.7);
abline(h=cutHeights[set], col = "red");
}
dev.off()

# remove outlier
for (set in 1:nSets)
{
# Find clusters cut by the line
labels = cutreeStatic(sampleTrees[[set]], cutHeight = cutHeights[set])
# Keep the largest one (labeled by the number 1)
keep = (labels==1)
multiExpr[[set]]$data = multiExpr[[set]]$data[keep, ]
}

collectGarbage()
# Check the size of the leftover data
exprSize = checkSets(multiExpr)
exprSize # n samples = 50 , 53 (1 removed from the RHA dataset)

# save
save(multiExpr, file = "multiExpr.RData")


#############################
######## TRAIT DATA #########
#############################

# traits of interest:
# 3=bench, 6=osmocote, 7=salt, 8=treatment, 14=accession, 17=group, 21=transplant date, 23=chlorophyll, 31= totAG, 32=totBG, 33=total, 42= reproductive, 43=sampleday, 
traits <- metadata_cult[,c(2,3,6:8,14,17,21,23,31:33,42,43)]
traits$Accession <- make.names(traits$Accession) # syntactically valid names

# Form a multi-set structure that will hold the clinical traits.
Traits = vector(mode="list", length = nSets)
for (set in 1:nSets)
{
setSamples = rownames(multiExpr[[set]]$data);
traitRows = match(setSamples, traits$Plant);
Traits[[set]] = list(data = traits[traitRows, -1]);
rownames(Traits[[set]]$data) = traits[traitRows, 1];
}

collectGarbage()

# Define data set dimensions
nGenes = exprSize$nGenes
nSamples = exprSize$nSamples

save(Traits, nGenes, nSamples, setLabels, shortLabels, exprSize, file = "Consensus-dataInput.RData")


################## Archived
# plot HA/RHA separately
y_groups <- split(y, as.factor(y$Group))
names(y_groups) <- c("HA", "RHA")
sampleTree_groups <- lapply(y_groups, function(x) hclust(dist(x), method = "average"))
y_model_plot_groups <- split(y_model_plot, as.factor(y_model_plot$GroupHA))
names(y_model_plot_groups) <- c("RHA", "HA")

#HA
pdf(file = "Consensus_Network/cultivated_clusterHA.pdf", width = 10, height = 8)
traitColors = numbers2colors(y_model_plot_groups$HA, signed = FALSE);
plotDendroAndColors(sampleTree_groups$HA, traitColors,
                        groupLabels = names(y_model_plot_groups$HA),
                        main = "Sample dendrogram HA")
dev.off()

#RHA
pdf(file = "Consensus_Network/cultivated_clusterRHA.pdf", width = 10, height = 8)
traitColors = numbers2colors(y_model_plot_groups$RHA, signed = FALSE);
plotDendroAndColors(sampleTree_groups$RHA, traitColors,
                        groupLabels = names(y_model_plot_groups$RHA),
                        main = "Sample dendrogram RHA")
dev.off()
