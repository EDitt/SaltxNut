#### WGCNA ANALYSES ####

#qsub -I -q s_interq -l walltime=8:00:00 -l nodes=1:ppn=8 -l mem=50gb
#module load R/3.6.1-foss-2018a-X11-20180131-GACRC
#R

setwd("/scratch/eld72413/Salty_Nut/CultivatedOnly/DE_Analyses_Inbred")

#########################
######## SETUP ##########
#########################

library(DESeq2)
library(tximport)
library(WGCNA)
#library(sva)

### To parallelize:
library("BiocParallel")
register(MulticoreParam(8)) 


#############################
##### LOAD DATA OBJECTS #####
#############################

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
keep <- rowSums(counts(dds_cult) >= 10) >= 80 ### at least 80 samples have a count of 10 or higher
length(which(keep==1)) #23236
dds_cult_red1 <- dds_cult[keep,]
dim(dds_cult_red1) #23236 x 104

vsd_cult <- vst(dds_cult_red1 , blind = FALSE)

head(assay(vsd_cult), 3)

### To extract matrix of normalized values
mat_norm <- assay(vsd_cult)
dim(mat_norm) #23236 x 104

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
Accession=metadata_cult$Accession
Accession.num = as.numeric(metadata_cult$Accession) #needs to be numeric for plotClusterTreeSamples function

#############################
###### QUALITY CONTROL ######
#############################

#determine mean expression per array
meanExpressionByArray=apply( datExpr,1,mean, na.rm=T)

pdf(file = "mean_expression_array_cultivated_red_no.outlier.pdf", width = 8, height = 8)
barplot(meanExpressionByArray,
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across samples (Cultivated)",
        names.arg = c(metadata_cult$Plant), cex.names = 0.7)
dev.off()

#remove genes with no variance
#first, calculate the variances of the probes (genes) and the number of present entries
variancedatExpr=as.vector(apply(as.matrix(datExpr),2,var, na.rm=T))
no.presentdatExpr=as.vector(apply(!is.na(as.matrix(datExpr)),2, sum) ) # number of entries=number of samples/genes

#keep only genes whose variance is non-zero and have at least 4 present entries
KeepGenes= variancedatExpr>0 & no.presentdatExpr>=4

table(KeepGenes) # true for all

#############################
##### CLUSTERING/PLOTS ######
#############################

#y should be a matrix in which each column corresponds to a different trait and each row to a sample
y = metadata_cult[,c(1,3:8,14,17,18,23,33,42,43)]
y$Treatment <- as.numeric(y$Treatment)
y$Accession <- as.numeric(y$Accession)
y$Group <- as.numeric(y$Group)
y$Reproductive <- as.numeric(y$Reproductive)
#(can automatically transform the other numbers to numeric)

# detection of outliers
y_model <- model.matrix(~ 0 + Group + Cross + Treatment + Reproductive, data=metadata_cult)
y_model_plot <- as.data.frame(y_model)

# Cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")

pdf(file = "cultivated_cluster.pdf", width = 10, height = 8)
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(y_model_plot, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(y_model_plot),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

#############################
###### REMOVE OUTLIERS ######
#############################

cutHeights = 160
# Re-plot the dendrograms including the cut lines
pdf(file = "SampleClustering_wline.pdf", width = 12, height = 12);
plot(sampleTree2, main = "Sample clustering on all genes", 
xlab="", sub="", cex = 0.7);
abline(h=cutHeights, col = "red");
dev.off()

# Find clusters cut by the line
clust = cutreeStatic(sampleTree2, cutHeight = 180)
table(clust)
keepSamples = (clust==1)
datExpr_red = datExpr[keepSamples, ] # removed 2 outliers
nGenes = ncol(datExpr_red)
nSamples = nrow(datExpr_red)

sampleNums = rownames(datExpr_red) #N=102

# outliers remove from metadata:
y_red <- y[which(y$Plant %in% sampleNums),]
rownames(y_red) <- y_red$Plant
y_red = y_red[, -1]

#check
length(match(rownames(y_red), rownames(datExpr_red))) # N=102

# Re-cluster samples
sampleTree3 = hclust(dist(datExpr_red), method = "average")
y_model2 <- model.matrix(~ 0 + Group + Cross + Treatment + Reproductive, data=y_red)
y_model_plot2 <- as.data.frame(y_model2)

pdf(file = "cultivated_cluster_no.outliers.pdf", width = 10, height = 8)
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(y_model_plot2, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree3, traitColors,
                    groupLabels = names(y_model_plot2),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

save(y_red, datExpr_red, file = "NetworkData_outliersRem.RData")


