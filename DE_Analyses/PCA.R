#srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=1 --time=12:00:00 --job-name=qlogin /bin/bash -l
#module load R/3.6.2-foss-2019b
#R

#########################
######## SETUP ##########
#########################

library(DESeq2)
library(tximport)

library(matrixStats)
library(ggfortify)
library(ggplot2)

#setwd("/scratch/eld72413/Salty_Nut/CultivatedOnly/DE_Analyses_Inbred")

library( "genefilter" ) #no package?
library("gplots")
library("RColorBrewer")
library("viridis")

#########################
######### DATA ##########
#########################

#load("ddsTxi_cult_red.RData")

load("dds_cult.RData")

# df "Sig_AllTrt" of genes significant in at least 1 treatment
load("SigGenes.RData") 
dim(Sig_AllTrt) #23789

# transform
vsd <- vst(dds_cult, blind = FALSE)

#########################
######## ORDER? #########
#########################

# all significant genes
SigGenes <- vsd[which(rownames(vsd) %in% Sig_AllTrt$Gene),] # 23,789

# save transformed significant genes
save(SigGenes, file = "SigGenes_transformed.RData") #put in results directory?
load("SigGenes_transformed.RData")

SigGenes_Mat_Transform <- assay(SigGenes)
save(SigGenes_Mat_Transform, file = "SigGenes_Mat_Transform.RData")


#########################
###### PCA - LOCAL ######
#########################

load("ResultsFiles/SigGenes_Mat_Transform.RData")

# load grouping variables:
design <- read.csv("DataFiles/StudyDesign_Inbred_noOut.csv", header=T)

### Use ggplot & all genes
Matrix_all <- as.data.frame(t(SigGenes_Mat_Transform))
pca_all <- prcomp(Matrix_all, center = TRUE, scale. = TRUE)
pca_all$sdev
pca_resultsAll_graph <- cbind(design[,c(2,8,14,17)], pca_all$x[,1:4])

p <- ggplot(data=pca_resultsAll_graph, aes(x=PC1, y=PC2))
p+geom_point(aes(color=Treatment, shape=Accession))

p2 <- ggplot(data=pca_resultsAll_graph, aes(x=PC3, y=PC4))
p2+geom_point(aes(color=Treatment, shape=Accession))

# ggfortify extension
Matrix_data <- cbind(design[,c(2,8,14,17)], Matrix_all)
pca_res <- prcomp(Matrix_data[,-c(1:4)], scale. = TRUE)
pca_res2 <- prcomp(Matrix_data[,-c(1:4)], center = TRUE) #Accessions much more clustered this way
pca_res3 <- prcomp(Matrix_data[,-c(1:4)], center = TRUE, scale. = TRUE)
summary(pca_res3) # PC1=19.04%, PC2=9.4%
autoplot(pca_res3, data=Matrix_data, 
         colour='Treatment', shape='Accession', frame = TRUE) +
  theme_minimal()
#,frame.type = 'norm')

ggsave("/Users/emilydittmar/Google Drive/Active Projects/Transcriptomics_Exp/Analyses/Figures/Inbred_PCA.png")

# customize color scheme
Matrix_data$Treatment <- factor(Matrix_data$Treatment, levels=c("Control", "LowNut", "HighSalt", "Combo"))
levels(Matrix_data$Treatment)

display.brewer.pal(n=8, name = "Spectral") # 1,4,
display.brewer.pal(n=8, name = "Dark2") 
brewer.pal(n=8, name = "Dark2")
viridis(8)
display.brewer.pal(n=8, name = "Set1") 
brewer.pal(n=8, name = "Set1")

# dark2 colors
autoplot(pca_res3, data=Matrix_data, 
         colour='Treatment', shape='Accession', frame = TRUE) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3")) +
  scale_colour_manual(values = c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3")) +
  theme_minimal()

# set1 colors with dark2 yellow for salt
autoplot(pca_res3, data=Matrix_data, 
         colour='Treatment', shape='Accession', frame = TRUE) +
  scale_fill_manual(values = c("#4DAF4A", "#E41A1C", "#E6AB02", "#377EB8")) +
  scale_colour_manual(values = c("#4DAF4A", "#E41A1C", "#E6AB02", "#377EB8")) +
  theme_minimal()
ggsave("/Users/emilydittmar/Google Drive/Active Projects/Transcriptomics_Exp/Manuscript/SaltxNut/Figures/PCA.png")

plotMDS(SigGenes_Mat_Transform, labels = design$Treatment) # need deseq2

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


#########################
####### HEATMAP #########
#########################

pdf(file = "Heatmap.pdf", width = 12, height = 12);
heatmap.2( assay(vsd)[ topVarGenes, ], scale="row",
trace="none", dendrogram="column",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()

pdf(file = "HeatmapSig.pdf", width = 12, height = 12);
heatmap.2( head(SigGenes,35), scale="row",
trace="none", dendrogram="column",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()

pdf(file = "HeatmapSig2.pdf", width = 12, height = 12);
heatmap.2( SigGenes[ topSigGenes, ], scale="row",
trace="none", dendrogram="column",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()

pdf(file = "HeatmapSigAll.pdf", width = 12, height = 12);
heatmap.2( SigGenes, scale="row",
trace="none", dendrogram="column",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()

###############################
##  PCA WITH DESEQ2 FUNCTION ##
###############################

pcaData <- plotPCA(vsd, intgroup = c("Accession", "Treatment"), ntop = 500, returnData = TRUE) 
pcaData <- plotPCA(SigGenes, intgroup = c("Accession", "Treatment"), ntop = 500, returnData = TRUE) 

p <- ggplot(pcaData, aes(PC1, PC2, colour=Accession))

png("PCA.png")
p + geom_point(aes(shape=Treatment))
dev.off()

png("PCA_wlabels.png")
p + geom_text(aes(label=name))
dev.off()

# outliers: 63, 335, 37, 261, 336

##### Using prcomp:

# load grouping variables:
design <- read.csv("DataFiles/StudyDesign_Inbred_noOut.csv", header=T)

topSigGenes <- head( order( rowVars( SigGenes_Mat_Transform ), decreasing=TRUE ), 500 )
top500Sig <- as.data.frame(SigGenes_Mat_Transform[topSigGenes,])

Matrix <- as.data.frame(t(top500Sig)) #transpose matrix
Matrix$Plant <- rownames(Matrix)

Matrix_data <- merge(design[,c(2,8,14,17)], Matrix, by="Plant")

pca_results <- prcomp(Matrix_data[,-c(1:4)], center = TRUE)
plot(pca_results$x[,1:2], col=Matrix_data$Treatment)

pca_results2 <- prcomp(scale(Matrix_data[,-c(1:4)]))
plot(pca_results2$x[,1:2], col=Matrix_data$Treatment) #exactly the same


pca.plot <- autoplot(pca_results)
DE.shared.df <- as.data.frame(pca_results$rotation)
plot(DE.shared.df$PC1 ~ DE.shared.df$PC2)

DE.shared.labels <- cbind (DE.shared.df, grouping)

#requires matrixStats
Make_PCA <- function(labeldf, vsddf, Num_genes, grouping) {
  topSigGenes <- head( order( rowVars( vsddf ), decreasing=TRUE ), Num_genes )
  topNSig <- as.data.frame(vsddf[topSigGenes,])
  transpose <- as.data.frame(t(topNSig))
  transpose$Plant <- rownames(transpose)
  transpose_data <- merge(labeldf, transpose, by="Plant")
  pca_results <- prcomp(transpose_data [,-c(1:4)], center = TRUE)
  plot(pca_results$x[,1:2], col=transpose_data[,grouping])
  #return (pca_results)
}

Make_PCA(design[,c(2,8,14,17)], SigGenes_Mat_Transform, 23789, "Treatment")
plot(test$x[,1:2]) 

biplot(pca_results, choices = 1:2, scale =1, pc.biplot = FALSE)