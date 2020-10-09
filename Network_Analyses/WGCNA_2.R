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
source("Functions.R")

setwd("/scratch/eld72413/Salty_Nut/CultivatedOnly/DE_Analyses_Inbred")

### To parallelize:
library("BiocParallel")
register(MulticoreParam(8)) #register cores so can specify parallel=TRUE when need to parallelize. Need to increase memory if increasing cores.
#here I used 8 cores and 50 gb memory/core


#############################
##### LOAD DATA OBJECTS #####
#############################

load("NetworkData_outliersRem.RData")
dim(y_red) # 102 x 13
dim(datExpr_red) # 102  23236


###############################
# CHOOSE SOFT-THRESHOLD POWER #
###############################
# library(cluster) do I need this?

### to multi-thread
enableWGCNAThreads(nThreads = 8)

# Need to choose soft thresholding power, beta 
# Co-expression similarity is transformed into adjacency (raised to power of beta for weighted network)
# soft-threshold power is a trade-off between scale free topology and mean connectivity

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
# Using biweight mid-correlation instead of default Pearson
sft = pickSoftThreshold(datExpr_red, 
                        powerVector = powers, 
                        corFnc = bicor,
                        networkType = "signed",
                        verbose = 5)

# Plot the results:
pdf(file = "Soft_Threshold_bicor.pdf", width = 8, height = 8)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"),
     ylim=c(0,1.0));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

#### Choose value of 12 - lowest power for which scale-free topology fit index reaches 0.9
#### topology index doesn't quite reach 0.9 however

#### Does this value show a network that exhibits scale free topology?
k=softConnectivity(datE=datExpr_red,power=12)
# softConnectivity: FYI: connecitivty of genes with less than 34 valid samples will be returned as NA.

# Plot a histogram of k and a scale free topology plot
pdf(file = "k and scale free topology.pdf", width = 8, height = 8)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")
dev.off()


### Define beta and calculate adjacency using the soft thresholding power 12
softPower = 12
adjacency = adjacency(datExpr_red,
                      type = "signed",
                      power = softPower,
                      corFnc = bicor)


################################
## TOPOLOGICAL OVERLAP MATRIX ##
################################

# transform adjacency into Topological Overlap Matrix to minimize effects of noise and spurious associations
TOM = TOMsimilarity(adjacency)
### this step takes awhile

save(TOM, file = "TOM_Inbred.RData")

#load("TOM_Inbred.RData")

# calculate corresponding dissimilarity
dissTOM = 1-TOM


#############################
## HIERARCHICAL CLUSTERING ##
#############################

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")

# Plot the resulting clustering tree (dendrogram)
pdf(file = "GeneTree.pdf", width = 8, height = 8)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

#save(geneTree, file = "GeneTree.RData")

#############################
### MODULE ID EXPLORATORY ###
#############################

# explore deepSplit and pamStage parameters

split = list(0,1,2,3,4)

true_list <- lapply(split, function(x) {Module_ID(datExpr_red, geneTree, dissTOM, x, TRUE, 30) } )
true_df <- do.call("rbind", true_list)
true_df$pam <- "TRUE"
true_df$split <- c(0,1,2,3,4)

false_list <- lapply(split, function(x) {Module_ID(datExpr_red, geneTree, dissTOM, x, FALSE, 30) } )
false_df <- do.call("rbind", false_list)
false_df$pam <- "FALSE"
false_df$split <- c(0,1,2,3,4)

df_module_results <- rbind(true_df, false_df)

write.csv(df_module_results, 'TreeCutVary.csv')

# highest average % var explained is with deepsplit = 0, and pamStage = FALSE (61%)
# >16k unassigned genes

# highest average $ var explained when pamStage = TRUE is with deepsplit = 4 (49.66%)
# pamStage = TRUE, deepsplit = 4 (ave: 197 genes per module, 118 modules)
# pamStage = TRUE, deepsplit = 3 explains 49.28% (~ 240 genes/module, 97 modules)

### because very little difference in % var explained, will use deepsplit = 3 (pamStage=TRUE)


#############################
### MODULE IDENTIFICATION ###
#############################

# deepSplit = 3, pamStage=TRUE

minModuleSize = 30
# branch cutting methods include constant-height or Dynamic Branch cut

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 3, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize,
                            pamStage = TRUE)
table(dynamicMods) # N=97


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
pdf(file = "GeneTree_Modules_pamStageTRUE.pdf", width = 8, height = 8)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


#############################
#### MODULE SIMILARITY ######
#############################

# Calculate eigengenes
MEList = moduleEigengenes(datExpr_red, colors = dynamicColors)
MEs = MEList$eigengenes
MEList$varExplained

# Mean % Var Explained = 0.4928
mean(data.frame(matrix(unlist(MEList$varExplained), nrow=length(MEList$varExplained), byrow=T))[,1])

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

# Plot the result
pdf(file = "ModuleCluster_PAMstageTRUE_split3.pdf", width = 16, height = 8)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25  ### height cut of 0.25, corresponds to correlation of 0.75 (figure)
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
abline(h=0.2, col = "blue", lty=2)
abline(h=0.1, col = "green")
abline(h=0.05, col = "blue")
dev.off()


#############################
#### TEST MODULE CLUSTER ####
#############################

# how does the tree cut height value influence the average number of genes in each module &
# the variance explained by eigengenes?

#  threshold between 0.1 - 1.0, by 0.1
Threshold = seq(from = 0.05, to = 1.0, by = 0.05)

df_module_results <- data.frame(Threshold, ModuleNumber = numeric(length(Threshold)), 
                MeanNumGenes=numeric(length(Threshold)), MinNumGenes=numeric(length(Threshold)), 
                MaxNumGenes=numeric(length(Threshold)), MeanVarExp=numeric(length(Threshold)), 
                MinVarExp=numeric(length(Threshold)), MaxVarExp=numeric(length(Threshold)))

for (i in 1:length(Threshold)) {
  mod_result <- Module_Cluster(df_module_results[i,1], datExpr_red, dynamicColors)
  df_module_results$ModuleNumber[i] = mod_result[1,2]
  df_module_results$MeanNumGenes[i] = mod_result[2,2]
  df_module_results$MinNumGenes[i] = mod_result[3,2]
  df_module_results$MaxNumGenes[i] = mod_result[4,2]
  df_module_results$MeanVarExp[i] = mod_result[5,2]
  df_module_results$MinVarExp[i] = mod_result[6,2]
  df_module_results$MaxVarExp[i] = mod_result[7,2]
}

write.csv(df_module_results, 'VariableCutThreshold_PAMTRUE.csv')

### cutting at 0.20 maximizes the minimum % variance explained by a module
## average % var explained goes from 0.492 -> 0.482. Ave # genes ~ 327

#############################
### MERGE SIMILAR MODULES ###
#############################

#try new threshold
MEDissThres = 0.2
# Call an automatic merging function
merge = mergeCloseModules(datExpr_red, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

### See what merging did to module colors
pdf(file = "ModuleCluster_merged_0.2Cut.pdf", width = 8, height = 8)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

table(dynamicColors) #N=97
table(mergedColors) # N=71 (at 0.2)


### To use merged module colors, save relevant variables

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
dim(MEs) # 102 x 71

# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "ModuleFiles_100820.RData")

#############################
#### OUTPUT MODULE INFO #####
#############################

### module expression and eigengene info
write.csv(MEs, file = "ModEigengenes_pamTRUE.csv")

MExp0 = moduleEigengenes(datExpr_red, moduleColors)$averageExpr

MExp = orderMEs(MExp0) #not sure how these are ordered?
#write.csv(MExp, file = "ModExpression_pamFALSE.csv")
write.csv(MExp, file = "ModExpression_pamTRUE.csv")


#############################
# SUMMARIZE MODULE PROFILES #
#############################

### Module Eigengenes:
signif(cor(datME, use="p"), 2) #how correlated are modules

