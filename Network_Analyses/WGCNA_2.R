#### WGCNA ANALYSES ####

#srun --pty  -p inter_p  --mem=100G --nodes=1 --ntasks-per-node=8 --time=6:00:00 --job-name=qlogin /bin/bash -l
#module load R/4.0.0-foss-2019b
#R

#########################
######## SETUP ##########
#########################

library(WGCNA)
source("Functions.R")

setwd("/scratch/eld72413/Salty_Nut/CultivatedOnly/DE_Analyses_Inbred")

### To parallelize:
library("BiocParallel")
register(MulticoreParam(8)) #register cores so can specify parallel=TRUE when need to parallelize. Need to increase memory if increasing cores.
#here I used 8 cores and 50 gb memory/core - > later changed to 100 gb memory

library(pryr) # to check how much memory is being used

#############################
##### LOAD DATA OBJECTS #####
#############################

load("multiExpr.RData")
checkSets(multiExpr)

nSets = checkSets(multiExpr)$nSets

load("Consensus-dataInput.RData")
nGenes
nSamples
setLabels
shortLabels
exprSize

head(Traits)


###############################
# CHOOSE SOFT-THRESHOLD POWER #
###############################

# library(cluster) do I need this?

### to multi-thread
enableWGCNAThreads(nThreads = 8)

# Need to choose soft thresholding power, beta 
# Co-expression similarity is transformed into adjacency (raised to power of beta for weighted network)
# soft-threshold power is a trade-off between scale free topology and mean connectivity

powers = c(seq(4,10,by=1), seq(12,20, by=2))

# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets)
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, networkType = "signed", verbose = 2)[[2]]) # didn't have networkType = "signed" at first
collectGarbage()

# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
"Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
for (col in 1:length(plotCols))
{
ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
}
}

# Plot the quantities in the chosen columns vs. the soft thresholding power
pdf(file = "Consensus_Network/Network_thresholds_signed.pdf", width = 8, height = 8)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
if (set==1)
{
plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
main = colNames[col]);
addGrid();
}
if (col==1)
{
text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
labels=powers,cex=cex1,col=colors[set]);
} else
text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
labels=powers,cex=cex1,col=colors[set]);
if (col==1)
{
legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
} else
legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off()

################################
#### CALCULATE ADJACENCIES #####
################################

### Define beta and calculate adjacency using the soft thresholding power 12
softPower = 5

# Initialize an appropriate array to hold the adjacencies
adjacencies = array(0, dim = c(nSets, nGenes, nGenes))

# Calculate adjacencies in each individual data set
#for (set in 1:nSets)
#adjacencies[set, , ] = adjacency(multiExpr[[set]]$data,
#                      type = "signed",
#                      power = softPower,
#                      corFnc = bicor)

# alternative? Does this do the same thing? It is more similar to code in manual
# In adjacency function documentation:
# for type = "unsigned", adjacency = |cor| ^ power
# for type = "signed", adjacency = (0.5 * (1 + cor)) ^ power

for (set in 1:nSets)
adjacencies[set, , ] = (0.5 * (1 + bicor(multiExpr[[set]]$data, 
                                             maxPOutliers = 0.1,
                                             pearsonFallback = "individual",
                                             nThreads = 8))) ^softPower

# code from manual to calculate adjacencies:
#for (set in 1:nSets)
#adjacencies[set, , ] = abs(cor(multiExpr[[set]]$data, use = "p"))^softPower

# Recommendations from FAQs:
# "signed is preferred". Use `type = "signed"`
# default correlation is Pearson. In general, they recommend the biweight mid-correlation

#  Warning messages:
#  1: In (function (x, y = NULL, robustX = TRUE, robustY = TRUE, use = "all.obs",  :
#    bicor: zero MAD in variable 'x'. Pearson correlation was used for individual columns with zero (or missing) MAD.
#  2: In (function (x, y = NULL, robustX = TRUE, robustY = TRUE, use = "all.obs",  :
#    bicor: zero MAD in variable 'x'. Pearson correlation was used for individual columns with zero (or missing) MAD.

# in bicor documentation:
# pearsonFallback - specifies whether the bicor calculation should revert to Pearson when median absolute deviation (mad) is zero. (Can be "none", "individual", "all"). 
# If set to none, zero mad will result in NA for the corresponding correlation. If set to "individual", Pearson calculation will be used only for columns that have zero mad.

# save(adjacencies, file="Consensus_Network/adjacencies.RData") # probably not worth saving because it takes longer to save than to calculate

mem_used() # 19.2 GB

################################
## TOPOLOGICAL OVERLAP MATRIX ##
################################

# transform adjacency into Topological Overlap Matrix to minimize effects of noise and spurious associations

# Initialize an appropriate array to hold the TOMs
TOM = array(0, dim = c(nSets, nGenes, nGenes))

mem_used() # 38.2 GB

# Calculate TOMs in each individual data set
for (set in 1:nSets)
TOM[set, , ] = TOMsimilarity(adjacencies[set, , ])
### this step takes awhile
# crashed with only 50 gb memory
# took ~ 90 min

object.size(TOM) # 18942567648 bytes
format(object.size(TOM), units = "auto") # "17.6 Gb"

save(TOM, file = "Consensus_Network/TOM_Inbred.RData")

#load("Consensus_Network/TOM_Inbred.RData")

ls()
format(object.size(adjacencies), units = "auto") # "17.6 Gb"

################################
########## SCALE TOMs ##########
################################

# scale TOMs to make them comparable across sets
# (different datasets have different statistical properties)

# Define the reference percentile
scaleP = 0.95

# Set RNG seed for reproducibility of sampling
set.seed(12345)

# Sample sufficiently large number of TOM entries
nSamples = as.integer(1/(1-scaleP) * 1000)

# Choose the sampled TOM entries
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
TOMScalingSamples = list()

# These are TOM values at reference percentile
scaleQuant = rep(1, nSets)

# Scaling powers to equalize reference TOM values
scalePowers = rep(1, nSets)

# Loop over sets
for (set in 1:nSets)
{
# Select the sampled TOM entries
TOMScalingSamples[[set]] = as.dist(TOM[set, , ])[scaleSample]
# Calculate the 95th percentile
scaleQuant[set] = quantile(TOMScalingSamples[[set]],
probs = scaleP, type = 8);
# Scale the male TOM
if (set>1)
{
scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
TOM[set, ,] = TOM[set, ,]^scalePowers[set];
}
}

# the array TOM now contains the scaled TOMs. Form a quantile-quantile plot of the TOMs before and after scaling:

# For plotting, also scale the sampled TOM entries
scaledTOMSamples = list();
for (set in 1:nSets)
scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set]

pdf(file = "Consensus_Network/TOMScaling-QQPlot.pdf", width = 6, height = 6);
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[2]], plot.it = TRUE, cex = 0.6,
xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[2]),
main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[2]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20);
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off()


#############################
## CALCULATE CONSENSUS TOM ##
#############################

consensusTOM = pmin(TOM[1, , ], TOM[2, , ])

format(object.size(consensusTOM), units = "auto") # "8.8 Gb"
rm(adjacencies)
mem_used() # 28.7 GB

save(consensusTOM, file = "Consensus_Network/ConsensusTOM.RData")
load("Consensus_Network/ConsensusTOM.RData")

#############################
## HIERARCHICAL CLUSTERING ##
#############################

# Clustering
consTree = hclust(as.dist(1-consensusTOM), method = "average")
format(object.size(consTree), units = "auto") # "674.3 Kb"

# Plot the resulting clustering tree (dendrogram)
pdf(file = "Consensus_Network/ConsensusTree.pdf", width = 8, height = 8)
plot(consTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

save(consTree, file = "Consensus_Network/ConsensusTree.RData")

load("Consensus_Network/ConsensusTree.RData")

##################### stopping point

#############################
### MODULE ID EXPLORATORY ###
#############################

# explore deepSplit and pamStage parameters

dissTOM = 1 - consensusTOM

split = list(0,1,2,3,4)

true_list <- lapply(split, function(x) {Module_ID(multiExpr, consTree, dissTOM, x, TRUE, 30) } )
#Error in moduleEigengenes(dataframe, colors = dynamicColors) : 
#  moduleEigengenes: Error: expr must be two-dimensional.


true_df <- do.call("rbind", true_list)
true_df$pam <- "TRUE"
true_df$split <- c(0,1,2,3,4)

false_list <- lapply(split, function(x) {Module_ID(multiExpr, consTree, dissTOM, x, FALSE, 30) } )
false_df <- do.call("rbind", false_list)
false_df$pam <- "FALSE"
false_df$split <- c(0,1,2,3,4)

df_module_results <- rbind(true_df, false_df)

write.csv(df_module_results, 'TreeCutVary.csv')


# previous results:

# highest average % var explained is with deepsplit = 0, and pamStage = FALSE (61%)
# >16k unassigned genes

# highest average $ var explained when pamStage = TRUE is with deepsplit = 4 (49.66%)
# pamStage = TRUE, deepsplit = 4 (ave: 197 genes per module, 118 modules)
# pamStage = TRUE, deepsplit = 3 explains 49.28% (~ 240 genes/module, 97 modules)

### because very little difference in % var explained, will use deepsplit = 3 (pamStage=TRUE)


#############################
### MODULE IDENTIFICATION ###
#############################

# deepSplit = 3 (4), pamStage=TRUE

minModuleSize = 30
# branch cutting methods include constant-height or Dynamic Branch cut

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 4, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize,
                            pamStage = TRUE)
table(dynamicMods) # 118
# (N=97 with deepSplit = 3)


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

# Mean % Var Explained = 0.4966 (0.4928 deepSplit3)
mean(data.frame(matrix(unlist(MEList$varExplained), nrow=length(MEList$varExplained), byrow=T))[,1])

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

# Plot the result
pdf(file = "ModuleCluster_PAMstageTRUE_split4.pdf", width = 16, height = 8)
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

write.csv(df_module_results, 'VariableCutThreshold_PAMTRUE_split4.csv')

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

table(dynamicColors) #N=118 (97 with deepSplit3)
table(mergedColors) # N=86 (71 with deepSplit3) (at 0.2)


### To use merged module colors, save relevant variables

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
dim(MEs) # 102 x 71/86

# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "ModuleFiles_100820.RData") #deepSplit3
save(MEs, moduleLabels, moduleColors, geneTree, file = "ModuleFiles_101020.RData") #deepSplit4

#############################
###### IF NOT MERGING #######
#############################

moduleColors = dynamicColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
dim(MEs) # 102 x 118  (already defined)

#############################
#### OUTPUT MODULE INFO #####
#############################

### module expression and eigengene info
write.csv(MEs, file = "ModEigengenes_pamTRUE.csv")
# module eigengene is the first principal component of the expression matrix. The eigengene 
# can be thought of as a weighted average expression profile

MExp0 = moduleEigengenes(datExpr_red, moduleColors)$averageExpr
# averageExpr is a df containing average normalized expression in each module

MExp = orderMEs(MExp0) #not sure how these are ordered?
#write.csv(MExp, file = "ModExpression_pamFALSE.csv")
write.csv(MExp, file = "ModExpression_pamTRUE.csv")


#############################
# SUMMARIZE MODULE PROFILES #
#############################

### Module Eigengenes:
signif(cor(datME, use="p"), 2) #how correlated are modules




###### Archive (before I made separate networks for HA & RHA)
#####
# Choose a set of soft-thresholding powers
#powers = c(c(1:10), seq(from = 12, to=20, by=2))
powers = c(seq(4,10,by=1), seq(12,20, by=2))

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


################################
#### CALCULATE ADJACENCIES #####
################################

# for one network:
adjacency = adjacency(datExpr_red,
                      type = "signed",
                      power = softPower,
                      corFnc = bicor)

TOM = TOMsimilarity(adjacency)

# calculate corresponding dissimilarity
dissTOM = 1-TOM