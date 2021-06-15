
# On MSI:

#############################
###### SETUP WORKSPACE ######
#############################

# tmux new -s wgcna
# srun -N 1 -n 1 -c 1 --mem=50gb -t 12:00:00 -p interactive --pty bash
# module load R/3.6.0
# cd /scratch.global/edittmar/WGCNA
# R

setwd("/scratch.global/edittmar/WGCNA")
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

load("multiExpr.RData")
checkSets(multiExpr)
nSets = checkSets(multiExpr)$nSets

load("Consensus-dataInput.RData")
load("ConsensusTOM_signed.RData")

# To continue where I left off (in 2.a.6):

#############################
## HIERARCHICAL CLUSTERING ##
#############################

# Clustering
consTree = hclust(as.dist(1-consensusTOM), method = "average")
format(object.size(consTree), units = "auto") # "674.3 Kb"

# Plot the resulting clustering tree (dendrogram)
pdf(file = "Results/ConsensusTree_signed.pdf", width = 8, height = 8)
plot(consTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

save(consTree, file = "ConsensusTree_signed.RData")

load("ConsensusTree_signed.RData")

#############################
### MODULE ID EXPLORATORY ###
#############################

# explore deepSplit and pamStage parameters

#dissTOM = 1 - consensusTOM

split = list(0,1,2,3,4)

true_list <- lapply(split, function(x) {Module_ID_multiExpr(multiExpr, consTree, 1-consensusTOM, x, TRUE, 30) } )
#Error in moduleEigengenes(dataframe, colors = dynamicColors) : 
#  moduleEigengenes: Error: expr must be two-dimensional.

true_df <- do.call("rbind", true_list)
true_df$pam <- "TRUE"
true_df$split <- c(0,1,2,3,4)

false_list <- lapply(split, function(x) {Module_ID_multiExpr(multiExpr, consTree, 1-consensusTOM, x, FALSE, 30) } )
false_df <- do.call("rbind", false_list)
false_df$pam <- "FALSE"
false_df$split <- c(0,1,2,3,4)

df_module_results <- rbind(true_df, false_df)

write.csv(df_module_results, 'TreeCutVary.csv')

# highest mean variance explained (0.554) is for split=2, pam=FALSE, 17k genes un-assigned
# for pam=TRUE, split=3 looks the best (4 has highest mean var explained but very little difference)

#############################
### MODULE IDENTIFICATION ###
#############################

minModuleSize = 30
# Module identification using dynamic tree cut:
unmergedLabels = cutreeDynamic(dendro = consTree, distM = 1-consensusTOM,
     method = "hybrid",
     deepSplit = 2,
     pamRespectsDendro = FALSE,
     minClusterSize = minModuleSize,
     pamStage = FALSE,
     verbose = 4)

unmergedColors = labels2colors(unmergedLabels)

table(unmergedLabels)
# 93 modules

pdf(file = "Results/ConsensusTree_Split2_PAMfalse.pdf", width = 8, height = 6)
plotDendroAndColors(consTree, unmergedColors, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

#############################
#### MODULE SIMILARITY ######
#############################

# Calculate module eigengenes
unmergedMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedColors, verbose = 5)
# the default is align = "along average" which controls whether eigengenes are aligned with average expression

### Grey module (set 2 only) had errors:
#   moduleEigengenes : Working on ME for module grey
#    ... 17160 genes
#    ..principal component calculation for module grey failed with the following error:
#         Error in svd(datModule, nu = min(n, p, nPC), nv = min(n, p, nPC)) : 
#  infinite or missing values in 'x'
#     ..hub genes will be used instead of principal components.

# Calculate consensus dissimilarity of consensus module eigengenes
consMEDiss = consensusMEDissimilarity(unmergedMEs)
# Cluster consensus modules
consMETree = hclust(as.dist(consMEDiss), method = "average")
# Plot the result
pdf(file = "Results/ConsensusCluster_EigengenesSplit2_pamFalse.pdf", width = 9, height = 6)
par(mfrow = c(1,1))
plot(consMETree, main = "Consensus clustering of consensus module eigengenes",
xlab = "", sub = "")
abline(h=0.25, col = "red")
dev.off()

save(unmergedMEs, file = "unmergedMEs.RData")
load("unmergedMEs.RData")

#############################
# TEST CUT HEIGHT FOR MERGE #
#############################

# how does the tree cut height value influence the average number of genes in each module &
# the variance explained by eigengenes?

#  threshold between 0.1 - 1.0, by 0.1
Threshold = seq(from = 0.05, to = 1.0, by = 0.05)

df_module_results <- data.frame(Threshold, ModuleNumber = numeric(length(Threshold)), 
     MeanVarExp=numeric(length(Threshold)), MinVarExp=numeric(length(Threshold)), 
     MaxVarExp=numeric(length(Threshold)))

for (i in 1:length(Threshold)) {
  mod_result <- Module_Cluster_multiExpr(df_module_results[i,1], multiExpr, unmergedLabels)
  df_module_results$ModuleNumber[i] = mod_result[1,2]
  df_module_results$MeanVarExp[i] = mod_result[2,2]
  df_module_results$MinVarExp[i] = mod_result[3,2]
  df_module_results$MaxVarExp[i] = mod_result[4,2]
}

write.csv(df_module_results, 'VariableCutThreshold_PAMFALSE_split2.csv')

### cutting at 0.15 reduces % variance explained to 0.5355 and 89 modules (from 0.55 and 94 modules)

#############################
### MERGE SIMILAR MODULES ###
#############################

merge = mergeCloseModules(multiExpr, unmergedColors, cutHeight = 0.15, verbose = 3)
#
names(merge$newMEs[[1]])
# see 'moduleEigengenes' function for explanation
# data - not sure what this is? (for each individual in each module) <- is this the eigengene?
# average expression (for each individual in each module)
# variance explained - for each module
# nPC = 1 (not sure what that means), validMEs (all TRUE), validColors
# allOK, allPC, isPC, isHub (length = number of modules, all are false),
# validAEs, allAEOK 

# Numeric module labels
moduleLabels = merge$colors

# Convert labels to colors
#moduleColors = labels2colors(moduleLabels) # no longer need since I used 'unmergedColors' instead of 'unmergedLabels' in merge function

# Eigengenes of the new merged modules:
consMEs = merge$newMEs

pdf(file = "Results/ConsensusTree_MergedModsSplit2_Cut0.15.pdf", width = 9, height = 6)
plotDendroAndColors(consTree, cbind(unmergedColors, moduleColors),
c("Unmerged", "Merged"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

# save
save(merge, file = "Merge_Split2_pamFALSE_cut0.15.RData")
save(consMEs, moduleLabels, consTree, file = "Consensus-NetworkConstruction-man.RData")


# Calculate consensus dissimilarity of merged consensus module eigengenes
consMEDiss_merge = consensusMEDissimilarity(consMEs)
# Cluster merged consensus modules
consMETree_merge = hclust(as.dist(consMEDiss_merge), method = "average")

# plot new modules
pdf(file = "Results/MergedModsSplit2_Cut0.15.pdf", width = 9, height = 6)
par(mfrow = c(1,1))
plot(consMETree_merge, main = "Consensus clustering of consensus (merged) module eigengenes",
xlab = "", sub = "")
abline(h=0.15, col = "red")
dev.off()

#############################
#### OUTPUT MODULE INFO #####
#############################

### module eigengene info
Eigengenes_Group1 <- consMEs[[1]]$data
Eigengenes_Group1$Sample <- rownames(Eigengenes_Group1)
Eigengenes_Group2 <- consMEs[[2]]$data
Eigengenes_Group2$Sample <- rownames(Eigengenes_Group2)
Eigengenes_all <- rbind(Eigengenes_Group1, Eigengenes_Group2)

write.csv(Eigengenes_all, file = "ModEigengenes_Consensus.csv", row.names = FALSE, quote = FALSE)
# module eigengene is the first principal component of the expression matrix. The eigengene 
# can be thought of as a weighted average expression profile

### module expression info
Expression_Group1 <- consMEs[[1]]$averageExpr
Expression_Group1$Sample <- rownames(Expression_Group1)

Expression_Group2 <- consMEs[[2]]$averageExpr
Expression_Group2$Sample <- rownames(Expression_Group2)

Expression_all <- rbind(Expression_Group1, Expression_Group2)
write.csv(Expression_all, file = "ModExpression_Consensus.csv", row.names = FALSE, quote = FALSE)
# averageExpr is a df containing average normalized expression in each module

### variance explained by each module
VarExplained_1 <- consMEs[[1]]$varExplained
VarExplained_2 <- consMEs[[2]]$varExplained
VarExplained_1$Group <-c("HA")
VarExplained_2$Group <-c("RHA")

VarExplained_all <- rbind(VarExplained_1, VarExplained_2)
write.csv(VarExplained_all, file = "ModVarExplained.csv", row.names = FALSE, quote = FALSE)

####
# do I need this?
#MExp = orderMEs(MExp0) #not sure how these are ordered?

