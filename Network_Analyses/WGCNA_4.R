
# To obtain information about each gene

load("Consensus-NetworkConstruction-man.RData")
# objects: consMEs, moduleLabels, consTree

load("multiExpr.RData")
nSets = checkSets(multiExpr)$nSets

#############################
# GET GENE MODULE ASSIGNMENTS 
#############################

length(colnames(multiExpr[[1]]$data)) #34,408
length(moduleLabels) #34,408
modGenes <- as.data.frame(cbind(moduleLabels, colnames(multiExpr[[1]]$data)))
colnames(modGenes)[2] <- c("Gene")
#check
aggregate(modGenes$Gene, by=list(modGenes$moduleLabels), length)
table(moduleLabels)

################################
###### MODULE MEMBERSHIP #######
################################
### For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile.
# this allows us to quantify the similarity of all genes to every module

# gene names
probes <- colnames(multiExpr[[1]]$data)

# recalculate module eigengenes in the "alphabetic" order and calculate gene significances and module memberships in each dataset
# I don't understand why ^ they do this, so I will take original module assignments

kME = list();
for (set in 1:nSets)
{
kME[[set]] = corAndPvalue(multiExpr[[set]]$data, consMEs[[set]]$data);
}

# combine the Z scores of correlations from each set to form a "meta-Z" score and corresponding p-value
kME.metaZ = (kME[[1]]$Z + kME[[2]]$Z)/sqrt(2)
kME.metaP = 2*pnorm(abs(kME.metaZ), lower.tail = FALSE)

# form matrices to hold the kME. Re-shaping trick to put the values and associated p-values and meta-analysis results next to each other
kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[1]]$p, kME[[2]]$p, kME.metaZ, kME.metaP);
MEnames = colnames(consMEs[[1]]$data);
nMEs = checkSets(consMEs)$nGenes
dim(kMEmat) = c(nGenes, 6*nMEs)
rownames(kMEmat) = probes;
colnames(kMEmat) = spaste(
c("kME.set1.", "kME.set2.", "p.kME.set1.", "p.kME.set2.", "Z.kME.meta.", "p.kME.meta"),
rep(MEnames, rep(6, nMEs)))

# put together full information and write into a plain text CSV file (probes not sorted in any particular way)

info = data.frame(Probe = probes,
ModuleLabel = moduleLabels,
kMEmat)

head(info[,1:3])
#check
aggregate(info$Probe, by=list(info$ModuleLabel), length)
table(moduleLabels)

write.csv(info, file = "Gene_ModuleMembership.csv", row.names = FALSE, quote = FALSE)

################################
## MODULE MEMBERSHIP CONSENSUS #
################################

# discovered functions: hierarchicalConsensusKME and consensusKME for calculating consensus module membership

# tmux new -s wgcna #lm0006
# srun -N 1 -n 1 -c 1 --mem=50gb -t 12:00:00 -p interactive --pty bash
# module load R/3.6.0
# cd /scratch.old/edittmar/WGCNA # msi changed global scratch path
# R

setwd("/scratch.old/edittmar/WGCNA")

library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

load("Consensus-NetworkConstruction-man.RData")
# objects: consMEs, moduleLabels, consTree
length(moduleLabels)

load("multiExpr.RData")
nSets = checkSets(multiExpr)$nSets

setLabels = c("HA", "RHA")

HierConKME <- consensusKME(
multiExpr,
moduleLabels,
multiEigengenes = NULL,
consensusQuantile = 0,
signed = TRUE,
useModules = NULL,
metaAnalysisWeights = NULL,
corAndPvalueFnc = corAndPvalue, corOptions = list(), corComponent = "cor",
getQvalues = FALSE,
useRankPvalue = TRUE,
#rankPvalueOptions = list(calculateQvalue = getQvalues, pValueMethod = "scale"),
setNames = setLabels,
excludeGrey = TRUE,
greyLabel = if (is.numeric(moduleLabels)) 0 else "grey")

dim(HierConKME)
head(colnames(HierConKME))
ConskME <- HierConKME # not a hierarchical consensus
save(ConskME, file = "/panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Consensus_kME.RData")

# save df of relevant info:
colnames(ConskME)[1:20]
Modules <- unique(moduleLabels)

### kME values
consensusCols <- c(paste0("consensus.kME", Modules))
Mod_kMEs <- ConskME [names(ConskME) %in% consensusCols]
Mod_kMEs$ID <- row.names(Mod_kMEs)
write.csv(Mod_kMEs, "/panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Module_kMEs.csv",
	row.names=FALSE)

### p values
pvalCols <- c(paste0("meta.p.equalWeights.kME", Modules))
Mod_pVals <- ConskME [names(ConskME) %in% pvalCols]
Mod_pVals$ID <- row.names(Mod_pVals)
write.csv(Mod_pVals, "/panfs/roc/groups/9/morrellp/shared/Projects/Sunflower/Module_kME_Pvals.csv",
	row.names=FALSE)


################ couldn't get hierarchical consensus kME functions to work

#consTree$labels = c(colnames(multiExpr[[1]]$data))

consOptions <- newConsensusOptions(
calibration = "single quantile",
# Simple quantile scaling options
calibrationQuantile = 0.95,
sampleForCalibration = TRUE,
sampleForCalibrationFactor = 1000,
# Consensus definition
consensusQuantile = 0,
useMean = FALSE,
setWeights = NULL,
suppressNegativeResults = FALSE,
# Name to prevent files clashes
analysisName = "new")

consTree_new <- newConsensusTree(
	consTree)

######
HierConKME <- hierarchicalConsensusKME(
multiExpr,
moduleLabels,
multiWeights = NULL,
multiEigengenes = NULL,
consTree,
signed = TRUE,
useModules = NULL,
metaAnalysisWeights = NULL,
corAndPvalueFnc = corAndPvalue, corOptions = list(),
corComponent = "cor", getFDR = FALSE,
useRankPvalue = TRUE,
rankPvalueOptions = list(calculateQvalue = getFDR, pValueMethod = "scale"),
setNames = setLabels, excludeGrey = TRUE,
greyLabel = if (is.numeric(moduleLabels)) 0 else "grey",
reportWeightType = NULL,
getOwnModuleZ = TRUE,
getBestModuleZ = TRUE,
getOwnConsensusKME = TRUE,
getBestConsensusKME = TRUE,
getAverageKME = TRUE,
getConsensusKME = TRUE,
getMetaColsFor1Set = FALSE,
getMetaP = FALSE,
getMetaFDR = getMetaP && getFDR,
getSetKME = TRUE,
getSetZ = FALSE,
getSetP = FALSE,
getSetFDR = getSetP && getFDR,
includeID = TRUE,
additionalGeneInfo = NULL,
includeWeightTypeInColnames = TRUE)

#Error in names(consensusTree$inputs) <- spaste("Level.", level, ".Input.",  : 
#  attempt to set an attribute on NULL




#############################
##### GENE CONNECTIVITY #####
#############################

#The function intramodularConnectivity computes:

# ktotal = whole network connectivity
# kWithin = within module connectivity
# kOut = kTotal - kWithin
# kDiff = kIn - kOut = 2*kIn-kTotal

Alldegrees1 = list()
for (set in 1:nSets)
{
Alldegrees1[[set]] = intramodularConnectivity.fromExpr(multiExpr[[set]]$data, moduleLabels,
	networkType = "signed", ignoreColors = "grey",
	getWholeNetworkConnectivity = TRUE)
}


# softConnectivity: FYI: connecitivty of genes with less than 17 valid samples will be returned as NA.
# ..calculating connectivities....100% 
# softConnectivity: FYI: connecitivty of genes with less than 18 valid samples will be returned as NA.
# ..calculating connectivities....100% 

# how to summarize across the 2 datasets?

colnames(Alldegrees1[[1]]) <- paste0("HA_", colnames(Alldegrees1[[1]]))
colnames(Alldegrees1[[2]]) <- paste0("RHA_", colnames(Alldegrees1[[2]]))
All_Connectivity <- cbind(probes, Alldegrees1[[1]], Alldegrees1[[2]])

write.csv(All_Connectivity, 'Gene_Connectivity.csv', row.names = FALSE, quote = FALSE)

#############################
####### TOP HUB GENES ####### 
#############################

softPower = 12

top_hubs = list()
for (set in 1:nSets)
{
top_hubs[[set]] <- chooseTopHubInEachModule(
multiExpr[[set]]$data,
moduleLabels,
omitColors = "grey",
power = softPower,
type = "signed",
corFnc = bicor)
}

names(top_hubs[[1]])
top_hub_1 <- data.frame(matrix(unlist(top_hubs[[1]]), nrow=length(top_hubs[[1]]), byrow=T), stringsAsFactors=FALSE)
colnames(top_hub_1) <- "TopHub_HA"
top_hub_1$Modules <- names(top_hubs[[1]])

top_hub_2 <- data.frame(matrix(unlist(top_hubs[[2]]), nrow=length(top_hubs[[2]]), byrow=T), stringsAsFactors=FALSE)
colnames(top_hub_2) <- "TopHub_RHA"
top_hub_2$Modules <- names(top_hubs[[2]])

top_hubs_all <- merge(top_hub_1, top_hub_2, by="Modules")

write.csv(top_hubs_all , file = "topHubGenes.csv", row.names = FALSE, quote = FALSE)

