# 

#############################
###### SETUP WORKSPACE ######
#############################

# tmux new -s wgcna
# srun -N 1 -n 1 -c 1 --mem=50gb -t 12:00:00 -p interactive --pty bash
# srun -N 1 -n 1 -c 1 --mem=22gb -t 12:00:00 -p interactive --pty bash
# module load R/3.6.0
# cd /scratch.global/edittmar/WGCNA
# R

setwd("/scratch.global/edittmar/WGCNA")
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# Load the data saved in the first part
lnames = load(file = "Consensus-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames # Traits, nGenes, nSamples, setLabels, shortLabels, exprSize
# Load the results of network analysis, tutorial part 2.a
lnames = load("Consensus-NetworkConstruction-man.RData")
lnames # consMEs, moduleLabels, consTree

load("multiExpr.RData")
exprSize = checkSets(multiExpr)
nSets = exprSize$nSets

#############################
### MODULE CORRESPONDANCE ###
#############################

# Create a variable for treatment 
Tot_Biomass = vector(mode = "list", length = nSets);
for (set in 1:nSets)
{
Tot_Biomass[[set]] = list(data = as.data.frame(Traits[[set]]$data$Tot_Biomass))
names(Tot_Biomass[[set]]$data) = "Tot_Biomass"
}

# We add the weight trait to the eigengenes and order them by consesus hierarchical clustering:
MET = consensusOrderMEs(addTraitToMEs(consMEs, Tot_Biomass))

pdf(file = "EigengeneNetworks.pdf", width= 8, height = 10);
par(cex = 0.9)
plotEigengeneNetworks(MET, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
dev.off()


#############################
## MODULE TRAIT CORRELATION ##
#############################

Model1 <- ~ 0 + Treatment + Accession
datTreatment <- lapply(Traits, function(x) {as.data.frame(model.matrix(Model1, data=x$data))})

# Set up variables to contain the module-trait correlations
moduleTraitCor = list()
moduleTraitPvalue = list()

# head(Traits[[1]]$data[,c(14,15)]) using only the salt and nutrient levels

# Calculate the correlations
for (set in 1:nSets)
{
moduleTraitCor[[set]] = cor(consMEs[[set]]$data, datTreatment[[set]], use = "p")
moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set]);
}

MEColors = names(consMEs[[1]]$data)
MEColorNames = names(consMEs[[1]]$data)


###########################
##### PLOT HA AND RHA #####
###########################

####
pdf(file = "ModuleTraitRelationships-HA.pdf", wi = 7, he = 10);
# Plot the module-trait relationship table for set number 1
set = 1
textMatrix = paste(signif(moduleTraitCor[[set]], 2), "\n(",
signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
xLabels = names(datTreatment[[set]]),
yLabels = MEColorNames,
ySymbols = MEColorNames,
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module--trait relationships in", setLabels[set]))
dev.off()

pdf(file = "ModuleTraitRelationships-RHA.pdf", wi = 7, he = 10);
# Plot the module-trait relationship table for set number 1
set = 2
textMatrix = paste(signif(moduleTraitCor[[set]], 2), "\n(",
signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
xLabels = names(datTreatment[[set]]),
yLabels = MEColorNames,
ySymbols = MEColorNames,
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module--trait relationships in", setLabels[set]))
dev.off()


#############################
###### CONSENSUS CORRS ######
#############################

# For each module-trait pair, take the correlation that has the lower absolute value in the two sets 
# if the two correlations have the same sign, and zero relationship if the two correlations have opposite signs

# Initialize matrices to hold the consensus correlation and p-value
consensusCor = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))
consensusPvalue = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))

# Find consensus negative correlations
negative = moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0
consensusCor[negative] = pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative])
consensusPvalue[negative] = pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative])

# Find consensus positive correlations
positive = moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0
consensusCor[positive] = pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive])
consensusPvalue[positive] = pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive])

textMatrix = paste(signif(consensusCor, 2), "\n(",
signif(consensusPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])

pdf(file = "ModuleTraitRelationships-consensus.pdf", wi = 7, he = 10)
labeledHeatmap(Matrix = consensusCor,
xLabels = names(datTreatment[[set]]),
yLabels = MEColorNames,
ySymbols = MEColorNames,
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Consensus module--trait relationships across\n",
paste(setLabels, collapse = " and ")))
dev.off()



#######################################
# SEPARATE MODEL FACTORS FOR SALT/NUT #
#######################################

# need to make salt and nut numeric
for (set in 1:nSets)
{
Traits[[set]]$data$LowNut <- ifelse(Traits[[set]]$data$Osmocote=="High", 0, 1)
Traits[[set]]$data$HighSalt <- ifelse(Traits[[set]]$data$Salt=="Salt", 1, 0)
}

Model2 <- ~ 0 + LowNut + HighSalt + LowNut:HighSalt + Accession
datTreatment2 <- lapply(Traits, function(x) {as.data.frame(model.matrix(Model2, data=x$data))})

moduleTraitCor2 = list()
moduleTraitPvalue2 = list()

# Calculate the correlations
for (set in 1:nSets)
{
moduleTraitCor2[[set]] = cor(consMEs[[set]]$data, datTreatment2[[set]], use = "p")
moduleTraitPvalue2[[set]] = corPvalueFisher(moduleTraitCor2[[set]], exprSize$nSamples[set])
}


pdf(file = "ModuleTraitRelationships-HA.pdf", wi = 7, he = 10);
# Plot the module-trait relationship table for set number 1
set = 1
textMatrix = paste(signif(moduleTraitCor2[[set]], 2), "\n(",
signif(moduleTraitPvalue2[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor2[[set]])
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor2[[set]],
xLabels = names(datTreatment2[[set]]),
yLabels = MEColorNames,
ySymbols = MEColorNames,
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module--trait relationships in", setLabels[set]))
dev.off()

pdf(file = "ModuleTraitRelationships-RHA.pdf", wi = 7, he = 10);
# Plot the module-trait relationship table for set number 1
set = 2
textMatrix = paste(signif(moduleTraitCor2[[set]], 2), "\n(",
signif(moduleTraitPvalue2[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor2[[set]])
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor2[[set]],
xLabels = names(datTreatment2[[set]]),
yLabels = MEColorNames,
ySymbols = MEColorNames,
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module--trait relationships in", setLabels[set]))
dev.off()


#############################
###### CONSENSUS CORRS ######
#############################

# For each module-trait pair, take the correlation that has the lower absolute value in the two sets 
# if the two correlations have the same sign, and zero relationship if the two correlations have opposite signs

# Initialize matrices to hold the consensus correlation and p-value
consensusCor2 = matrix(NA, nrow(moduleTraitCor2[[1]]), ncol(moduleTraitCor2[[1]]))
consensusPvalue2 = matrix(NA, nrow(moduleTraitCor2[[1]]), ncol(moduleTraitCor2[[1]]))

# Find consensus negative correlations
negative = moduleTraitCor2[[1]] < 0 & moduleTraitCor2[[2]] < 0
consensusCor2[negative] = pmax(moduleTraitCor2[[1]][negative], moduleTraitCor2[[2]][negative])
consensusPvalue2[negative] = pmax(moduleTraitPvalue2[[1]][negative], moduleTraitPvalue2[[2]][negative])

# Find consensus positive correlations
positive = moduleTraitCor2[[1]] > 0 & moduleTraitCor2[[2]] > 0
consensusCor2[positive] = pmin(moduleTraitCor2[[1]][positive], moduleTraitCor2[[2]][positive])
consensusPvalue2[positive] = pmax(moduleTraitPvalue2[[1]][positive], moduleTraitPvalue2[[2]][positive])

textMatrix = paste(signif(consensusCor2, 2), "\n(",
signif(consensusPvalue2, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor2[[set]])

pdf(file = "ModuleTraitRelationships-consensus.pdf", wi = 7, he = 10)
labeledHeatmap(Matrix = consensusCor2,
xLabels = names(datTreatment2[[set]]),
yLabels = MEColorNames,
ySymbols = MEColorNames,
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Consensus module--trait relationships across\n",
paste(setLabels, collapse = " and ")))
dev.off()


#############################
###### GENOTYPE CORRS #######
#############################

# to remove modules that have a significant effect of genotype, will take the opposite approach-
# will remove based on the highest correlations across the two sets
# if the highest correlations are 

# first, save correlations locally:
save(moduleTraitCor2, moduleTraitPvalue2, file = "Module-Trait-Correlations.RData")
load("ResultsFiles/Coexpression/Module-Trait-Correlations.RData")

TraitCor <- lapply(moduleTraitCor2, function(x) {as.data.frame(x)})
TraitCor[[1]]$Network <- "HA"
TraitCor[[2]]$Network <- "RHA"

TraitP <- lapply(moduleTraitPvalue2, function(x) {as.data.frame(x)})
TraitP[[1]]$Network <- "HA"
TraitP[[2]]$Network <- "RHA"

# maximum correlation for accession
MaxAccCor <- lapply(TraitCor, function(x) {
  apply(x[,c(3:5)], 1, function(x) {max(abs(x))})
})

# get corresponding p-value
apply(MaxAccCor[[1]], 2, function(x) {which(TraitCor[[1]][,c(3:5)] == x)})


MaxAccCor <- lapply(TraitCor, function(x) {
  apply(x[,c(3:5)], 1, function(x) {max(abs(x))})
})

##############
apply(TraitP[[1]][,c(3:5)], 1, FUN=min)

TraitCor[[1]]$MinAccP <- pmin(TraitP[[1]][,c(3:5)])
# minimum pvalues + maximum Correlation coefficients for accession

moduleTraitPvalue2[[1]]$AccP <- pmin(moduleTraitPvalue2[[1]][,c(3:5)])
moduleTraitCor2[[1]]$Accession



# save only treatment values:
moduleTreatCor <- lapply(moduleTraitCor2, function(x) {x[,c("LowNut", "HighSalt", "LowNut:HighSalt")]})
moduleTreatP <- lapply(moduleTraitPvalue2, function(x) {x[,c("LowNut", "HighSalt", "LowNut:HighSalt")]})

# correlation for accession
# not requiring the same direction
moduleAccCor <- lapply(moduleTraitCor2, function(x) {x[,c(3:5)]})
moduleAccP <- lapply(moduleTraitPvalue2, function(x) {x[,c(3:5)]})



SigAccession = moduleAccP[[1]] < 0.05 | moduleAccP[[2]] < 0.05




SigAccessionTest = ifelse(moduleAccP[[1]] < 0.05 | moduleAccP[[2]] < 0.05,
                          abs(pmax(moduleAccCor)), "NA")

#############
consensusTreatCorSameDir <- negativeTreat == "TRUE" | positiveTreat == "TRUE"


###### Consensus results for treatment effects using the conservative approach
# Initialize matrices to hold the consensus correlation and p-value
consensusTreatCor = matrix(NA, nrow(moduleTreatCor[[1]]), ncol(moduleTreatCor[[1]])) #89 x 3
consensusTreatP = matrix(NA, nrow(moduleTreatCor[[1]]), ncol(moduleTreatCor[[1]])) # 89 x 3



# Find consensus negative correlations
negativeTreat = moduleTreatCor[[1]] < 0 & moduleTreatCor[[2]] < 0
consensusTreatCorNeg = pmax(moduleTreatCor[[1]][negativeTreat], moduleTreatCor[[2]][negativeTreat])

consensusTreatCor[negativeTreat] = pmax(moduleTreatCor[[1]][negativeTreat], moduleTreatCor[[2]][negativeTreat])
consensusTreatP[negativeTreat] = pmax(moduleTreatP[[1]][negativeTreat], moduleTreatP[[2]][negativeTreat])

# Find consensus positive correlations
positiveTreat = moduleTreatCor[[1]] > 0 & moduleTreatCor[[2]] > 0
consensusTreatCor[positiveTreat] = pmin(moduleTreatCor[[1]][positiveTreat], moduleTreatCor[[2]][positiveTreat])
consensusTreatP[positiveTreat] = pmax(moduleTreatP[[1]][positiveTreat], moduleTreatP[[2]][positiveTreat])

hist(consensusTreatP)
length(which(consensusTreatP < 0.05)) #54

str(consensusTreatCor)


consensusTreatCorDF <- as.data.frame(as.table(consensusTreatCor)) #doesn't work

# save only accession values:
head(consensusTreatCor, 2)

textMatrix = paste(signif(consensusCor2, 2), "\n(",
                   signif(consensusPvalue2, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor2[[set]])


####################

consensusTreatCorNeg = pmax(moduleTreatCor[[1]][negativeTreat], moduleTreatCor[[2]][negativeTreat])

consensusTreatCor[negativeTreat] = pmax(moduleTreatCor[[1]][negativeTreat], moduleTreatCor[[2]][negativeTreat])
consensusTreatP[negativeTreat] = pmax(moduleTreatP[[1]][negativeTreat], moduleTreatP[[2]][negativeTreat])




####################

# Initialize matrices to hold the consensus correlation and p-value
consensusCor2 = matrix(NA, nrow(moduleTraitCor2[[1]]), ncol(moduleTraitCor2[[1]]))
consensusPvalue2 = matrix(NA, nrow(moduleTraitCor2[[1]]), ncol(moduleTraitCor2[[1]]))

# Find consensus negative correlations
negative = moduleTraitCor2[[1]] < 0 & moduleTraitCor2[[2]] < 0
consensusCor2[negative] = pmax(moduleTraitCor2[[1]][negative], moduleTraitCor2[[2]][negative])
consensusPvalue2[negative] = pmax(moduleTraitPvalue2[[1]][negative], moduleTraitPvalue2[[2]][negative])

# Find consensus positive correlations
positive = moduleTraitCor2[[1]] > 0 & moduleTraitCor2[[2]] > 0
consensusCor2[positive] = pmin(moduleTraitCor2[[1]][positive], moduleTraitCor2[[2]][positive])
consensusPvalue2[positive] = pmax(moduleTraitPvalue2[[1]][positive], moduleTraitPvalue2[[2]][positive])
