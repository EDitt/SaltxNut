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

load("ModuleFiles_100820.RData")
dim(MEs) #102 x 71 -> eigengenes values for all samples across all 71 modules

length(moduleLabels) #23236 (assignments for each gene?)
length(moduleColors) #23236 (module color assignments for each gene?)

str(geneTree) # co-expression network

load("NetworkData_outliersRem.RData")
dim(y_red) # 102 x 13
dim(datExpr_red) # 102  23236

# associated metadata
metadata_cult <- read.csv("StudyDesign_Inbred_noOut.csv", header=T)
head(metadata_cult)
dim(metadata_cult) 

sampleNums = rownames(datExpr_red) #N=102
metadata_cult_red <- metadata_cult[which(metadata_cult$Plant %in% sampleNums),]
rownames(metadata_cult_red) <- metadata_cult_red$Plant


#############################
####### TOP HUB GENES ####### 
#############################

softPower = 12

topHub <- chooseTopHubInEachModule(
datExpr_red,
moduleColors,
omitColors = "grey",
power = softPower,
type = "signed",
corFnc = bicor)

str(topHub)
class(topHub)

topHubdf <- data.frame(matrix(unlist(topHub ), nrow=length(topHub), byrow=T), stringsAsFactors=FALSE)
topHubdf$Modules <- names(topHub)
colnames(topHubdf)[1] <- c("Gene")

write.csv(topHub, file = "topHubGenes.csv")

#############################
# MODULE-TRAIT CORRELATIONS #
#############################

### Can correlate eigengenes with external traits to look for most significant associations:

# Define numbers of genes and samples
nGenes = ncol(datExpr_red)
nSamples = nrow(datExpr_red)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr_red, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

length(metadata_cult$Plant)

#TreatmentMod <- ~ 0 + Treatment
FullMod <- ~ 0 + Treatment + Group + Group:Accession + Reproductive + SampleDay

Treatments <- Module_Trait_Corr(metadata_cult_red, FullMod, MEs, nSamples)

pdf(file = "Treatments_pamStageTrue.pdf", width = 8, height = 8);
# Will display correlations and their p-values
textMatrix = paste(signif(Treatments$ModTraitCorr, 2), "\n(",
                   signif(Treatments$ModTraitP, 1), ")", sep = "");
dim(textMatrix) = dim(Treatments$ModTraitCorr)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = Treatments$ModTraitCorr,
               xLabels = names(Treatments$Traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#############################
# GET GENE MODULE ASSIGNMENTS 
#############################

modGenes <- as.data.frame(cbind(moduleColors, colnames(datExpr_red)))
colnames(modGenes)[2] <- c("Gene")
#check
aggregate(modGenes$Gene, by=list(modGenes$moduleColors), length)
table(dynamicColors)

write.csv(modGenes, 'GeneMod_Assignments_pamTRUE.csv')

#to get info for individual modules:
#modGenespink <- moduleColors=="pink"
#pinkGenes = datExpr_red[modGenespink]

################################
###### MODULE MEMBERSHIP #######
################################

### For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile.
# this allows us to quantify the similarity of all genes to every module

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr_red, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

write.csv(geneModuleMembership, 'GeneModMembership.csv')
write.csv(MMPvalue, 'GeneModMembership_Pvalues.csv')

#######################################
########## GENE SIGNIFICANCE ##########
#######################################

### Can quantify associations with traits of interest by defining Gene Significance (GS) as
# the absolute value of the correlation between the gene and trait

#make numeric
metadata_cult_red$Nut <- as.factor(ifelse(metadata_cult_red$Treatment == "LowNut" | metadata_cult_red$Treatment == "Combo", "1", "0"))
metadata_cult_red$Salt <- as.factor(ifelse(metadata_cult_red$Treatment == "HighSalt" | metadata_cult_red$Treatment == "Combo", "1", "0"))

## Nutrient
geneNutSignificance = as.data.frame(cor(datExpr_red, metadata_cult_red$Nut, use = "p"))
GSPvalue_Nut = as.data.frame(corPvalueStudent(as.matrix(geneNutSignificance), nSamples))

Nutrient = as.data.frame(metadata_cult_red$Nut)
names(Nutrient) = "Nutrient"

names(geneNutSignificance) = paste("GS.", names(Nutrient), sep="")
names(GSPvalue_Nut) = paste("p.GS.", names(Nutrient), sep="")

## Salt
geneSaltSignificance = as.data.frame(cor(datExpr_red, metadata_cult_red$Salt, use = "p"))
GSPvalue_Salt = as.data.frame(corPvalueStudent(as.matrix(geneSaltSignificance), nSamples))

Salt = as.data.frame(metadata_cult_red$Salt)
names(Salt) = "Salt"

names(geneSaltSignificance) = paste("GS.", names(Salt), sep="")
names(GSPvalue_Salt) = paste("p.GS.", names(Salt), sep="")

geneTreatmentSig = data.frame(geneNutSignificance, GSPvalue_Nut,
                      		geneSaltSignificance, GSPvalue_Salt)

write.csv(geneTreatmentSig, 'Gene_TreatmentSig.csv')

#############################
# GENES WITH HIGH GS and MM # ### Example
#############################

module = "bisque4"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
pdf(file = "Nutrient_bisque4_pamFALSE.pdf", width = 8, height = 8);
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneNutSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Low Nutrient",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

#############################
##### GENE CONNECTIVITY #####
#############################

ADJ1=abs(cor(datExpr_red,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ1, moduleColors)
head(Alldegrees1)

#The function intramodularConnectivity computes the whole network connectivity kTotal,
#the within module connectivity kWithin, kOut=kTotal-kWithin, and kDiff=kIn-kOut=2*kIN-kTotal

### See III. 7: Module membership, intramodular connectivity, and screening for intramocular hub genes for more

write.csv(Alldegrees1, 'Gene_Connectivity.csv')

#############################
###### OUTPUT GENE INFO #####
#############################

geneInfo = data.frame(geneModuleMembership, MMPvalue, 
                      Alldegrees1,
                      geneNutSignificance, GSPvalue_Nut,
                      geneSaltSignificance, GSPvalue_Salt)

#write.csv(geneInfo, file = "geneModuleInfo_pamFALSE.csv")
write.csv(geneInfo, file = "geneModuleInfo_pamTRUE.csv")


###############################
# HETEROTIC GROUP ASSOCIATION #
###############################

GroupMod <- ~ 0 + Group

Groups <- Module_Trait_Corr(metadata_cult_red, GroupMod, MEs, nSamples)

pdf(file = "Groups_pamStageTrue.pdf", width = 8, height = 8);
# Will display correlations and their p-values
textMatrix = paste(signif(Groups$ModTraitCorr, 2), "\n(",
                   signif(Groups$ModTraitP, 1), ")", sep = "");
dim(textMatrix) = dim(Groups$ModTraitCorr)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = Groups$ModTraitCorr,
               xLabels = names(Groups$Traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-Heterotic Group"))
dev.off()
