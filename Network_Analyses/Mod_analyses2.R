# likelihood ratio tests

#########################
######## SETUP ##########
#########################

library(dplyr)
library(car)
source("Functions.R")
library(tidyr)
library(emmeans)
library(UpSetR)

#########################
###### READ DATA ########
#########################

ModEGs <- read.csv("/Users/emilydittmar/Google Drive/Active Projects/Transcriptomics_Exp/Analyses/WGCNA/ConsensusNetwork_June2021/ModEigengenes_Consensus.csv", header=T)
ModEGs <- read.csv("/Users/eld72413/Google Drive/Active Projects/Transcriptomics_Exp/Analyses/WGCNA/ConsensusNetwork_June2021/ModEigengenes_Consensus.csv", header=T)
colnames(ModEGs)[length(ModEGs)] <- "Plant" # last column is plant ID

ModExp <- read.csv("/Users/emilydittmar/Google Drive/Active Projects/Transcriptomics_Exp/Analyses/WGCNA/ConsensusNetwork_June2021/ModExpression_Consensus.csv", header=T)
ModExp <- read.csv("/Users/eld72413/Google Drive/Active Projects/Transcriptomics_Exp/Analyses/WGCNA/ConsensusNetwork_June2021/ModExpression_Consensus.csv", header=T)
colnames(ModExp)[length(ModExp)] <- "Plant"

Design <- read.csv("DataFiles/StudyDesign_Inbred_noOut.csv", header=T)

Design$Treatment <- relevel(as.factor(Design$Treatment), ref="Control")

levels(Design$Treatment)
Design$Accession <- as.factor(make.names(Design$Accession))
levels(Design$Accession)
Design$Group <- as.factor(Design$Group)
levels(Design$Group)
Design$SampleDay <- factor(Design$SampleDay)
levels(Design$SampleDay)
Design$Bench <- factor(Design$Bench)
levels(Design$Bench)
Design$Reproductive <- factor(Design$Reproductive)
levels(Design$Reproductive)

## samples per accssion:
aggregate(Design$Plant, by=list(Design$Accession, Design$Treatment), length)

#########################
###### DATA SETUP #######
#########################

AllData <- merge(Design[,c(2:8, 14, 17, 18, 23, 31:33, 42, 43)], ModEGs, by="Plant")
str(AllData)

################################
# LIKELIHOOD RATIO INTERACTION #
################################

### make treatment 2 factors with 2 levels each

AllData$Nut <- as.factor(ifelse(AllData$Treatment == "LowNut" | AllData$Treatment == "Combo", "1", "0"))
AllData$Salt <- as.factor(ifelse(AllData$Treatment == "HighSalt" | AllData$Treatment == "Combo", "1", "0"))

#check
aggregate(AllData$Plant, by=list(AllData$Nut, AllData$Salt, AllData$Treatment), length)

LR_Mod <- function(Yvar, dataset) {
  mod <- lm(Yvar ~ Nut + Salt + Nut:Salt +
               Accession + 
               Bench, data = dataset)
  result <- as.data.frame(drop1(mod, test = "Chi"))
  SigValues <-
  data.frame("AccessionSig" = result["Nut","Pr(>Chi)"], 
             "BenchSig" = result["Bench", "Pr(>Chi)"],
             "NutxSaltSig" = result["Nut:Salt", "Pr(>Chi)"])
  return(SigValues)
}

### change model to type iii to evaluate main effect in presence of interaction?
LR_Mod2 <- lm(MEdarkmagenta ~ Osmocote + Salt +
              Osmocote:Salt +
              Accession +
              Osmocote:Accession +
              Salt:Accession +
              Osmocote:Salt:Accession +
              Bench, data = AllData,
            contrast=list(Accession=contr.sum, Bench=contr.sum))

LR_mod_results <- lapply (AllData[,c(17:104)], function(x) {LR_Mod(x, AllData) })
str(LR_mod_results)

critp <- 0.05/length(AllData[,c(17:104)])

### modules that are significant for Bench:
LRBench <- lapply (LR_mod_results, function(x) {which(x$BenchSig < critp) })
length(which(LRBench==1)) #20
LRBenchSigMods <- names(which(LRBench==1))

### modules that are significant for Accession:
Accession <- lapply (LR_mod_results, function(x) {which(x$AccessionSig < critp) })
length(which(Accession==1)) #27
AccessionSigMods <- names(which(Accession==1))

### modules that are significant for Nut x Salt:
NutxSalt <- lapply (LR_mod_results, function(x) {which(x$NutxSaltSig < critp) })
length(which(NutxSalt==1)) #27
NutxSaltSigMods <- names(which(NutxSalt==1))

# how many modules significant for Nut x Salt are significant for accession and/or bench?
VennNumbers(NutxSaltSigMods, AccessionSigMods, BenchSigMods)
# 14 are significant for all, 13 are significant for Nut x Salt & Accession, 6 significant for bench only


################################
# LIKELIHOOD RATIO MAIN EFFECTS #
################################

### main effects
# take out Nut:Salt interaction

LR_Mod2 <- function(Yvar, dataset) {
  mod <- lm(Yvar ~ Nut + Salt + 
               Accession + 
               Bench, data = dataset)
  result <- drop1(mod, test = "Chi")
  main_result <- data.frame("LowNutSig" = result["Nut","Pr(>Chi)"], "HighSaltSig" = result["Salt","Pr(>Chi)"])
  return(main_result)
}

LR_mod_results2 <- lapply (AllData[,c(17:104)], function(x) {LR_Mod2(x, AllData) })
str(LR_mod_results2)

Nut <- lapply (LR_mod_results2, function(x) {which(x$LowNutSig < critp) })
length(which(Nut==1)) #54

Salt <- lapply (LR_mod_results2, function(x) {which(x$HighSaltSig < critp) })
length(which(Salt==1)) #13

SaltSigMods <- names(which(Salt==1))
NutSigMods <- names(which(Nut==1))

#########################
######## OVERLAP ########
#########################

AllTreatMods <- union(NutxSaltSigMods, union(NutSigMods, SaltSigMods))
length(AllTreatMods) #57

Overlap <- GeneSets(NutxSaltSigMods, NutSigMods, SaltSigMods)

lapply(Overlap, function(x) {length(x)})

# 4 in common all
# 21 NutxSalt & Nut
# 7 Nut & Salt
# 1 NutxSalt & Salt
# 1 NutxSalt only
# 22 Nut only
# 1 Salt only

# how many of the modules in these categories are also significant for Accession?

AccessionOverlap <- lapply(Overlap, function(x) {intersect(x, AccessionSigMods)})

lapply(AccessionOverlap, function(x) {length(x)})
# 4 in common all
# 21 NutxSalt & Nut
# 0 Nut & Salt
# 1 NutxSalt & Salt
# 1 NutxSalt
# 0 Nut only
# 0 Salt only

############################
######## UPSET PLOT ########
############################


LR_list <- list(NutxSalt = NutxSaltSigMods, 
                  NutMain = NutSigMods, 
                  SaltMain = SaltSigMods,
                  Accession = AccessionSigMods)

upset(fromList(LR_list),
      #keep.order = TRUE,
      order.by = "freq",
      #group.by = "sets", 
      nsets = 13,
      #empty.intersections = "on",
      nintersects = 30)



#######################


############################
# PLOT SIGNIFICANT MODULES #
############################

AllSigMods <- union(NutxSalt, union(NutMods, SaltMods))
length(AllSigMods) # 57

# transpose
SigModEGs <- t(ModEGs[,c(AllSigMods)])

consMETree = hclust(dist(SigModEGs), method = "average")
plot(consMETree, main = "Consensus clustering of consensus module eigengenes",
     xlab = "", sub = "")

# all modules
AllModEGs <- t(ModEGs[-c(89,90)])

AllconsMETree = hclust(dist(AllModEGs), method = "average")
plot(AllconsMETree, main = "Consensus clustering of all module eigengenes",
     xlab = "", sub = "")



