# likelihood ratio tests

#########################
######## SETUP ##########
#########################

library(dplyr)
library(car)
source("Functions.R")
library(tidyr)
library(emmeans)

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

AllData <- merge(Design[,c(2:8, 14, 17, 18, 23, 31:33, 42, 43)], 
                 ModEGs, by="Plant")
str(AllData)

levels(AllData$Osmocote)
levels(AllData$Salt)

# save spreadsheet
write.csv(AllData, file="ResultsFiles/Coexpression/Module_EGs.csv", row.names = FALSE)

################################
### MODEL 2 TREATMENT FACTORS ###
################################

### treatment 2 factors with 2 levels each

LR_Mod <- function(Yvar, dataset) {
  mod <- lm(Yvar ~ Osmocote + Salt +
                Osmocote:Salt +
                Accession +
                Osmocote:Accession +
                Salt:Accession +
                Osmocote:Salt:Accession +
                Bench, data = dataset,
              contrast=list(Accession=contr.sum, Bench=contr.sum))
  return(mod)
}


##############################
###### ANALYZE MODULES #######
##############################

LR_mod_results <- lapply (AllData[,c(17:104)], function(x) {LR_Mod(x, AllData) })
str(LR_mod_results)

### ANOVA
LR_mod_ANOVA <- lapply (LR_mod_results, function(x) {as.data.frame(Anova(x, test="F", type=2))})


##############################
#### SAVE F and P-VALUES #####
##############################

# select relevant info (F values & p values for factors of interest)
Anova_columns <- lapply(LR_mod_ANOVA, function(x) {x[c(1:3,5),c(3,4)]})
#Anova_columns <- lapply(LR_mod_ANOVA, function(x) {cbind(rownames(x[c(1:3,5),]), 
#                                                         x[c(1:3,5),c(3,4)])})

# transpose
Anova_wide <- lapply(Anova_columns, function(x) {as.data.frame(t(x))})

# add column for module identity
Anova_widewLabels <-
  lapply(names(Anova_wide), function(x) { Anova_wide[[x]]$Module <- x;return(Anova_wide[[x]])})

# F values
Anova_Fvals <- lapply(Anova_widewLabels, function(x) {x[1,]})

#combine into 1 df
All_Fvals <- do.call("rbind", Anova_Fvals)
colnames(All_Fvals) <- c("Nut_F", "Salt_F", "Accession_F", "NutxSalt_F", "Module")

# p values
Anova_pvals <- lapply(Anova_widewLabels, function(x) {x[2,]})

#combine into 1 df
All_pvals <- do.call("rbind", Anova_pvals)
colnames(All_pvals) <- c("Nut_p", "Salt_p", "Accession_p", "NutxSalt_p", "Module")

# merge
Anova_results <- merge(All_Fvals, All_pvals, by="Module")

##############################
########## LS MEANS ##########
##############################

# going to calculate for all modules
Mod_means <- lapply (LR_mod_results, 
                             function(x) {emmeans(x, ~ Osmocote*Salt, type = "response")})
# dataframe of all the means
Mod_means_df <- lapply(Mod_means, function(x) {as.data.frame(x)[,c(1:4)]})

Mod_meanswLabels <-
  lapply(names(Mod_means_df), function(x) { Mod_means_df[[x]]$Module <- x;return(Mod_means_df[[x]])})

All_ModMeans <- do.call("rbind", Mod_meanswLabels)

# add treatment column
All_ModMeans$Treatment <- ifelse(All_ModMeans$Osmocote == "High" &
                                   All_ModMeans$Salt == "NoSalt",
                                 "Control", ifelse(All_ModMeans$Osmocote == "High" &
                                                     All_ModMeans$Salt == "Salt",
                                                   "Salt", ifelse(All_ModMeans$Osmocote == "Low" &
                                                                    All_ModMeans$Salt == "NoSalt",
                                                                  "LowNut", "Combo")))


# check
aggregate(All_ModMeans$Module, by=list(All_ModMeans$Osmocote,
                                       All_ModMeans$Salt, All_ModMeans$Treatment), length)

# wide format
All_ModMeans_wide <- reshape(All_ModMeans[,c(3:6)], idvar = "Module", timevar = "Treatment",
                             direction = "wide")


##############################
###### PAIRWISE DIFFS ########
##############################

Mod_means_pairwisediffs <- lapply (Mod_means, 
                                     function(x) {as.data.frame(pairs(x))})

# better contrast labels
Contrast_labels <- c("DE_Nut", "DE_Salt", "DE_Combo", "Nut-Salt", "Nut-Combo", "Salt-Combo")

# dataframe of pairwise diffs
Mod_means_pairwisediffs_df <- lapply(Mod_means_pairwisediffs, function(x) {
  cbind(as.data.frame(x)[,c(2,3,6)], Contrast_labels)})

# add module labels
PairwiseDiffs_wModLabels <-
  lapply(names(Mod_means_pairwisediffs_df), function(x) { 
    Mod_means_pairwisediffs_df[[x]]$Module <- x;return(Mod_means_pairwisediffs_df[[x]])})

# combine into 1 dataframe
All_PairwiseDiffs <- do.call("rbind", PairwiseDiffs_wModLabels)
colnames(All_PairwiseDiffs) <- c("Difference", "Difference_SE", "Difference_p", "Contrast", "Module")

# wide format
All_PairwiseDiffs_wide <- reshape(All_PairwiseDiffs, idvar = "Module", timevar = "Contrast",
                             direction = "wide")

##############################
## COMBINE AND SAVE RESULTS ##
##############################

# df to save:
#Anova_results, All_ModMeans_wide, All_PairwiseDiffs_wide

LM_Results <- merge(Anova_results, merge(All_ModMeans_wide, All_PairwiseDiffs_wide, by="Module"),
                     by="Module")

# save
write.csv(LM_Results, file = "ResultsFiles/Coexpression/Module_Anova.csv", row.names = FALSE)

##############################
### LS MEANS PER ACCESSION ###
##############################



##############################
###### HIGHEST F VALUES ######
##############################


LR_highestF <- lapply (LR_mod_ANOVA, function(x) {rownames(x[which(x$`F value` == 
                                                         max(na.omit(x$`F value`))),])})

length(which(LR_highestF=="Osmocote")) #40
length(which(LR_highestF=="Accession")) #7 (34 additional which will be removed due to high R^2 for accession)
length(which(LR_highestF=="Salt")) #4
length(which(LR_highestF=="Bench")) #1 (MEgreen)
length(which(LR_highestF=="Osmocote:Salt")) #2 (salmon, greenyellow)
length(which(LR_highestF=="Osmocote:Accession")) #0
length(which(LR_highestF=="Salt:Accession")) #0
length(which(LR_highestF=="Osmocote:Salt:Accession ")) #0


