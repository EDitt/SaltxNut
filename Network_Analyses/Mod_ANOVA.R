# ANOVA with Treatment as 4-levels for posthoc (Wald?) tests

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

AllData$Osmocote <- factor(AllData$Osmocote)
levels(AllData$Osmocote) # high is baseline

AllData$Salt <- factor(AllData$Salt)
levels(AllData$Salt) # no salt is baseline

##############################
########### MODEL ############
##############################

lm_Mod <- function(Yvar, dataset) {
  mod <- lm(Yvar ~ Treatment +
               Accession + 
               Treatment:Accession +
               Bench, data = dataset,
             contrast=list(Treatment=contr.treatment, Accession=contr.sum, Bench=contr.sum))
  #result <- as.data.frame(Anova(mod, test="F", type=2))
  #SigValues <-
  #  data.frame("BenchSig" = result["Bench", "Pr(>Chi)"],
  #             "TreatmentxAccesionSig" = result["Treatment:Accession", "Pr(>Chi)"])
  return(mod)
}


##############################
###### ANALYZE MODULES #######
##############################

lm_mod_results <- lapply (AllData[,c(17:104)], function(x) {lm_Mod(x, AllData) })

critp <- 0.05/length(AllData[,c(17:104)])

### ANOVA
lm_mod_ANOVA <- lapply (lm_mod_results, function(x) {as.data.frame(Anova(x, test="F", type=2))})


### modules that are significant for Bench:
lmBench <- lapply (lm_mod_ANOVA, function(x) {which(x["Bench", "Pr(>F)"] < critp) })
length(which(lmBench==1)) #9
lmBenchSigMods <- names(which(lmBench==1))

### modules that are significant for Treatment x Accession:
lmTreatmentxAccession <- lapply (lm_mod_ANOVA, function(x) {which(x["Treatment:Accession", "Pr(>F)"] < critp) })
length(which(lmTreatmentxAccession==1)) #11
lmTreatmentxAccessionSigMods <- names(which(lmTreatmentxAccession==1))

### modules that are significant for Treatment:
lmTreatment <- lapply (lm_mod_ANOVA, function(x) {which(x["Treatment", "Pr(>F)"] < critp) })
length(which(lmTreatment==1)) #55
lmTreatmentSigMods <- names(which(lmTreatment==1))

### modules that are significant for Accession:
lmAccession <- lapply (lm_mod_ANOVA, function(x) {which(x["Accession", "Pr(>F)"] < critp) })
length(which(lmAccession==1)) #62
lmAccessionSigMods <- names(which(lmAccession==1))


VennNumbers(lmTreatmentSigMods, lmAccessionSigMods, lmTreatmentxAccessionSigMods)
# 11 shared by all
# 22 Treatment & Accession significant
# 22 Treatment significant only
# 29 Accession significant only

##############################
### TREATMENT SIGNIFICANCE ###
##############################

# pick modules for which Treatment F value is greater than Accession F value?

TreatmentF <- lapply (lm_mod_ANOVA, function(x) {which(x["Accession", "F value"] < x["Treatment", "F value"]) })
length(which(TreatmentF==1)) #41
Treatment_highFSigMods <- names(which(TreatmentF==1))

AccessionF <- setdiff(lmTreatmentSigMods, Treatment_highFSigMods) #18
## a few of these (which had similar F-values between Treatment and Accession, just slightly higher in Accession appeared to potentially be interesting)

# some still look important (and have relatively equal F-values)
# Accession F-value greater than 2x Treatment F-value?
AccessionF2 <- lapply (lm_mod_ANOVA, function(x) {which(x["Accession", "F value"] > 2*(x["Treatment", "F value"])) })
length(which(AccessionF2==1)) #41
AccessionF2Mods <- names(which(AccessionF2==1))

length(intersect(lmTreatmentSigMods, AccessionF2Mods)) #12

##############################
###### PAIRWISE DIFFS ########
##############################

# using lsmeans

mod_results_means <- lapply (lm_mod_results[lmTreatmentSigMods], 
                             function(x) {emmeans(x, ~ Treatment, type = "response")})
mod_results_pairwisediffs <- lapply (mod_results_means, 
                                     function(x) {as.data.frame(pairs(x))})

critp2 <- 0.05/length(lmTreatmentSigMods) ## do I want to correct again?

# Differences between Combo and Control
DE_ComboP <- lapply(mod_results_pairwisediffs, function(x) {x[1,"p.value"]})
#DE_Combo_sig <- names(which(DE_ComboP < critp2)) #30
DE_Combo_sig <- names(which(DE_ComboP < 0.05)) #40

# Differences between Salt and Control
DE_SaltP <- lapply(mod_results_pairwisediffs, function(x) {x[2,"p.value"]})
#DE_Salt_sig <- names(which(DE_SaltP < critp2)) #9
DE_Salt_sig <- names(which(DE_SaltP < 0.05)) #26

# Differences between Nutrient and Control
DE_NutP <- lapply(mod_results_pairwisediffs, function(x) {x[3,"p.value"]})
#DE_Nut_sig <- names(which(DE_NutP < critp2)) #42
DE_Nut_sig <- names(which(DE_NutP < 0.05)) #53

length(union(DE_Combo_sig, union(DE_Salt_sig, DE_Nut_sig))) #55

DE_Overlaps <- GeneSets(DE_Combo_sig,
                     DE_Salt_sig,
                     DE_Nut_sig)
lapply(DE_Overlaps, function(x) {length(x)})
# 20: in common all
# 0 combo-salt only
# 4 salt-nut only
# 20 combo-nut only
# 0 combo only
# 2 salt only
# 9 nut only

# which overlap with ones that have a stronger effect of Treatment than Accession?
lapply(DE_Overlaps, function(x) {length(intersect(x, Treatment_highFSigMods))})
# 16/20 in common all -4
# 4/4 salt-nut only
# 12/20 combo-nut only -8
# 2/2 salt only
# 3/9 nut only -6

# which ones have a significant effect of Accession at all?
lapply(DE_Overlaps, function(x) {length(intersect(x, lmAccessionSigMods))})
# 11/20 in common all
# 2/4 salt-nut only
# 13/20 combo-nut only
# 1/2 salt only
# 6/9 nut only


##############################
##### DIFF B/W TREATMENTS ####
##############################

Combo_SaltP <- lapply(mod_results_pairwisediffs, function(x) {x[4,"p.value"]})
#Combo_Salt_sig <- names(which(Combo_SaltP < critp2)) #10
Combo_Salt_sig <- names(which(Combo_SaltP < 0.05)) #30

Combo_NutP <- lapply(mod_results_pairwisediffs, function(x) {x[5,"p.value"]})
#Combo_Nut_sig <- names(which(Combo_NutP < critp2)) #24
Combo_Nut_sig <- names(which(Combo_NutP < 0.05)) #33

Salt_NutP <- lapply(mod_results_pairwisediffs, function(x) {x[6,"p.value"]})
#Salt_Nut_sig <- names(which(Salt_NutP < critp2)) #38
Salt_Nut_sig <- names(which(Salt_NutP < 0.05)) #48


lm_list <- list(Nut = DE_Nut_sig, 
                Salt = DE_Salt_sig, 
                Combo = DE_Combo_sig,
                #Accession = lmAccessionSigMods,
                Combo_Salt = Combo_Salt_sig,
                Combo_Nut = Combo_Nut_sig,
                Salt_Nut = Salt_Nut_sig)

upset(fromList(lm_list),
      #keep.order = TRUE,
      order.by = "freq",
      #group.by = "sets", 
      nsets = 13,
      #empty.intersections = "on",
      nintersects = 30)

##############################
######### CATEGORIES #########
##############################

PairwiseDiffs <- GeneSets(Combo_Salt_sig,
                          Combo_Nut_sig,
                          Salt_Nut_sig)

InCommonDiffs <- lapply(TreatDiffs, function(x) {intersect(x, Overlaps$InCommonAll)})

Allpairwisediffs <- union(Combo_Salt_sig,
                          union(Combo_Nut_sig,
                                Salt_Nut_sig))

### Nutrient only- shared with combo/diff from combo
### Salt only- shared with combo/diff from combo

### Nutrient & Salt- shared with combo (no diffs), cancel out

#########################
### NUTRIENT-SPECIFIC ###
#########################

# DE in Nutrient but not salt
Nut_specific <- setdiff(DE_Nut_sig, DE_Salt_sig) # N=29


NutCombo_Cats <- GeneSets(Nut_specific,
                          DE_Combo_sig,
                          Combo_Nut_sig)
NutCombo_Cats[c(1,2,4,5)]

NutCatNums <- lapply(NutCombo_Cats, function(x) {length(x)})
# Significant in both + pairwise difference: 12
# Nut & Combo, no pairwise difference: 8
# Nut & Combo-Nut difference: 3
# Nut only: 6 <- but no difference between combo-nut?

# which are also independent of accession effects?
NutCatNums_noAcc <- lapply(NutCombo_Cats, function(x) {setdiff(x, lmAccessionSigMods)})
# 5 shared + diff
# 2 shared + same
# 2 not shared, diff
# 1 not shared, same

# any patterns caused by accession significance?
lapply(NutCombo_Cats, function(x) {length(intersect(x, lmAccessionSigMods))})
# 7/12 shared + diff
# 6/8 shared + same
# 1/3 not shared, diff
# 5/6 not shared, no diff

lapply(NutCombo_Cats, function(x) {length(intersect(x, AccessionF))})
# 2/12 shared + diff
# 6/8 shared, same
# 1/3 not shared, diff
# 5/6 not shared, no diff

# intersect with modules that have a large (2x F-value) of Treatment F-value
lapply(NutCombo_Cats, function(x) {length(intersect(x, AccessionF2Mods))})
# 6/8 shared, same
# 4/6 Nut only, no difference

#########################
##### SALT-SPECIFIC #####
#########################

# DE in Salt but not nutrient
Salt_specific <- setdiff(DE_Salt_sig, DE_Nut_sig) # N=2

SaltCombo_Cats <- GeneSets(Salt_specific,
                          DE_Combo_sig,
                          Combo_Salt_sig)
SaltCombo_Cats[c(1,2,4,5)]
# 1 is not shared & diff
# 1 is not shared & same

lapply(SaltCombo_Cats, function(x) {length(intersect(x, lmAccessionSigMods))})
# shared & same has significant accession effect


#########################
#### NUT-SALT SHARED ####
#########################

Nut_Salt_shared <- intersect(DE_Salt_sig, DE_Nut_sig) #24

Nut_Salt_Cats <- GeneSets(Nut_Salt_shared,
                           DE_Combo_sig,
                          Salt_Nut_sig)
# 16: DE in all, diff between salt-nut
# 4: DE in all, *not* different with salt-nut
# 4: DE in nut/salt only, different between them
# 0: DE in nut/salt only, *not* different between them

lapply(Nut_Salt_Cats, function(x) {length(intersect(x, lmAccessionSigMods))})
# 9: DE in all, diff between salt-nut
# 2: DE in all, *not* diff b/w salt-nut
# 2: DE in nut/salt only, different between them

lapply(Nut_Salt_Cats, function(x) {length(intersect(x, AccessionF2Mods))})
# 2: DE in all, *not* diff b/w salt-nut

#####
test1Means <- emmeans(test1, ~ Treatment, type = "response")
pairwisediffs <- pairs(test1Means)

lm_mod_ANOVA$MEdarkmagenta
lm_mod_ANOVA$MEskyblue4
lm_mod_ANOVA$MEthistle1

#####



#### graphs to test results
AllData_ordered <- AllData[order(AllData$Treatment, AllData$Accession, decreasing = FALSE),]
barplot(AllData_ordered$MEplum3, names.arg = AllData_ordered$Treatment)

# Modules that would be excluded with the high F filter that look like an important treatment effect:
# thistle2
# navajowhite1
# blue
# lightsteelblue1

lm_mod_ANOVA$MEthistle2
lm_mod_ANOVA$MEnavajowhite1
lm_mod_ANOVA$MEblue
lm_mod_ANOVA$MElightsteelblue1
# for all of these the F values are high for Treatment and Accession


### modules I do want to exclude:

lm_mod_ANOVA$MEdarkmagenta


####
ls_means <- function(Yvar, dataset) {
  mod <- lm(Yvar ~ Treatment +
              Accession + 
              Bench, data = dataset)
  Means <- emmeans(mod, ~ Treatment, type = "response")
  pairwisediffs <- pairs(Means)
  return(as.data.frame(pairwisediffs))
}

mod_results_means <- lapply (AllData[,c(AllTreatmentSigMods)], 
                             function(x) {ls_means(x, AllData) })
# multiple comparisons
Critp2 <- 0.05/length(AllData[,c(AllTreatmentSigMods)])

DE_Combo <- lapply(mod_results_means, function(x) {x[1,"p.value"]})
#DE_Combo_sig <- names(which(DE_Combo < Critp2)) #31
DE_Combo_sig <- names(which(DE_Combo < 0.05))

DE_Salt <- lapply(mod_results_means, function(x) {x[2,"p.value"]})
#DE_Salt_sig <- names(which(DE_Salt < Critp2)) #8
DE_Salt_sig <- names(which(DE_Salt < 0.05))

DE_Nut <- lapply(mod_results_means, function(x) {x[3,"p.value"]})
#DE_Nut_sig <- names(which(DE_Nut < Critp2)) #47
DE_Nut_sig <- names(which(DE_Nut < 0.05))

OverlapNums <- VennNumbers(DE_Combo_sig,
                           DE_Salt_sig,
                           DE_Nut_sig)



#####
AllData_ordered <- AllData[order(AllData$Treatment, AllData$Accession, decreasing = FALSE),]
barplot(AllData_ordered$MElightcyan, names.arg = AllData_ordered$Treatment)


