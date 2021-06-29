# Module categories

#########################
######## SETUP ##########
#########################

source("Functions.R")
library(tidyr)
library(UpSetR)

# load ANOVA results
LM_Results <- read.csv("ResultsFiles/Coexpression/Module_Anova.csv", header=T) # 88 modules

## R^2 categories for simple models:
load("ResultsFiles/Coexpression/QC/Module_RsquaredCats.RData")

# # Will remove modules with R^2 values higher than 0.95 in models with Accession as only explanatory factor
LM_ResultsSub <- LM_Results[which(!LM_Results$Module %in% ModuleRsquaredList$AccessionR2_more95),] # now 54 modules

##############################
# SIG MODULES BY MODEL FACTOR #
##############################

critp <- 0.05/length(LM_ResultsSub$Module)

### modules that are significant for Nutrient Main Effect:
LR_Nut_Sig <- LM_ResultsSub[which(LM_ResultsSub$Nut_p < critp), "Module"] #40 (52 before removing 34 modules)

### modules that are significant for Salt Main Effect:
LR_Salt_Sig <- LM_ResultsSub[which(LM_ResultsSub$Salt_p < critp), "Module"] #15

### modules that are significant for Nut x Salt Interaction:
LR_NutxSalt_Sig <- LM_ResultsSub[which(LM_ResultsSub$NutxSalt_p < critp), "Module"] #26

VennNumbers(LR_Nut_Sig, LR_Salt_Sig, LR_NutxSalt_Sig)
# 6 sig for all
# 7 Nut & Salt Main Effects
# 0 Salt & Nut x Salt
# 19 Nut & Nut x Salt
# 8 Nut (was 20 before removing 34 modules from dataset)
# 2 Salt
# 1 NutxSalt

### modules that are significant for Accession:
LR_Accession_Sig <- LM_ResultsSub[which(LM_ResultsSub$Accession_p < critp), "Module"] #29 (was 62 before removing 34 modules)

######
### modules that are significant for Bench: (changed code & did not keep this factor in spreadsheet)
#LR_bench <- lapply (LR_mod_ANOVA, function(x) {which(x["Bench", "Pr(>F)"] < critp) })
lapply (LR_mod_ANOVA, function(x) {x["Bench", "F value"]})
#length(which(LR_bench==1)) #9: darkviolet, lightcoral, turquoise, navajowhite1, thistle2, green, navajowhite2, steelblue, darkorange
# navajowhite2, green > 8, turqoise, lightcoral, darkviolet >7
#### MEgreen: Bench is the highest F-value

############################
# MODEL FACTOR MOD OVERLAP #
############################

AllSigTreatMods <- union(LR_NutxSalt_Sig, union(LR_Nut_Sig, LR_Salt_Sig))
length(AllSigTreatMods) #43 (was 55 before removing 34 modules)

SigModOverlap <- GeneSets(LR_Nut_Sig, LR_Salt_Sig, LR_NutxSalt_Sig)

lapply(SigModOverlap, function(x) {length(x)})

# how many of the modules in these categories are also significant for Accession?
AccessionOverlap <- lapply(SigModOverlap, function(x) {intersect(x, LR_Accession_Sig)})
lapply(AccessionOverlap, function(x) {length(x)})
# 4 sig for all
# 4 Nut & Salt Main Effects
# 0 Salt & Nut x Salt
# 10 Nut & Nut x Salt
# 4 Nut (was 15)
# 1 Salt
# 0 NutxSalt

############################
######## UPSET PLOT ########
############################

LR_list <- list(NutxSalt = LR_NutxSalt_Sig, 
                NutMain = LR_Nut_Sig, 
                SaltMain = LR_Salt_Sig,
                Accession = LR_Accession_Sig
)

upset(fromList(LR_list),
      #keep.order = TRUE,
      order.by = "freq",
      #group.by = "sets", 
      nsets = 13,
      #empty.intersections = "on",
      nintersects = 30)

#################################
### HOLM/HOCHBERG CORRECTION ? ###
#################################

critp2 <- 0.05/length(AllSigTreatMods) 
# sidak correction 
critp2 <- 1-(1-0.05)^(1/length(AllSigTreatMods)) # 0.00119

# take only modules for which there are significant treatment effects
LM_Results_Sig <- LM_ResultsSub[which(LM_ResultsSub$Module %in%
                                        AllSigTreatMods),]
# different power for Nutrient and Salt (diff distribution of p-values)

# Holm's Step-Down Procedure
NutDE_Sig <- Holm_Hochberg("Difference_p.DE_Nut", LM_Results_Sig)
# Holm's or Hochberg would give same result here (depends on which way sorted)
# first n.s. is 0.2

SaltDE_Sig <- Holm_Hochberg("Difference_p.DE_Salt", LM_Results_Sig)
# first n.s. is 0.00175

#################################
# SIG MODULES BY PAIRWISE DIFFS #
#################################

# Differences between Nutrient and Control
#DE_Nut_sig <- LM_Results_Sig[which(LM_Results_Sig$Difference_p.DE_Nut < critp2),"Module"] #38 (was 43)
#DE_Nut_sig <- LM_Results_Sig[which(LM_Results_Sig$Difference_p.DE_Nut < 0.05),"Module"] #41
DE_Nut_sig <- LM_Results_Sig[which(LM_Results_Sig$Difference_p.DE_Nut < 0.01),"Module"] #40

# Differences between Salt and Control
#DE_Salt_sig <- LM_Results_Sig[which(LM_Results_Sig$Difference_p.DE_Salt < critp2),"Module"] #10
#DE_Salt_sig <- LM_Results_Sig[which(LM_Results_Sig$Difference_p.DE_Salt < 0.05),"Module"] #24
DE_Salt_sig <- LM_Results_Sig[which(LM_Results_Sig$Difference_p.DE_Salt < 0.01),"Module"] #18

# Differences between Combo and Control
#DE_Combo_sig <- LM_Results_Sig[which(LM_Results_Sig$Difference_p.DE_Combo < critp2),"Module"] #25
#DE_Combo_sig <- LM_Results_Sig[which(LM_Results_Sig$Difference_p.DE_Combo < 0.05),"Module"] #32
DE_Combo_sig <- LM_Results_Sig[which(LM_Results_Sig$Difference_p.DE_Combo < 0.01),"Module"] #28

length(union(DE_Combo_sig, union(DE_Salt_sig, DE_Nut_sig))) #41 (#43 with p=0.05 & p=0.01)

DE_Overlaps <- GeneSets(DE_Combo_sig,
                        DE_Salt_sig,
                        DE_Nut_sig)
lapply(DE_Overlaps, function(x) {length(x)})

# in common all: 6 (was 18)
# combo-salt only: 1 (was 0)
# salt-nut only: 1 (was 4)
# combo-nut only: 18 (was 22)
# combo only: 0 (was 0)
# salt only: 2 (was 2)
# nut only: 13 (was 9)

# for p=0.05
# in common all: 18
# combo-salt only: 0
# salt-nut only: 4
# combo-nut only: 14
# combo only: 0
# salt only: 2
# nut only: 5

# for p=0.01
# in common all: 12
# combo-salt only: 1
# salt-nut only: 3
# combo-nut only: 15
# combo only: 0
# salt only: 2
# nut only: 10

### use p < 0.01

##############################
##### DIFF B/W TREATMENTS ####
##############################

Salt_Nut_sig <- LM_Results_Sig[which(LM_Results_Sig$Difference_p.Nut.Salt < 0.01),"Module"] #37

Combo_Nut_sig <- LM_Results_Sig[which(LM_Results_Sig$Difference_p.Nut.Combo < 0.01),"Module"] #30

Combo_Salt_sig <- LM_Results_Sig[which(LM_Results_Sig$Difference_p.Salt.Combo < 0.01),"Module"] #20

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

SigDiffOverlap <- GeneSets(Combo_Nut_sig, Combo_Salt_sig, Salt_Nut_sig)

############################
#### SAVE MODULE LISTS #####
############################

AllSigMods <- union(DE_Combo_sig, union(DE_Salt_sig, DE_Nut_sig))
length(AllSigMods) # 43

save(SigModOverlap, DE_Overlaps, SigDiffOverlap, 
     LR_list, lm_list, DE_Combo_sig, DE_Salt_sig, DE_Nut_sig,
     file="ResultsFiles/Coexpression/LR_SigModOverlap.RData")

############################
### INTERESTING MODULES ###
############################

### Shared responses:

# are there any modules that are DE in all & have no significant pairwise differences?
AllSigDiffs <- union(Combo_Nut_sig, union(Combo_Salt_sig, Salt_Nut_sig))

setdiff(DE_Overlaps$InCommonAll, AllSigDiffs) # blue2, darkorange2, mediumpurple2, salmon: All have significant model effects Nutrient + Nut x Salt
# blue2/darkorange2/salmon are all more highly correlated with each other than mediumpurple2

# what pairwise responses are present?
lapply(SigDiffOverlap, function(x) {intersect(x, DE_Overlaps$InCommonAll)})

### Nutrient responses that are not affected by the presence of salt:
#Nut_unconditional <- 
setdiff(DE_Overlaps$DE_Combo_sigDE_Nut_sigOnly, Combo_Nut_sig) # only 1: yellow4 (nutrient main effect significant)

### Nutrient responses that *are* affected by the presence of salt:
intersect(DE_Overlaps$DE_Combo_sigDE_Nut_sigOnly, Combo_Nut_sig) # N=14

### Salt responses not affected by nutrients:
setdiff(DE_Overlaps$DE_Combo_sigDE_Salt_sigOnly, Combo_Salt_sig) # 1: lightcyan1 (only salt/combo one)
intersect(DE_Overlaps$DE_Combo_sigDE_Salt_sigOnly, Combo_Salt_sig) # none-
intersect(DE_Overlaps$DE_Salt_sigOnly, Combo_Salt_sig)  # greenyellow (black is the other one)



######## Plot raw data
AllData_ordered <- AllData[order(AllData$Treatment, AllData$Accession, decreasing = FALSE),]
barplot(AllData_ordered$MEwhite, names.arg = AllData_ordered$Treatment)

######### Scratch

LM_ResultsSub[which(LM_ResultsSub$Module %in%
                      DE_Overlaps$DE_Salt_sigDE_Nut_sigOnly), 
              c(1,18,20,21,23,24,26)]
# darkseagreen4: Nut (p<0.0001), Salt (p=0.007)
# darkslateblue: Nut (p<0.0001), Salt (p=0.03)
# grey60: Nut (p<0.0001), Salt (p=0.0066)
# plum1: Nut (p<0.0001), Salt (p=0.0007)

############################
# PLOT SIGNIFICANT MODULES #
############################

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


