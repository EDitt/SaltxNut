# Module categories

#########################
######## SETUP ##########
#########################

source("Functions.R")
library(tidyr)
library(UpSetR)
library(ggfortify)

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
#critp <- 0.05/length(LM_Results$Module)

### modules that are significant for Nutrient Main Effect:
LR_Nut_Sig <- LM_ResultsSub[which(LM_ResultsSub$Nut_p < critp), "Module"] #40 (52 before removing 34 modules)

### modules that are significant for Salt Main Effect:
LR_Salt_Sig <- LM_ResultsSub[which(LM_ResultsSub$Salt_p < critp), "Module"] #15

### modules that are significant for Nut x Salt Interaction:
LR_NutxSalt_Sig <- LM_ResultsSub[which(LM_ResultsSub$NutxSalt_p < critp), "Module"] #27

VennNumbers(LR_Nut_Sig, LR_Salt_Sig, LR_NutxSalt_Sig)
# 6 sig for all
# 7 Nut & Salt Main Effects
# 1 Salt & Nut x Salt
# 19 Nut & Nut x Salt
# 8 Nut (was 20 before removing 34 modules from dataset)
# 1 Salt
# 1 NutxSalt

### modules that are significant for Accession:
LR_Accession_Sig <- LM_ResultsSub[which(LM_ResultsSub$GroupxAccession_p < critp), "Module"] #30 (was 62 before removing 34 modules)

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

# take only modules for which there are significant treatment effects
LM_Results_Sig <- LM_ResultsSub[which(LM_ResultsSub$Module %in%
                                         AllSigTreatMods),]

# check that p-value is significant for HA & RHA AND they are in the same direction
SigDiff <- function(Df, Treatment, Pcol, Pcrit){
   #meanColumnHA <- paste0("emmean.HA_", Treatment)
   #meanColumnRHA <- paste0("emmean.RHA_", Treatment)
   meanColumnHA <- paste0("Difference.HA_", Treatment)
   meanColumnRHA <- paste0("Difference.RHA_", Treatment)
   pColumnHA <- paste0("Difference_p.HA_", Pcol)
   pColumnRHA <- paste0("Difference_p.RHA_", Pcol)
   SameDir <- Df[which((Df[,meanColumnHA] > 0 &
                           Df[,meanColumnRHA] > 0) |
                          (Df[,meanColumnHA] < 0 &
                              Df[,meanColumnRHA] < 0)), "Module"]
   Sig <- Df[which(Df[,pColumnHA] < Pcrit &
                      Df[,pColumnRHA] < Pcrit),"Module"]
   DE_sig <- intersect(Sig, SameDir)
   return(DE_sig)
}


#DE_Salt_sig <- SigDiff(LM_Results_Sig, "DE_Salt", "DE_Salt", 0.05) #7
#DE_Nut_sig <- SigDiff(LM_Results_Sig, "DE_Nut", "DE_Nut", 0.05) #33 (31 at p=0.01)
#DE_Combo_sig <- SigDiff(LM_Results_Sig, "DE_Combo", "DE_Combo", 0.05) #17 (7 at p=0.01)

DE_Salt_sig <- LM_Results_Sig[which(LM_Results_Sig$Difference_p.DE_Salt < 0.01), "Module"]
DE_Nut_sig <- LM_Results_Sig[which(LM_Results_Sig$Difference_p.DE_Nut < 0.01), "Module"]
DE_Combo_sig <- LM_Results_Sig[which(LM_Results_Sig$Difference_p.DE_Combo < 0.01), "Module"]

length(union(DE_Combo_sig, union(DE_Salt_sig, DE_Nut_sig))) #43

DE_Overlaps <- GeneSets(DE_Combo_sig,
                        DE_Salt_sig,
                        DE_Nut_sig)
lapply(DE_Overlaps, function(x) {length(x)})

# lsmeans across groups
# in common all: 17
# combo-salt only: 1
# salt-nut only: 5
# combo-nut only: 15
# combo only: 0
# salt only: 1
# nut only: 4

# requiring same patterns in HA & RHA
# in common all: 1
# combo-salt only: 1
# salt-nut only: 1
# combo-nut only: 9
# combo only: 0
# salt only: 1 
# nut only: 22

##############################
##### DIFF B/W TREATMENTS ####
##############################


#Salt_Nut_sig <- SigDiff(LM_Results_Sig, "Nut.Salt", "Nut.Salt", 0.05) #27
#Combo_Nut_sig <- SigDiff(LM_Results_Sig, "Nut.Combo","Nut.Combo", 0.05) #15
#Combo_Salt_sig <- SigDiff(LM_Results_Sig, "Salt.Combo", "Salt.Combo", 0.05) #2

Salt_Nut_sig <- LM_Results_Sig[which(LM_Results_Sig$Difference_p.Nut.Salt < 0.01), "Module"] #37
Combo_Nut_sig <- LM_Results_Sig[which(LM_Results_Sig$Difference_p.Nut.Combo < 0.01), "Module"] #30
Combo_Salt_sig <- LM_Results_Sig[which(LM_Results_Sig$Difference_p.Salt.Combo < 0.01), "Module"] #21

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

setdiff(DE_Overlaps$InCommonAll, AllSigDiffs) 
# blue2, darkorange2, mediumpurple2, salmon

# what pairwise responses are present?
lapply(SigDiffOverlap, function(x) {intersect(x, DE_Overlaps$InCommonAll)}) 
# in common all: antiquewhite2, blue, cyan, floralwhite, orangered3
# nut-salt only: lightcyan1, orangered4, skyblue4, white

### Nutrient responses that are not affected by the presence of salt:
#Nut_unconditional <- 
setdiff(DE_Overlaps$DE_Combo_sigDE_Nut_sigOnly, Combo_Nut_sig) 
# yellow4

### Nutrient responses that *are* affected by the presence of salt:
intersect(DE_Overlaps$DE_Combo_sigDE_Nut_sigOnly, Combo_Nut_sig) # N=14

### Salt responses not affected by nutrients: no combo-salt only modules
setdiff(DE_Overlaps$DE_Combo_sigDE_Salt_sigOnly, Combo_Salt_sig) # 0
intersect(DE_Overlaps$DE_Combo_sigDE_Salt_sigOnly, Combo_Salt_sig) # 0
intersect(DE_Overlaps$DE_Salt_sigOnly, Combo_Salt_sig)  # greenyellow (other module is black)

### antagonistic effects
intersect(DE_Overlaps$DE_Salt_sigDE_Nut_sigOnly, Salt_Nut_sig) # plum1


############################
# PCA OF SIGNIFICANT MODS ##
############################

# load eigengenes
ModEGs <- read.csv("ResultsFiles/Coexpression/Module_EGs.csv", header=T)
ModEGs$Treatment <- factor(ModEGs$Treatment, levels=c("Control", "LowNut", "HighSalt", "Combo"))
levels(ModEGs$Treatment)

pca <- prcomp(ModEGs[,AllSigTreatMods], center = TRUE, scale. = TRUE)
summary(pca) # PC1=43.91%, PC2=15.10%

autoplot(pca, data=ModEGs, 
         colour='Treatment', shape='Accession', frame = TRUE) +
   scale_fill_manual(values = c("#4DAF4A", "#E41A1C", "#E6AB02", "#377EB8")) +
   scale_colour_manual(values = c("#4DAF4A", "#E41A1C", "#E6AB02", "#377EB8")) +
   theme_minimal()
ggsave("/Users/emilydittmar/Google Drive/Active Projects/Transcriptomics_Exp/Manuscript/SaltxNut/Figures/ModsPCA.png")


############################
######### BARPLOT ##########
############################

#my_dataSig <- lapply(colnames(LM_Results_Sig[,c(21,23,25,27,29,31,33,35,37,39,41,43)]), 
#                     function(x) {SigDEdf(LM_Results_Sig, PvaluesCol=x, CritP=0.05)})
#names(my_dataSig) <- colnames(LM_Results_Sig[,c(21,23,25,27,29,31,33,35,37,39,41,43)])

# unchanged nutrient-combo
Nut_Uncond <- setdiff(DE_Overlaps$DE_Combo_sigDE_Nut_sigOnly, Combo_Nut_sig) #7
# no longer significant
Nut_Cond_NS <- DE_Overlaps$DE_Nut_sigOnly #22

# conditional
Nut_Cond <- intersect(DE_Overlaps$DE_Combo_sigDE_Nut_sigOnly, Combo_Nut_sig) #2
# combo is reduced relative to nutrient
Nut_Cond_reduced <- intersect(Nut_Cond,
                             LM_Results_Sig[which(abs(LM_Results_Sig$Difference.HA_DE_Nut) >
                                                           abs(LM_Results_Sig$Difference.HA_DE_Combo)), "Module"]) # 2
Nut_Cond_increased <- intersect(Nut_Cond,
                             LM_Results_Sig[which(abs(LM_Results_Sig$Difference.HA_DE_Nut) <
                                                     abs(LM_Results_Sig$Difference.HA_DE_Combo)), "Module"]) #0
Nut_Cond_opposite <- setdiff(Nut_Cond, union(Nut_Cond_reduced, Nut_Cond_increased)) #0

# unchanged nutrient-combo
Salt_Uncond <- setdiff(DE_Overlaps$DE_Combo_sigDE_Salt_sigOnly, Combo_Salt_sig) #1
# no longer significant
Salt_Cond_NS <- DE_Overlaps$DE_Salt_sigOnly #1
# salt conditional
Salt_Cond <- intersect(DE_Overlaps$DE_Combo_sigDE_Salt_sigOnly, Combo_Salt_sig) #0

### shared responses
Nut_Unspec_Uncond <- setdiff(DE_Overlaps$InCommonAll, Combo_Nut_sig)
Salt_Unspec_Uncond <- setdiff(DE_Overlaps$InCommonAll, Combo_Salt_sig)

### Salt-Nut onkly
Nut_Unspec_NS <- DE_Overlaps$DE_Salt_sigDE_Nut_sigOnly
Salt_Unspec_NS <- DE_Overlaps$DE_Salt_sigDE_Nut_sigOnly

### make df
Labels <- c("Unchanged", "Decreased_mag",
            "Increased_mag", "Opposite_dir", "Not_DE")
NutSpecificDf <- data.frame(Labels, "Number" = c(length(Nut_Uncond),
                                                            length(Nut_Cond_reduced),
                                                            length(Nut_Cond_increased),
                                                            length(Nut_Cond_opposite),
                                                            length(Nut_Cond_NS)))
NutSpecificDf$Category <- "Nut_specific"
NutUnSpecificDf <- data.frame(Labels, "Number" = c(length(Nut_Unspec_Uncond),
                                                 0,
                                                 0,
                                                 0,
                                                 length(Nut_Unspec_NS)))
NutUnSpecificDf$Category <- "Nut_Unspecific"

SaltSpecificDf <- data.frame(Labels, "Number" = c(length(Salt_Uncond),
                                                 0,
                                                 0,
                                                 0,
                                                 length(Salt_Cond_NS)))
SaltSpecificDf$Category <- "Salt_specific"

SaltUnSpecificDf <- data.frame(Labels, "Number" = c(length(Salt_Unspec_Uncond),
                                                   0,
                                                   0,
                                                   0,
                                                   length(Salt_Unspec_NS)))
SaltUnSpecificDf$Category <- "Salt_Unspecific"


df_all_full <- rbind(NutSpecificDf, NutUnSpecificDf,
                     SaltSpecificDf, SaltUnSpecificDf,
                     make.row.names=FALSE)
df_all_full$Category <- factor(df_all_full$Category, levels=c("Nut_specific",
                                                              "Nut_Unspecific",
                                                              "Salt_Unspecific",
                                                              "Salt_specific"))

Opposite_col <- c("olivedrab4")
Unchanged_col <- c("#4A6990FF") #periwinkle blue
Increased_col <- c("#A73030FF")
Decreased_col <- c("#EFC000FF")
NotDE_col <- c("grey55")

str(df_all)

all_p <- ggplot(df_all_full, aes(fill=Labels, y=Number, x=Category))
all_p + geom_bar(position="stack", stat="identity", color="black", size=0.2) +
   ylab("Number of Modules") +
   scale_fill_manual(values=c(Unchanged_col, Decreased_col, Increased_col,
                              Opposite_col, NotDE_col),
                     name="Addition of Second Stress",
                     breaks=c("Unchanged", "Decreased_mag",
                              "Increased_mag", "Opposite_dir", "Not_DE"),
                     labels=c("Un-changed",
                              "Decreased Magnitude", "Increased Magnitude",
                              "Opposite Direction", "No Longer Differentially Expressed")) +
   theme_bw(base_size = 14)

ggsave("/Users/emilydittmar/Google Drive/Active Projects/Transcriptomics_Exp/Manuscript/SaltxNut/Figures/Modsbarplot.png")
ggsave("/Users/eld72413/Google Drive/Active Projects/Transcriptomics_Exp/Manuscript/SaltxNut/Figures/Modsbarplot.png")




######## Plot raw data
AllData_ordered <- AllData[order(AllData$Treatment, AllData$Accession, decreasing = FALSE),]
barplot(AllData_ordered$MEwhite, names.arg = AllData_ordered$Treatment)
barplot(AllData_ordered$MEgreenyellow, names.arg = AllData_ordered$Treatment)


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


