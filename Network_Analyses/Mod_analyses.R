
#########################
######## SETUP ##########
#########################

library(dplyr)
library(ggplot2)
library(ggthemes)
library(car)
source("Functions.R")
library(tidyr)
library(emmeans)
library(UpSetR)
library("viridis")
library(RColorBrewer)


#########################
###### READ DATA ########
#########################

ModEGs <- read.csv("ResultsFiles/Coexpression/deepSplit4/ModEigengenes_pamTRUE.csv", header=T)
colnames(ModEGs)[1] <- "Plant"

ModExp <- read.csv("ResultsFiles/Coexpression/deepSplit4/ModExpression_pamTRUE.csv", header=T)
colnames(ModExp)[1] <- "Plant"

Design <- read.csv("DataFiles/StudyDesign_Inbred_noOut.csv", header=T)

Design$Treatment <- relevel(Design$Treatment, ref="Control")

levels(Design$Treatment)
levels(Design$Accession)
levels(Design$Group)
Design$Cross <- factor(Design$Cross)
Design$SampleDay <- factor(Design$SampleDay)

## samples per accssion:
aggregate(Design$Plant, by=list(Design$Accession, Design$Treatment), length)

#########################
###### DATA SETUP #######
#########################

AllData <- merge(Design[,c(2:8, 14, 17, 18, 23, 33, 42, 43)], ModEGs, by="Plant")
str(AllData)

# convert to wide format:
AllData_long <- gather(AllData, Module, eigengeneExp, 
                       MEgreen4:MEthistle1, factor_key = TRUE)
head(AllData_long)

#########################
#### LINEAR MODEL #######
#########################


EigenExp <- lm(eigengeneExp ~ Treatment + Group + 
                 Group:Accession + 
                 Module + 
                 Treatment:Module +
                 Reproductive +
                 SampleDay, 
               data = AllData_long)
Anova(EigenExp) 
drop1(EigenExp, test = "F") # reproductive is n.s.
hist(resid(EigenExp))

EigenMean <- emmeans(EigenExp, ~ Treatment | Module, type = "response")
EigenMean
pairwisediffs <- pairs(EigenMean)
pairwisediffs_df <- as.data.frame(pairwisediffs)
length(pairwisediffs_df$contrast)

pairwisediffs_df$sig <- ifelse(pairwisediffs_df$p.value < 0.05, "Sig", "NS")
length(which (pairwisediffs_df$sig == "Sig")) #108 (208 before merging similar modules)
aggregate(pairwisediffs_df$contrast, 
          by=list(pairwisediffs_df$contrast, pairwisediffs_df$sig), length)

str(pairwisediffs_df)

### split into lists based on contrasts then use functions below:

Contrasts <- split(pairwisediffs_df, pairwisediffs_df$contrast)
Sig_Contrasts <- lapply(Contrasts, SigDEdf, PvaluesCol=7, CritP=0.05)

Mod_overlap <- GeneSets(Sig_Contrasts$`Control - Combo`$Module,
                        Sig_Contrasts$`Control - HighSalt`$Module,
                        Sig_Contrasts$`Control - LowNut`$Module)
lapply(Mod_overlap, function(x) {length(x)})
# In common all: 4 (no pairwise differences)
# DE Combo + DE High Salt: 6
# DE High Salt + DE Low Nut: 2 (not different)
# DE Combo + DE Low Nut: 5
# DE Combo only: 1
# DE High Salt only: 3
# DE Low Nut only: 19

Mod_overlap2 <- GeneSets(Sig_Contrasts$`Combo - HighSalt`$Module,
                        Sig_Contrasts$`Combo - LowNut`$Module,
                        Sig_Contrasts$`HighSalt - LowNut`$Module)
lapply(Mod_overlap2, function(x) {length(x)})
# In common all: 1
# Combo-High Salt & Combo-Low Nut: 0
# Combo-Low Nut & High Salt - Low Nut: 15
# Combo-High Salt & High Salt-Low Nut: 3
# Combo-High Salt only: 1
# Combo-Low Nut only: 0
# High Salt-Low Nut only: 7

#########################
#### MULTIVARIATE #######
#########################
# https://data.library.virginia.edu/getting-started-with-multivariate-multiple-regression/
#mods <- list(colnames(AllData[,15:100]))

mlm1 <- lm(do.call(cbind, AllData[,c(15:100)]) ~ Treatment + Group + 
             Group:Accession + 
             Reproductive +
             SampleDay, data = AllData)
summary(mlm1)
Anova(mlm1)
head(vcov(mlm1))

#########################
###### UPSET PLOT #######
#########################

# all the modules that are different from control
DE_Mods <- union(Sig_Contrasts$`Control - Combo`$Module, 
                 union(Sig_Contrasts$`Control - HighSalt`$Module, Sig_Contrasts$`Control - LowNut`$Module))

# only take the modules that intersect with the above (don't care about pairwise contrasts if not diff from control)
Contrasts_DE <- lapply(Sig_Contrasts, function(x) {intersect(x$Module, DE_Mods)})

upset(fromList(Contrasts_DE),
      #keep.order = TRUE,
      order.by = "freq",
      #group.by = "sets", 
      nsets = 13,
      #empty.intersections = "on",
      nintersects = 30)

#########################
###### CATEGORIES ######
#########################

Sig_pairwise <- GeneSets(Contrasts_Graph$`Combo - HighSalt`,
                         Contrasts_Graph$`Combo - LowNut`,
                         Contrasts_Graph$`HighSalt - LowNut`)


lapply(Sig_pairwise, function(x) {length(x)})

##### Pairwise Diffs that are Significant:
# In common all: 9
# Combo-Salt + Combo-Nut: 0
# Combo-Nut + Salt-Nut: 18
# Combo-Salt + Salt-Nut: 6
# Combo-Salt only: 1
# Combo-Nut only: 0
# Salt-Nut only: 10

EigenMean_df <- as.data.frame(EigenMean)

# significant modules
Sig_Mods_df <- subset(EigenMean_df, Module %in% DE_Mods)

## wide based on treatment

# convert to wide based on treatment
Sig_Mods_df_treatwide <- spread(Sig_Mods_df[,1:3], Treatment, emmean)
head(Sig_Mods_df_treatwide)

# All shared- shared by all treatments, not different
Pairwise_diffs <- union(Sig_Contrasts$`Combo - HighSalt`$Module, 
                        union(Sig_Contrasts$`Combo - LowNut`$Module,
                              Sig_Contrasts$`HighSalt - LowNut`$Module)) # N=27

All_shared <- setdiff(Mod_overlap$InCommonAll, Pairwise_diffs) #4 of 4 modules shared are not different from one another
#All_shared_NutDiff <- intersect(Mod_overlap$InCommonAll, Pairwise_diffs) #0

# Unconditional nutrient/salt - shared by nutrient/salt + combo, not different
Uncond_Nut <- setdiff(Mod_overlap$"Sig_Contrasts$`Control - Combo`$ModuleSig_Contrasts$`Control - LowNut`$ModuleOnly",
                        Sig_Contrasts$`Combo - LowNut`$Module) # 2 of the 5 Nut-Combo Mods are not different

Uncond_Salt <- setdiff(Mod_overlap$"Sig_Contrasts$`Control - Combo`$ModuleSig_Contrasts$`Control - HighSalt`$ModuleOnly",
                      Sig_Contrasts$`Combo - HighSalt`$Module) # 6 of the 6 Salt-Combo Mods are not different

# shared by nutrients + combo, different
Cond_Nut <- intersect(Mod_overlap$"Sig_Contrasts$`Control - Combo`$ModuleSig_Contrasts$`Control - LowNut`$ModuleOnly",
                      Sig_Contrasts$`Combo - LowNut`$Module) # N=3

# only significant in nutrient
Nut_only <- Mod_overlap$"Sig_Contrasts$`Control - LowNut`$ModuleOnly" # N=19
# only significant in salt
Salt_only <- Mod_overlap$"Sig_Contrasts$`Control - HighSalt`$ModuleOnly" #N=3
# only significant in combo
Combo_only <- Mod_overlap$"Sig_Contrasts$`Control - Combo`$ModuleOnly" #N=1

Nut_Salt_same <- setdiff(Mod_overlap$"Sig_Contrasts$`Control - HighSalt`$ModuleSig_Contrasts$`Control - LowNut`$ModuleOnly",
                         Sig_Contrasts$`HighSalt - LowNut`$Module) #2

Nut_Salt_diff <- intersect(Mod_overlap$"Sig_Contrasts$`Control - HighSalt`$ModuleSig_Contrasts$`Control - LowNut`$ModuleOnly",
                         Sig_Contrasts$`HighSalt - LowNut`$Module) #0

Sig_Mods_df_treatwide$cat <- ifelse(Sig_Mods_df_treatwide$Module %in% All_shared, "All_shared",
                                    ifelse(Sig_Mods_df_treatwide$Module %in% All_shared_NutDiff, "All_shared_NutDiff", 
                                           ifelse(Sig_Mods_df_treatwide$Module %in% Uncond_Nut, "Unconditional_Nut",
                                                  ifelse(Sig_Mods_df_treatwide$Module %in% Uncond_Salt, "Unconditional_Salt",
                                                         ifelse(Sig_Mods_df_treatwide$Module %in% Cond_Nut, "Conditional_Nut",
                                                                ifelse(Sig_Mods_df_treatwide$Module %in% Nut_only, "Nutrient_only",
                                                                       ifelse(Sig_Mods_df_treatwide$Module %in% Salt_only, "Salt_only",
                                                                              ifelse(Sig_Mods_df_treatwide$Module %in% Combo_only, "Combo_only",
                                                                                     ifelse(Sig_Mods_df_treatwide$Module %in% Nut_Salt_same, "Nut-Salt_shared",
                                                                                          "Nut-Salt_diff")))))))))
# check
aggregate(Sig_Mods_df_treatwide$Module, by=list(Sig_Mods_df_treatwide$cat), length)


#########################
#### REGRESSION PLOT ####
#########################

head(Sig_Mods_df_treatwide)
plot(Sig_Mods_df_treatwide$HighSalt ~ Sig_Mods_df_treatwide$LowNut)

pmod <- ggplot(data=Sig_Mods_df_treatwide, aes(LowNut, HighSalt))
pmod + geom_point()
pmod + geom_point(aes(fill=cat, color = cat), size = 3, stroke = 1)

pmod + geom_point(aes(fill=cat, shape=cat, color = cat), size = 3, stroke = 1) +
  scale_shape_manual(values=c(19, 19, 24,
                              19, 19, 19,
                              24, 24,
                              19, 19))

# Lots of nutrient (or conditional nutrient) only that appear in opposite direction in salt


#########################
#### PCA EXPRESSION #####
#########################

head(ModExp)

rownames(ModExp) <- ModExp$Plant
head(ModExp[,-1])
PCA <- prcomp(ModExp[,-1])
PCA_df <- as.data.frame(PCA$x)
PCA_df$Plant <- rownames(PCA_df)

plot(PCA_df$PC2 ~ PCA_df$PC1)

AllExpPCA <- merge(Design[,c(2, 8, 14, 17)], PCA_df[,c(1,2,103)], by="Plant")

p1 <- ggplot(data=AllExpPCA, aes(x=PC1, y=PC2))

p1 + geom_point(aes(color = Treatment, shape = Group))


#########################
#### PLOT EXPRESSION ####
#########################

head(DE_Mods)
AllExp <- merge(Design[,c(2, 8, 14, 17)], ModExp, by="Plant")
str(AllExp)



#####################################
#SCRATCH BELOW

#########################
#### LINEAR MODELS ######
#########################

lm( Sepal.Width ~ C(Species,contr.treatment(3, base=2)), data=iris )

lm_FUN <- function(Yvar, dataset) {
  mod1 <- lm(Yvar ~ Group + Group:Cross + SampleDay + Treatment, data = dataset)
  pvals <- Anova(mod1, type = "II") 
  return(pvals[3,4])
}

result <- lapply (AllData[,c(15:85)], function(x) {lm_FUN(x, AllData) })
str(result)
critp <- 0.05/length(result)
length(which(result < critp)) # N=46
Inter <- names(which(result < critp))

# Inter is a list of modules that significantly affect treatment



##########################
#### LIKELIHOOD RATIO ####
##########################

### make treatment 2 factors with 2 levels each

AllData$Nut <- as.factor(ifelse(AllData$Treatment == "LowNut" | AllData$Treatment == "Combo", "1", "0"))
AllData$Salt <- as.factor(ifelse(AllData$Treatment == "HighSalt" | AllData$Treatment == "Combo", "1", "0"))

#check
aggregate(AllData$Plant, by=list(AllData$Nut, AllData$Salt, AllData$Treatment), length)

lm_FUN2 <- function(Yvar, dataset) {
  mod1 <- lm(Yvar ~ Nut + Salt + Nut:Salt +
               Group + Group:Cross + SampleDay, data = dataset)
  result <- drop1(mod1, test = "Chi")
  return(result[3,5])
}

result2 <- lapply (AllData[,c(15:85)], function(x) {lm_FUN2(x, AllData) })
str(result2)
length(which(result2 < critp)) # N=21
which(result2 < critp)
NutxSalt <- names(which(result2 < critp))  ### NutxSalt modules

### main effects
# take out Nut:Salt interaction

lm_FUN3 <- function(Yvar, dataset) {
  mod1 <- lm(Yvar ~ Nut + Salt + 
               Group + Group:Cross + SampleDay, data = dataset)
  result <- drop1(mod1, test = "Chi")
  main_result <- data.frame("LowNutSig" = result[2,5], "HighSaltSig" = result[3,5])
  return(main_result)
}

result3 <- lapply (AllData[,c(15:85)], function(x) {lm_FUN3(x, AllData) })
str(result3)

length(which(result3$LowNutSig < critp)) #0
length(which(result3$HighSaltSig < critp)) #0

Nut <- lapply (result3, function(x) {which(x$LowNutSig < critp) })
length(which(Nut==1)) #40

Salt <- lapply (result3, function(x) {which(x$HighSaltSig < critp) })
length(which(Salt==1)) #19

SaltMods <- names(which(Salt==1))
NutMods <- names(which(Nut==1))

both <- intersect(names(which(Salt==1)), names(which(Nut==1)))
length(both) #N = 11
Main_Inter <- intersect(Inter, both)
length(Main_Inter) #N=11

### Nutrient only
Nut_only <- setdiff(names(which(Nut==1)), both)
length(Nut_only) #29

### Salt only
Salt_only <- setdiff(names(which(Salt==1)), both)
length(Salt_only) #8




