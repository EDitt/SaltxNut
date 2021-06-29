#### relationship between modules + phenotypes?

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

ModInfo <- read.csv("ResultsFiles/Coexpression/Module_EGs.csv", header=T)

# change numbers/characters to factors
ModInfo[,2:10] <- lapply(ModInfo[,2:10], factor) # doesn't work

ModInfo$Treatment <- relevel(as.factor(ModInfo$Treatment), ref="Control")
levels(ModInfo$Treatment)
levels(ModInfo$Accession)
levels(ModInfo$Bench)


################################
############ MODEL #############
################################

### ANCOVA to investigate relationship between modules and traits

Pheno_Mod <- function(Yvar, module, dataset) {
  mod <- lm(Yvar ~ module + Treatment +
              #module:Treatment +
              Accession +
              Treatment:Accession +
              Bench, data = dataset,
            contrast=list(Accession=contr.sum, Bench=contr.sum))
  return(mod)
}

##############################
######## CHLOROPHYLL #########
##############################

hist(ModInfo$Chlorophyll)
# chlorophyll is lower in nutrient, higher in salt, and intermediate with control & combo

Chlorophyll_mod <- lapply (ModInfo[,c(17:104)], function(x) {Pheno_Mod(ModInfo$Chlorophyll, x, ModInfo) })
str(Chlorophyll_mod)

Chlorophyll_mod_ANOVA <- lapply (Chlorophyll_mod, function(x) {as.data.frame(Anova(x, test="F", type=2))})

critp <- 0.05/length(ModInfo[,c(17:104)]) # 0.00057

Chlorophyll_mod_Sig <- lapply(Chlorophyll_mod_ANOVA, function(x) which(x["module", "Pr(>F)"] < critp))
length(which(Chlorophyll_mod_Sig==1)) #1
names(which(Chlorophyll_mod_Sig==1)) # magenta4

summary(Chlorophyll_mod$MEmagenta4) # slope=-12.3, p=0.006
plot(ModInfo$Chlorophyll ~ ModInfo$MEmagenta4)
abline(Chlorophyll_mod$MEmagenta4)
# low nutrient has lower chlorophyll - according to this relationship, lower chlorophyll=higher expression of this module, but nutrient has lower expression of this module
# !!! bargraph is showing higher expression in nutrient....???

Anova(Chlorophyll_mod$MEmagenta4) # treatment has a much higher F value
which(lapply(Chlorophyll_mod, function(x) summary(x)$coefficients[2,4] < critp)=="TRUE") # magenta4

summary(Chlorophyll_mod$MEmagenta4)$coefficients[2,4]
summary(Chlorophyll_mod$MEthistle1)$coefficients[2,4]

# at p=0.01:
Chlorophyll_mod_Sig2 <- lapply(Chlorophyll_mod_ANOVA, function(x) which(x["module", "Pr(>F)"] < 0.001))
length(which(Chlorophyll_mod_Sig2==1)) # 2 (12 at p<0.01)
names(which(Chlorophyll_mod_Sig2==1)) # thistle1 also

# thistle1 has significantly higher expression in the nutrient module
##### !!!! wait- after plotting raw data it looks more like expression in thistle model is reduced... need to check the ls means & barplots...
summary(Chlorophyll_mod$MEthistle1) # slope=-7.4, p=0.0975
plot(ModInfo$Chlorophyll ~ ModInfo$MEthistle1)

# interactions between treatment and module?
Chlorophyll_mod_Sig_Inter <- lapply(Chlorophyll_mod_ANOVA, function(x) which(x["module:Treatment", "Pr(>F)"] < 0.05))
length(which(Chlorophyll_mod_Sig_Inter==1)) # darkturquoise is only significant at p<0.05 - this is one that is highly correlated with accession
names(which(Chlorophyll_mod_Sig_Inter==1)) 

##############################
##### ROOT MASS FRACTION #####
##############################
ModInfo$RootMassFraction <- ModInfo$Tot_BG / ModInfo$Tot_Biomass
hist(ModInfo$RootMassFraction)

Root_mod <- lapply (ModInfo[,c(17:104)], function(x) {Pheno_Mod(ModInfo$RootMassFraction, x, ModInfo) })

Root_mod_ANOVA <- lapply (Root_mod, function(x) {as.data.frame(Anova(x, test="F", type=2))})

which(lapply(Root_mod_ANOVA, function(x) which(x["module", "Pr(>F)"] < 0.001))==1) #0 at critical p value, 1 at p<0.001 (MEbisque4)
# this module is strongly affected by accession- which makes sense considering this trait influenced by heterotic group

Root_mod_ANOVA$MEbisque4 # module has highest F value
drop1(Root_mod$MEbisque4, test = "F")




######
Chlorophyll_mod_Sig <- lapply(Chlorophyll_mod_ANOVA, function(x) which(x["module", "Pr(>F)"] < critp))
length(which(Chlorophyll_mod_Sig==1)) #1
names(which(Chlorophyll_mod_Sig==1)) # magenta4


##############################
### ANOVA
Anova(Chlorophyll_mod$MEbisque4)

Chlorophyll_mod_ANOVA$MEbisque4[1,4]
Chlorophyll_mod_ANOVA$MEbisque4["module", "Pr(>F)"]


### scratch

# p-values for module in summary:
summary(Chlorophyll_mod$MEbisque4)$coefficients[2,4]

# plot raw values:

ModInfo_ordered <- ModInfo[order(ModInfo$Treatment, ModInfo$Accession, decreasing = FALSE),]
ModInfo_ordered2 <- ModInfo[order(ModInfo$Accession, ModInfo$Treatment, decreasing = FALSE),]
barplot(ModInfo_ordered$MEbisque4, names.arg = ModInfo_ordered$Treatment)

