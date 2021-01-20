# RGR, AG Biomass, BG Biomass, chlorophyll SLA

#########################
######## SETUP ##########
#########################

library(dplyr)
library(ggplot2)
source("Functions.R")
library(car)
library(emmeans)

pheno <- read.csv("DataFiles/StudyDesign_Inbred_noOut.csv", header = T)
initial <- read.csv("DataFiles/InitialBiomass_070618.csv", header = T)
sla <- read.csv("DataFiles/SLA_data.csv", header = T)


#########################
####### CALC SLA ########
#########################

# Leaf Area / Dry Biomass

length(pheno$Plant) #104
length(sla$plantID) #350
length(unique(sla$plantID)) #336

# combine Area data
sla_area <- aggregate(sla$Area, by=list(sla$plantID), sum)
length(sla_area$Group.1) #336
colnames(sla_area) <- c("Plant", "LeafArea")

# which samples don't have leaf area?
pheno[which(!pheno$Plant %in% sla_area$Plant),] # plant 72

# combine SLA area data
pheno_sla <- merge(pheno, sla_area, by = "Plant", all.x = TRUE)
length(pheno_sla$Plant) #104

pheno_sla$SLA <- pheno_sla$LeafArea / pheno_sla$SLA.biomass
hist(pheno_sla$SLA)
length(which(is.na(pheno_sla$SLA))) #1

pheno_sla$Treatment <- factor(pheno_sla$Treatment, levels = c("Control", "LowNut",
                                                        "HighSalt", "Combo"))

#########################
##### SLA LS MEANS ######
#########################

modSLA_Acc <- lm(SLA ~ Treatment + Accession +
                   Treatment:Accession, data=pheno_sla)
SLA_meansAcc <- emmeans(modSLA_Acc, ~ Accession | Treatment, type = "response")
pairs(SLA_meansAcc) #only marginally significant: HA 445-RHA 274 (Low Nutrient)

# no pairwise differences among accessions, so can look at heterotic groups?

modSLA <- lm(SLA ~ Treatment + Group + Group:Accession +
               Treatment:Group +
               Treatment:Group:Accession, data=pheno_sla)
drop1(modSLA, test = "F") # treatment x group x accession is ns
Anova(modSLA, type = 2) # no treatment x group interactions significant

SLA_meansGroup <- emmeans(modSLA, ~ Group | Treatment, type = "response")
pairs(SLA_meansGroup) # different in low nutrient

SLA_means <- emmeans(modSLA, ~ Treatment | Group, type = "response")
pairwisediffs <- pairs(SLA_means) 
# combo-control different for HA; 
# combo-control, combo-LowNut, control-HighSalt, HighSalt-LowNut for RHA

SLA_df <- as.data.frame(SLA_means)
colnames(SLA_df)[3] <- c("SLA_mean")

#########################
####### SLA GRAPH #######
#########################

barplot(SLA_df[order(SLA_df$Treatment, SLA_df$Group)]$SLA_mean, 
        names.arg = SLA_df$Treatment)


SLA_plot <- ggplot(SLA_df, aes(y=SLA_mean, x=Treatment, color=Group))
SLA_plot + geom_col(position = position_dodge(), color="black", size=0.4) +  
  coord_cartesian(ylim=c(275, 340))

SLA_plot <- ggplot(SLA_df, aes(y=SLA_mean, x=Treatment, color=Group))
SLA_plot + geom_errorbar(aes(ymin=SLA_mean-SE, ymax=SLA_mean+SE), width=.1,
                         position=position_dodge(0.1)) +
  #geom_line(position=position_dodge(0.1)) +
  geom_point(position = position_dodge(0.1)) +
  theme_bw(base_size = 14)

# plot all data

SLA_plot_all <- ggplot(pheno_sla, aes(y=SLA, x=Treatment, fill=Group))
SLA_plot_all + geom_boxplot(outlier.shape = NA) +
  ylim(225,375) +
  theme_bw(base_size = 14)

#########################
### CHLOROPHYLL GRAPH ###
#########################

chlor_plot_all <- ggplot(pheno_sla, aes(y=Chlorophyll, x=Treatment, fill=Group))
chlor_plot_all + geom_boxplot(outlier.shape = NA) +
  ylim(10,27.5) +
  theme_bw(base_size = 14)

# is chlorophyll related to SLA?
mod <- lm(Chlorophyll ~ SLA + Treatment +
            Group + Group:Accession +
            Treatment:Group, data=pheno_sla)
summary(mod) #n.s SLA slope

plot(pheno_sla$SLA, pheno_sla$Chlorophyll)
abline(mod)

#########################
##### AG:BG BIOMASS #####
#########################

pheno_sla$AG_BG <- pheno_sla$Tot_AG / pheno_sla$Tot_BG
hist(pheno_sla$AG_BG)

ag_bg_plot <- ggplot(pheno_sla, aes(y=AG_BG, x=Treatment, fill=Group))
ag_bg_plot + geom_boxplot(outlier.shape = NA) +
  ylim(0,12.5) +
  theme_bw(base_size = 12) +
  ylab("Above:Below Biomass")

# check to see how mold influences

#########################
####### CALC RGR ########
#########################

# lsmeans for initial biomass for each accession

# for now will use full data
full <- read.csv("DataFiles/fulldata.csv", header=T)
full_inbred <- subset(full, full$Plant. %in% pheno$Plant)
length(full_inbred$Plant.) #104

rgr_plot <- ggplot(full_inbred, aes(y=RGR, x=Treatment, fill=Accession))
rgr_plot + geom_boxplot()

#mold?
ggplot(full_inbred, aes(y=Belowground.biomass, x=Treatment, fill=mold)) + geom_boxplot()

#########################
####### MODELS ########
#########################

# ask whether any traits are significantly influenced by Nutrients, Salt and/or interaction