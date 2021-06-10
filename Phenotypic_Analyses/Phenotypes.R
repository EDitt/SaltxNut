# RGR, AG Biomass, BG Biomass, chlorophyll SLA

#########################
######## SETUP ##########
#########################

library(tidyr)
library(dplyr)
library(ggplot2)
source("Functions.R")
library(car)
library(emmeans)

inbred_lines <- c("HA 412 HO", "HA 445", "HA 89",
                  "RHA 274", "RHA 373", "RHA 377",
                  "SAM 206", "SAM 227")

############################
# CALCULATE INITIAL WEIGHT #
############################

initial <- read.csv("DataFiles/InitialBiomass_070618.csv", header = T)

initial_inbred <- subset(initial, Accession %in%
                           inbred_lines)
initial_inbred$Accession <- droplevels(initial_inbred$Accession)
levels(initial_inbred$Accession)

initial_inbred$Tray <- factor(initial_inbred$Tray)
initial_inbred$Plate <- factor(initial_inbred$Plate)
str(initial_inbred)

# determine factors affecting initial weight
# transplant day is correlated with "Tray" so will use only "Tray" in model
Mod1 <- lm(Tot_biomass ~ Accession + Plate + Tray,
           data=initial_inbred)
drop1(Mod1, test = "F") # only accession is significant
summary(Mod1)

Means <- emmeans(Mod1, ~Accession)
pairs(Means)

# how different from raw means?
aggregate(initial_inbred$Tot_biomass, by=list(initial_inbred$Accession), 
          mean, na.rm=T)
# very small differences

Mod_AG <- lm(AboveGround_Biomass ~ Accession + Plate + Tray,
           data=initial_inbred)
drop1(Mod_AG, test = "F") # only accession is significant

Mod_BG <- lm(BelowGround_Biomass ~ Accession + Plate + Tray,
             data=initial_inbred)
drop1(Mod_BG, test = "F") # accession is significant, tray marginally significant (p=0.08)

Means_AG <- emmeans(Mod_AG, ~Accession)
Means_BG <- emmeans(Mod_BG, ~Accession)

Initial_df <- merge(Means, 
                    merge(Means_AG, Means_BG, by="Accession"), 
                    by="Accession")
Initial_df <- Initial_df[,-c(4:6,9:11,14:16)]
colnames(Initial_df) <- c("Accession", "Tot_Biomass", "Tot_Biomass_SE",
                          "AG_Biomass", "AG_Biomass_SE", "BG_Biomass",
                          "BG_Biomass_SE")

# initial ratio of AG / BG
Mod_AGBG_init <- lm(AG.BG ~ Accession + Plate + Tray,
             data=initial_inbred)
drop1(Mod_AGBG_init, test = "F")
emmeans(Mod_AGBG_init, ~Accession)

############################
###### SET UP FULL DF ######
############################

full_data <- read.csv("DataFiles/fulldata.csv", header=T)
full_inbred <- subset(full_data, Accession %in% inbred_lines)
length(full_inbred$Plant.)

full_inbred$Treatment <- factor(full_inbred$Treatment, levels = c("Control", "LowNut",
                                                                  "HighSalt", "Combo"))
full_inbred$Bench <- factor(full_inbred$Bench)

# add initial biomass
full_inbred$Initial_Biomass <- Initial_df[match(full_inbred$Accession, Initial_df$Accession),
                                          "Tot_Biomass"]
# add initial aboveground biomass
full_inbred$Initial_AGBiomass <- Initial_df[match(full_inbred$Accession, Initial_df$Accession),
                                                "AG_Biomass"]
# add initial belowground biomass
full_inbred$Initial_BGBiomass <- Initial_df[match(full_inbred$Accession, Initial_df$Accession),
                                                "BG_Biomass"]

############################
####### CALCULATIONS #######
############################

# RGR = LN(final biomass - initial biomass) / # days
full_inbred$RGR_new <- (log(full_inbred$Tot_Biomass) -
                          log(full_inbred$Initial_Biomass)) /
                        full_inbred$growdays

# Aboveground & Belowground biomass acquired (minus initial)
full_inbred$Acquired_AG <- full_inbred$Tot_AG - full_inbred$Initial_AGBiomass
full_inbred$Acquired_BG <- full_inbred$Tot_BG - full_inbred$Initial_BGBiomass
hist(full_inbred$Acquired_AG)
hist(full_inbred$Acquired_BG)

# Aboveground/Belowground biomass
full_inbred$AG_BG <- full_inbred$Tot_AG / full_inbred$Tot_BG

############################
########## MODEL ###########
############################

# Factors of interest: HA vs. RHA, Accession (nested within heterotic group), treatment
# Nuisance variables: Bench (which greenhouse bench plant was located on),
#                     Reproductive (whether or not plant had begun flowering at time of sampling),
#                     growdays (number of days between transplant and harvest)
#                     mold (some tissue didn't dry all the way and got moldy before getting weighed)

PhenotypeMod <- function(dataframe, yvar) {
  Mod <- lm(yvar ~ Hybrid_role + Hybrid_role:Accession +
              Hybrid_role:Treatment +
              Treatment + Bench + Reproductive + growdays + mold,
            data=dataframe)
  return (Mod)
}

############################
# TAKE OUT LINES NOT USED ##
############################

# take out SAM lines (not included in RNAseq exp)
full_inbred_sub <- subset(full_inbred, Hybrid_role!="")
full_inbred_sub$Hybrid_role <- droplevels(full_inbred_sub$Hybrid_role)
levels(full_inbred_sub$Hybrid_role)
aggregate(full_inbred_sub$Plant., by=list(full_inbred_sub$Hybrid_role,
                                          full_inbred_sub$Accession), length) # check

############################
######## RGR MEANS #########
############################
# factors that affect RGR
Mod_RGR <- PhenotypeMod(full_inbred_sub, full_inbred_sub$RGR_new)
drop1(Mod_RGR, test = "F")
Anova(Mod_RGR, test = "F")
summary(Mod_RGR)
hist(resid(Mod_RGR))
plot(Mod_RGR)

Means_RGR <- emmeans(Mod_RGR, ~Accession:Treatment) # by accession
Means_RGR <- emmeans(Mod_RGR, ~Hybrid_role|Treatment)
pairs(Means_RGR) # only significant in Low nutrient (p=0.0157)

# plot raw data
rgr_plot_raw <- ggplot(full_inbred_sub, aes(y=RGR_new, x=Treatment, fill=Hybrid_role))
rgr_plot_raw + geom_boxplot()

# plot lsmeans
rgr_plot <- ggplot(as.data.frame(Means_RGR), aes(y=emmean, x=Treatment, group=Hybrid_role))
rgr_plot + geom_bar(stat = "identity", position = position_dodge())
rgr_plot + geom_point(stat = "identity", aes(color=Hybrid_role, shape=Hybrid_role),
                      position=position_dodge(width = 0.2), size=4) +
  geom_line(stat = "identity", aes(color=Hybrid_role, linetype=Hybrid_role),
            position=position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE, color=Hybrid_role), width=.2,
                position=position_dodge(width = 0.2)) +
  scale_colour_discrete(name="", breaks=c("Female", "Male"), labels=c("HA", "RHA")) +
  scale_shape_discrete(name="", breaks=c("Female", "Male"), labels=c("HA", "RHA")) +
  scale_linetype_discrete(name="", breaks=c("Female", "Male"), labels=c("HA", "RHA")) +
  guides(fill=guide_legend(title=NULL)) +
  theme_minimal() + ylab("Relative Growth Rate")
ggsave("/Users/emilydittmar/Google Drive/Active Projects/Transcriptomics_Exp/Analyses/Figures/RGR.png")

############################
######## AG BIOMASS ########
############################

Mod_AG <- PhenotypeMod(full_inbred_sub, full_inbred_sub$Acquired_AG)
drop1(Mod_AG, test = "F") # mold, hybrid role x accession, hybrid role x treatment
Anova(Mod_AG, test = "F") # hybrid role, Treatment also
summary(Mod_AG)
hist(resid(Mod_AG))
plot(Mod_AG)

Means_AG <- emmeans(Mod_AG, ~Hybrid_role|Treatment)
pairs(Means_AG) # only sig. in control (p < 0.0001)

# plot raw data
#ag_plot_raw <- ggplot(full_inbred_sub, aes(y=Acquired_AG, x=Treatment, fill=Hybrid_role))
ag_plot_raw <- ggplot(full_inbred_sub, aes(y=Acquired_AG, x=Treatment, fill=Accession))
ag_plot_raw + geom_boxplot()

# plot lsmeans
ag_plot <- ggplot(as.data.frame(Means_AG), aes(y=emmean, x=Treatment, group=Hybrid_role))
ag_plot + geom_point(stat = "identity", aes(color=Hybrid_role, shape=Hybrid_role),
                      position=position_dodge(width = 0.2), size=4) +
  geom_line(stat = "identity", aes(color=Hybrid_role, linetype=Hybrid_role),
            position=position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE, color=Hybrid_role), width=.2,
                position=position_dodge(width = 0.2))

# do accessions differ?
Means_AG2 <- emmeans(Mod_AG, ~Accession:Treatment) # by accession
ag_plot2 <- ggplot(as.data.frame(Means_AG2), aes(y=emmean, x=Treatment, group=Accession))
ag_plot2 + geom_point(stat = "identity", aes(color=Accession, shape=Hybrid_role),
                     position=position_dodge(width = 0.2), size=4) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE, color=Accession), width=.2,
                position=position_dodge(width = 0.2))

############################
######## BG BIOMASS ########
############################

Mod_BG <- PhenotypeMod(full_inbred_sub, full_inbred_sub$Acquired_BG)
drop1(Mod_BG, test = "F") # bench (p=0.07), growdays (p=0.07), hybrid role x accession
Anova(Mod_BG, test = "F") # hybrid role, Treatment also
summary(Mod_BG)
hist(resid(Mod_BG))
plot(Mod_BG)

Means_BG <- emmeans(Mod_BG, ~Hybrid_role|Treatment)
pairs(Means_BG) # only sig. different in control (p=0.0010)

# plot raw data
bg_plot_raw <- ggplot(full_inbred_sub, aes(y=Acquired_BG, x=Treatment, fill=Hybrid_role))
bg_plot_raw <- ggplot(full_inbred_sub, aes(y=Acquired_BG, x=Treatment, fill=Accession))
bg_plot_raw + geom_boxplot()

# plot lsmeans
bg_plot <- ggplot(as.data.frame(Means_BG), aes(y=emmean, x=Treatment, group=Hybrid_role))
bg_plot + geom_point(stat = "identity", aes(color=Hybrid_role, shape=Hybrid_role),
                     position=position_dodge(width = 0.2), size=4) +
  geom_line(stat = "identity", aes(color=Hybrid_role, linetype=Hybrid_role),
            position=position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE, color=Hybrid_role), width=.2,
                position=position_dodge(width = 0.2))

# do accessions differ?
Means_BG2 <- emmeans(Mod_BG, ~Accession:Treatment) # by accession
bg_plot2 <- ggplot(as.data.frame(Means_BG2), aes(y=emmean, x=Treatment, group=Accession))
bg_plot2 + geom_point(stat = "identity", aes(color=Accession, shape=Hybrid_role),
                      position=position_dodge(width = 0.2), size=4) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE, color=Accession), width=.2,
                position=position_dodge(width = 0.2))

############################
#### AG/BG BIOMASS PLOT ####
############################

# how correlated?
plot(full_inbred_sub$Acquired_AG ~ full_inbred_sub$Acquired_BG)

# plot raw
AG_BG_raw <- ggplot(full_inbred_sub, aes(y=Acquired_AG/Acquired_BG, x=Treatment, fill=Hybrid_role))
AG_BG_raw + geom_boxplot()

Means_AG <- as.data.frame(Means_AG)
Means_AG$Tissue <- "Aboveground"

Means_BG <- as.data.frame(Means_BG)
Means_BG$Tissue <- "Belowground"

AG_BG_Means <- rbind(Means_AG, Means_BG)

# to graph heterotic groups separately
AG_BG_Means$Treat_Het <- paste0(AG_BG_Means$Treatment, "_", AG_BG_Means$Hybrid_role) 
# accessions?
#AG_BG_Means$Treat_Het <- paste0(AG_BG_Means$Treatment, "_", AG_BG_Means$Accession) 

agbg_plot <- ggplot(AG_BG_Means, aes(y=emmean, x=Treatment, fill=Tissue))
agbg_plot + geom_bar(stat="identity", position=position_dodge())

agbg_plot <- ggplot(AG_BG_Means, aes(y=emmean, x=Treat_Het, fill=Tissue))
agbg_plot + geom_bar(stat="identity", position=position_dodge())
agbg_plot + geom_bar(stat="identity", position="fill") #to create proportion


##### wide format (scratch)
AG_BG_Means <- merge(Means_AG, Means_BG, 
                     by=c("Hybrid_role", "Treatment"))
colnames(AG_BG_Means)[3] <- c("Mean_AG")
colnames(AG_BG_Means)[8] <- c("Mean_BG")

AG_BG_Means <- AG_BG_Means[,-c(5:7,10:12)]

AG_BG_long <- gather(AG_BG_Means, Measurement, Value, Mean_AG, Mean_BG, factor_key=TRUE)

############################
## AG/BG BIOMASS LSMEANS ###
############################
Mod_AG.BG <- PhenotypeMod(full_inbred_sub, full_inbred_sub$AG_BG)
drop1(Mod_AG.BG, test = "F") # bench, mold, hybrid role x accession, hybrid role x treatment (p=0.1)
Anova(Mod_AG.BG, test = "F") # hybrid role, Treatment also
summary(Mod_AG.BG)
hist(resid(Mod_AG.BG))
plot(Mod_AG.BG)

Means_AG.BG <- emmeans(Mod_AG.BG, ~Hybrid_role|Treatment)
Means_AG.BG <- emmeans(Mod_AG.BG, ~Treatment|Hybrid_role)
pairs(Means_AG.BG) # only sig. different in combo (p=0.0465)

# plot raw data
ag.bg_plot_raw <- ggplot(full_inbred_sub, aes(y=AG_BG, x=Treatment, fill=Hybrid_role))
ag.bg_plot_raw <- ggplot(full_inbred_sub, aes(y=AG_BG, x=Treatment, fill=Accession))
ag.bg_plot_raw + geom_boxplot()

# plot lsmeans
ag.bg_plot <- ggplot(as.data.frame(Means_AG.BG), aes(y=emmean, x=Treatment, group=Hybrid_role))
ag.bg_plot + geom_point(stat = "identity", aes(color=Hybrid_role, shape=Hybrid_role),
                      position=position_dodge(width = 0.2), size=4) +
  geom_line(stat = "identity", aes(color=Hybrid_role, linetype=Hybrid_role),
            position=position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE, color=Hybrid_role), width=.2,
                position=position_dodge(width = 0.2)) +
  scale_colour_discrete(name="", breaks=c("Female", "Male"), labels=c("HA", "RHA")) +
  scale_shape_discrete(name="", breaks=c("Female", "Male"), labels=c("HA", "RHA")) +
  scale_linetype_discrete(name="", breaks=c("Female", "Male"), labels=c("HA", "RHA")) +
  theme_minimal() + ylab("Aboveground / Belowground Biomass")
ggsave("/Users/emilydittmar/Google Drive/Active Projects/Transcriptomics_Exp/Analyses/Figures/AG_BG.png")

# add line for AG/BG for "initial"? # actually initial AG/BG is out of range

# look at effect of treatment only?
Means_AG.BG2 <- emmeans(Mod_AG.BG, ~Treatment)
ag.bg_plot2 <- ggplot(as.data.frame(Means_AG.BG2), aes(y=emmean, x=Treatment))
ag.bg_plot2 + geom_point(stat = "identity", size=4) +
  geom_line(stat = "identity") +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.2)

############################
######## CHLOROPHYLL #######
############################

Mod_Ch <- PhenotypeMod(full_inbred_sub, full_inbred_sub$Chlorophyll)
drop1(Mod_Ch, test = "F") # growdays, hybrid role x accession
Anova(Mod_Ch, test = "F") # Treatment also
summary(Mod_Ch)
hist(resid(Mod_Ch))
plot(Mod_Ch)

Means_Ch <- emmeans(Mod_Ch, ~Hybrid_role|Treatment)
pairs(Means_Ch) # no sig. differences (in LowNut p-value=0.08)

# plot raw data
ch_plot_raw <- ggplot(full_inbred_sub, aes(y=Chlorophyll, x=Treatment, fill=Hybrid_role))
ch_plot_raw <- ggplot(full_inbred_sub, aes(y=Chlorophyll, x=Treatment, fill=Accession))
ch_plot_raw <- ggplot(full_inbred_sub, aes(y=Chlorophyll, x=Treatment))
ch_plot_raw + geom_boxplot(outlier.shape = NA)


# plot lsmeans
ch_plot <- ggplot(as.data.frame(Means_Ch), aes(y=emmean, x=Treatment, group=Hybrid_role))
ch_plot + geom_point(stat = "identity", aes(color=Hybrid_role, shape=Hybrid_role),
                     position=position_dodge(width = 0.2), size=4) +
  geom_line(stat = "identity", aes(color=Hybrid_role, linetype=Hybrid_role),
            position=position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE, color=Hybrid_role), width=.2,
                position=position_dodge(width = 0.2)) +
  scale_colour_discrete(name="", breaks=c("Female", "Male"), labels=c("HA", "RHA")) +
  scale_shape_discrete(name="", breaks=c("Female", "Male"), labels=c("HA", "RHA")) +
  scale_linetype_discrete(name="", breaks=c("Female", "Male"), labels=c("HA", "RHA")) +
  guides(fill=guide_legend(title=NULL)) +
  theme_minimal() + ylab("Chlorophyll")
ggsave("/Users/emilydittmar/Google Drive/Active Projects/Transcriptomics_Exp/Analyses/Figures/Chlorophyll.png")

#########################
####### CALC SLA ########
#########################

sla <- read.csv("DataFiles/SLA_data.csv", header = T)

# Leaf Area / Dry Biomass

length(full_inbred$Plant.) #155
length(sla$plantID) #350
length(unique(sla$plantID)) #336

# combine Area data
sla_area <- aggregate(sla$Area, by=list(sla$plantID), sum)
length(sla_area$Group.1) #336
colnames(sla_area) <- c("Plant.", "LeafArea")

# which samples don't have leaf area?
pheno[which(!full_inbred$Plant. %in% sla_area$Plant.),] # 57, 247 (RHA 274, RHA 373)

# combine SLA area data
pheno_sla <- merge(full_inbred, sla_area, by = "Plant.", all.x = TRUE)
length(pheno_sla$Plant.) #155

pheno_sla$SLA <- pheno_sla$LeafArea / pheno_sla$SLA.biomass
hist(pheno_sla$SLA)
length(which(is.na(pheno_sla$SLA))) #2
pheno_sla[which(pheno_sla$SLA > 1500),] # SAM 227 - likely a typo

# take out SAM lines (not included in RNAseq exp)
pheno_sla_sub <- subset(pheno_sla, Hybrid_role!="")
pheno_sla_sub$Hybrid_role <- droplevels(pheno_sla_sub$Hybrid_role)
levels(pheno_sla_sub$Hybrid_role)
aggregate(pheno_sla_sub$Plant., by=list(pheno_sla_sub$Hybrid_role,
                                        pheno_sla_sub$Accession), length) # check


# is chlorophyll related to SLA?
mod <- lm(Chlorophyll ~ SLA + Treatment +
            Hybrid_role + Hybrid_role:Accession +
            Treatment:Hybrid_role, data=pheno_sla_sub)
summary(mod) # no- n.s SLA slope

plot(pheno_sla_sub$SLA, pheno_sla_sub$Chlorophyll)
abline(mod)


# plot raw data
sla_plot_raw <- ggplot(pheno_sla_sub, aes(y=SLA, x=Treatment, fill=Hybrid_role))
sla_plot_raw <- ggplot(pheno_sla_sub, aes(y=SLA, x=Treatment, fill=Accession))
sla_plot_raw <- ggplot(pheno_sla_sub, aes(y=SLA, x=Treatment))
sla_plot_raw + geom_boxplot()

#ls means
Mod_sla <- PhenotypeMod(pheno_sla_sub, pheno_sla_sub$SLA)
drop1(Mod_sla, test = "F") # bench, hybrid role x accession
Anova(Mod_sla, test = "F") # Treatment also
summary(Mod_sla)
hist(resid(Mod_sla))
plot(Mod_sla)

Means_sla <- emmeans(Mod_sla, ~Hybrid_role|Treatment)
#Means_sla <- emmeans(Mod_sla, ~Treatment|Hybrid_role)
pairs(Means_sla) # not sig. different

# plot lsmeans
sla_plot <- ggplot(as.data.frame(Means_sla), aes(y=emmean, x=Treatment, group=Hybrid_role))
sla_plot + geom_point(stat = "identity", aes(color=Hybrid_role, shape=Hybrid_role),
                     position=position_dodge(width = 0.2), size=4) +
  geom_line(stat = "identity", aes(color=Hybrid_role, linetype=Hybrid_role),
            position=position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE, color=Hybrid_role), width=.2,
                position=position_dodge(width = 0.2)) +
  scale_colour_discrete(name="", breaks=c("Female", "Male"), labels=c("HA", "RHA")) +
  scale_shape_discrete(name="", breaks=c("Female", "Male"), labels=c("HA", "RHA")) +
  scale_linetype_discrete(name="", breaks=c("Female", "Male"), labels=c("HA", "RHA")) +
  guides(fill=guide_legend(title=NULL)) +
  theme_minimal() + ylab("SLA")

ggsave("/Users/emilydittmar/Google Drive/Active Projects/Transcriptomics_Exp/Analyses/Figures/SLA.png")



########## old code:

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
barplot(SLA_df$SLA_mean, 
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
####### MODELS ########
#########################

# ask whether any traits are significantly influenced by Nutrients, Salt and/or interaction