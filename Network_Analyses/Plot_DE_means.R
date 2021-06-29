### graphs for interesting/important modules


#########################
######## SETUP ##########
#########################

library(dplyr)
source("Functions.R")
library(ggplot2)
library(tidyr)
#library(cowplot)
library(gridExtra)

#########################
###### READ DATA ########
#########################

AllData <- read.csv("ResultsFiles/Coexpression/Module_EGs.csv", header = T)
AllData$Treatment <- factor(AllData$Treatment, levels=c("Control", "LowNut", "HighSalt", "Combo"))
levels(AllData$Treatment)
AllData$Accession <- as.factor(make.names(AllData$Accession))
levels(Design$Accession)

# load ANOVA results
LM_Results <- read.csv("ResultsFiles/Coexpression/Module_Anova.csv", header=T) # 88 modules

# load gene lists
load("ResultsFiles/Coexpression/LR_SigModOverlap.RData")
# objects: SigModOverlap, DE_Overlaps, SigDiffOverlap, LR_list, DE_Combo_sig, DE_Salt_sig, DE_Nut_sig

###########################
## PLOT AVE TREATMENT DE ##
###########################

# DE (average pairwise difference between Control + the other 3 treatments)
# requires ANOVA results df (LM_Results)

Plot_ModDE <- function(ModuleName, dataframe, color, title){
  means <- as.data.frame(t(dataframe[which(dataframe$Module==ModuleName),c(1,18,21,24)]))
  se <- as.data.frame(t(dataframe[which(dataframe$Module==ModuleName),c(1,19,22,25)]))
  pvals <- as.data.frame(t(dataframe[which(dataframe$Module==ModuleName),c(1,20,23,26)]))
  Mod_df <- data.frame(Treatment=c("Nutrient", "Salt", "Combo"), 
                       DE=as.numeric(means[2:4,1]),
                       SE=as.numeric(se[2:4,1]),
                       pval=as.numeric(pvals[2:4,1]))
  Mod_df$Treatment <- factor(Mod_df$Treatment, levels=c("Nutrient", "Salt", "Combo"))
  Mod_df$plabs <- ifelse(Mod_df$pval < 0.0001, "***",
                        ifelse(Mod_df$pval < 0.001, "**",
                               ifelse(Mod_df$pval < 0.01, "*",
                                      paste0("p=",round(Mod_df$pval, digits=2)))))
  scale <- abs(max(Mod_df$DE) - min(Mod_df$DE))
  Mod_df$plabPositions <- ifelse(Mod_df$pval < 0.01,
                                 Mod_df$DE + sign(Mod_df$DE) * Mod_df$SE +
                                   sign(Mod_df$DE) * (scale/8),              
                                 Mod_df$DE + sign(Mod_df$DE) * Mod_df$SE +
                                   sign(Mod_df$DE) * (scale/5))
  ymin <- min(Mod_df$DE) - (2*max(Mod_df$SE) + (scale/5))
  #ymin <- min(Mod_df$DE) - scale/5
  ymin_plot <- ifelse(ymin > 0, 0, ymin)
  ymax <- max(Mod_df$DE) + (2*max(Mod_df$SE) + (scale/5))
  #ymax <- max(Mod_df$DE) + scale/5
  ymax_plot <- ifelse(ymax < 0, 0, ymax)
  p <- ggplot(data = Mod_df, aes(x=Treatment, y=DE)) +
    geom_bar(stat="identity", fill=color, alpha=0.7, color="black") +
    geom_errorbar(aes(ymin=DE-SE, ymax=DE+SE), width=0.05) +
    #geom_linerange(aes(ymin=DE-SE, ymax=DE+SE)) +
    theme_minimal() +
    ggtitle(title) +
    geom_text(aes(label=plabs, y=plabPositions)) +
              #vjust="outward") +
    ylim(ymin_plot, ymax_plot) 
  return(p)
}

#Mod_df$plabPositions <- ifelse(Mod_df$pval < 0.01,
#                               Mod_df$DE + sign(Mod_df$DE) * Mod_df$SE +
#                                 sign(Mod_df$DE) * (Mod_df$SE),                
#                               Mod_df$DE + sign(Mod_df$DE) * Mod_df$SE +
#                                sign(Mod_df$DE) * (0.75*Mod_df$SE))

#########################
# DE SALT + DE NUT ONLY #
#########################

lapply(SigModOverlap, function(x) {intersect(x, 
                                       DE_Overlaps$DE_Salt_sigDE_Nut_sigOnly)})

# all 3 overlap with Nut + Salt sig only (which has 7 in that category)

SaltNut <- lapply(DE_Overlaps$DE_Salt_sigDE_Nut_sigOnly, function(x) {
  Plot_ModDE(ModuleName=x,
             dataframe = LM_Results,
             color=gsub("ME", "", x), 
             title=gsub("ME", "", x))
}
)

grid.arrange(grobs=SaltNut, ncol=2)
ggsave(file = "ResultsFiles/Coexpression/barplots/SaltNut.png", 
       arrangeGrob(grobs = SaltNut, ncol = 2))

plot_grid(plotlist = SaltNut, ncol=2) # cowplot solution

#########################
# DE NUT + DE COMBO ONLY #
#########################

length(DE_Overlaps$DE_Combo_sigDE_Nut_sigOnly) #15
intersect(DE_Overlaps$DE_Combo_sigDE_Nut_sigOnly, lm_list$Combo_Nut) #14 (other than MEyellow4)
intersect(DE_Overlaps$DE_Combo_sigDE_Nut_sigOnly, LR_list$NutxSalt) #13
setdiff(DE_Overlaps$DE_Combo_sigDE_Nut_sigOnly, LR_list$NutxSalt) #plum2, yellow4
intersect(DE_Overlaps$DE_Combo_sigDE_Nut_sigOnly, LR_list$NutMain) #all
intersect(DE_Overlaps$DE_Combo_sigDE_Nut_sigOnly, LR_list$SaltMain) #6
# coral, indianred4, palevioletred, purple, skyblue, steelblue

NutCombo <- lapply(DE_Overlaps$DE_Combo_sigDE_Nut_sigOnly, function(x) {
  Plot_ModDE(ModuleName=x,
             dataframe = LM_Results,
             color=gsub("ME", "", x), 
             title=gsub("ME", "", x))
}
)
grid.arrange(grobs=NutCombo)
ggsave(file = "ResultsFiles/Coexpression/barplots/NutCombo.png", 
       arrangeGrob(grobs = NutCombo))

#########################
# DE SALT + DE COMBO ONLY #
#########################

length(DE_Overlaps$DE_Combo_sigDE_Salt_sigOnly) #1
intersect(DE_Overlaps$DE_Combo_sigDE_Salt_sigOnly, lm_list$Combo_Salt) #0
intersect(DE_Overlaps$DE_Combo_sigDE_Salt_sigOnly, LR_list$NutxSalt) #0
intersect(DE_Overlaps$DE_Combo_sigDE_Salt_sigOnly, LR_list$NutMain) #0

SaltCombo <- Plot_ModDE(DE_Overlaps$DE_Combo_sigDE_Salt_sigOnly,
                   LM_Results,
                   color=gsub("ME", "", DE_Overlaps$DE_Combo_sigDE_Salt_sigOnly), 
                   title=gsub("ME", "", DE_Overlaps$DE_Combo_sigDE_Salt_sigOnly))

ggsave(file = "ResultsFiles/Coexpression/barplots/SaltCombo.png")

#########################
##### DE SALT ONLY ######
#########################

length(DE_Overlaps$DE_Salt_sigOnly) #2

intersect(DE_Overlaps$DE_Salt_sigOnly, lm_list$Combo_Salt) #1 (greenyellow)
intersect(DE_Overlaps$DE_Salt_sigOnly, LR_list$NutxSalt) #1 (greenyellow)
intersect(DE_Overlaps$DE_Salt_sigOnly, LR_list$NutMain) #0
intersect(DE_Overlaps$DE_Salt_sigOnly, LR_list$SaltMain) #1 (black)

Salt <- lapply(DE_Overlaps$DE_Salt_sigOnly, function(x) {
  Plot_ModDE(ModuleName=x,
             dataframe = LM_Results,
             color=gsub("ME", "", x), 
             title=gsub("ME", "", x))
}
)

grid.arrange(grobs=Salt)

ggsave(file = "ResultsFiles/Coexpression/barplots/SaltOnly.png", 
       arrangeGrob(grobs = Salt))

#########################
###### DE NUT ONLY ######
#########################
length(DE_Overlaps$DE_Nut_sigOnly) #10
intersect(DE_Overlaps$DE_Nut_sigOnly, lm_list$Combo_Nut) #5
intersect(DE_Overlaps$DE_Nut_sigOnly, LR_list$NutxSalt) #0
intersect(DE_Overlaps$DE_Nut_sigOnly, LR_list$SaltMain)
intersect(DE_Overlaps$DE_Nut_sigOnly, LR_list$Accession) #7

Nut <- lapply(DE_Overlaps$DE_Nut_sigOnly, function(x) {
  Plot_ModDE(ModuleName=x,
             dataframe = LM_Results,
             color=gsub("ME", "", x), 
             title=gsub("ME", "", x))
}
)

grid.arrange(grobs=Nut)

#########################
## PLOT MODEL EFFECTS? ##
#########################

SigModOverlap
# 6 in common all, 7 nut + salt only, 0 salt + nutxsalt, 19 nut + nutxsalt,
# 8 nut only, 2 salt only, 1 nutxsalt only

lapply(SigModOverlap, function(x) {intersect(x, DE_Overlaps$InCommonAll)})
# all 12 overlap with nut+nutxsalt

lapply(DE_Overlaps, function(x) {intersect(x, SigModOverlap$LR_Nut_SigLR_NutxSalt_SigOnly)})
# remaining 7 are in the combo+nut only

#########################
##### LR ALL FACTORS ####
#########################

lapply(DE_Overlaps, function(x) {intersect(x, SigModOverlap$InCommonAll)})
# all 6 fall into "DE in combo and nutrient only"

lapply(SigDiffOverlap, function(x) {intersect(x, SigModOverlap$InCommonAll)})
# all 6 are different between all pairs

# sort by Salt p-values (ascending order)
LM_Results_SaltpSort <- LM_Results[order(LM_Results$Difference_p.DE_Salt,
                                           decreasing = FALSE),]

AllModEffects <- lapply(LM_Results_SaltpSort[which(LM_Results_SaltpSort$Module %in%
                                                     SigModOverlap$InCommonAll),"Module"], 
                        function(x) {
            Plot_ModDE(ModuleName=x,
             dataframe = LM_Results_SaltpSort,
             color=gsub("ME", "", x), 
             title=gsub("ME", "", x))
}
)

grid.arrange(grobs=AllModEffects, top = "Nutrient, Salt, and Nutrient x Salt Significant")

ggsave(file = "ResultsFiles/Coexpression/barplots/AllModEffectsSig.png", 
       arrangeGrob(grobs = AllModEffects,
                   top = "Nutrient, Salt, and Nutrient x Salt Significant"))

# all pairwise diffs are significant for this category
# only significantly different from control for nutrient & combo

##########################
# IN COMMON NUT+NUTxSALT #
##########################

SigModOverlap$LR_Nut_SigLR_NutxSalt_SigOnly # N=19
lapply(DE_Overlaps, function(x) {intersect(x, SigModOverlap$LR_Nut_SigLR_NutxSalt_SigOnly)})
# 12 are significant in all 3; 7 are significant only in combo + nutrient

diffs_overlap <- lapply(SigDiffOverlap, function(x) {intersect(x, 
                                              SigModOverlap$LR_Nut_SigLR_NutxSalt_SigOnly)})
# 5 are different among all pairs
# 10 are different between combo-nut & salt-nut (not diff between combo-salt)

no_diffs <- setdiff(SigModOverlap$LR_Nut_SigLR_NutxSalt_SigOnly, union(Combo_Nut_sig, 
                                                 union(Combo_Salt_sig, Salt_Nut_sig)))
# no pairwise differences: blue2, darkorange2, mediumpurple2, salmon
# these ^ 4 are all DE

# which of the ones that are n.s. in salt intersect w/ pairwise diffs?
lapply(SigDiffOverlap, function(x) {intersect(x, 
                                              intersect(SigModOverlap$LR_Nut_SigLR_NutxSalt_SigOnly,
                                                        DE_Overlaps$DE_Combo_sigDE_Nut_sigOnly))})
# 1 (skyblue3) is diff between all,
# 6 are diff between nut-salt only

# sort by:
#  use salt p-value sorted
# add sort (decreasing=TRUE) for pairwise difference p-values:
LM_Results_SaltpDiffpSort <- LM_Results_SaltpSort[order(LM_Results_SaltpSort$Difference_p.Nut.Combo,
                                                        LM_Results_SaltpSort$Difference_p.Nut.Salt,
                                                        LM_Results_SaltpSort$Difference_p.Salt.Combo,
                                         decreasing = TRUE),]


LRNutInter <- lapply(LM_Results_SaltpDiffpSort[which(LM_Results_SaltpDiffpSort$Module %in%
                                                   SigModOverlap$LR_Nut_SigLR_NutxSalt_SigOnly),"Module"], 
                      function(x) {
            Plot_ModDE(ModuleName=x,
             dataframe = LM_Results,
             color=gsub("ME", "", x), 
             title=gsub("ME", "", x))
}
)

# or sort by pairwise diffs?
ModOrder <- union(no_diffs, union(diffs_overlap$Combo_Nut_sigSalt_Nut_sigOnly,
                                  diffs_overlap$InCommonAll))

LRNutInter <- lapply(ModOrder, 
                     function(x) {
                       Plot_ModDE(ModuleName=x,
                                  dataframe = LM_Results,
                                  color=gsub("ME", "", x), 
                                  title=gsub("ME", "", x))
                     }
)


grid.arrange(grobs=LRNutInter)
grid.arrange(grobs=LRNutInter[1:12], 
             top = "Nutrient + Nutrient x Salt significant; DE in all 3")

grid.arrange(grobs=LRNutInter[13:19], 
             top = "Nutrient + Nutrient x Salt significant; DE in Nut+Combo")

#ggsave(file = "ResultsFiles/Coexpression/barplots/NutInterEffectsSig.png", 
#       arrangeGrob(grobs = LRNutInter))

# DE in all 3:
which(ModOrder %in% DE_Overlaps$InCommonAll)
grid.arrange(grobs=LRNutInter[which(ModOrder %in% DE_Overlaps$InCommonAll)])
ggsave(file = "ResultsFiles/Coexpression/barplots/NutInterEffectsSig_AllDE.png", 
       arrangeGrob(grobs = LRNutInter[which(ModOrder %in% DE_Overlaps$InCommonAll)],
                   top = "Nutrient + NutxSalt Interaction significant"))

ggsave(file = "ResultsFiles/Coexpression/barplots/NutInterEffectsSig_ComboNutDE.png", 
       arrangeGrob(grobs = LRNutInter[which(!ModOrder %in% DE_Overlaps$InCommonAll)],
                   top = "Nutrient + NutxSalt Interaction significant"))


###################################
# IN COMMON NUT+SALT MAIN EFFECTS #
###################################

SigModOverlap$LR_Nut_SigLR_Salt_SigOnly # N=7
lapply(DE_Overlaps, function(x) {intersect(x, SigModOverlap$LR_Nut_SigLR_Salt_SigOnly)})
# Salt + Nut only: 3, Nut only: 4

lapply(SigDiffOverlap, function(x) {intersect(x, 
                                              SigModOverlap$LR_Nut_SigLR_Salt_SigOnly)})
# 4 in common all, 3 diff b/w combo/nut + salt/nut (darkgrey, lightsteelblue1, plum1)

# sort by Nutrient + Salt p-values
LM_Results_SaltNutComboPsort <- LM_Results[order(LM_Results$Difference_p.DE_Salt,
                                             LM_Results$Difference_p.DE_Nut,
                                             LM_Results$Difference_p.DE_Combo,
                                         decreasing = FALSE),]

LRNutSaltMain <- lapply(LM_Results_SaltNutComboPsort[which(LM_Results_SaltNutComboPsort$Module %in%
                                                             SigModOverlap$LR_Nut_SigLR_Salt_SigOnly),"Module"], 
                     function(x) {
                       Plot_ModDE(ModuleName=x,
                                  dataframe = LM_Results,
                                  color=gsub("ME", "", x), 
                                  title=gsub("ME", "", x))
                     }
)

grid.arrange(grobs=LRNutSaltMain, 
             top = "Nutrient + Salt Main effects significant")

ggsave(file = "ResultsFiles/Coexpression/barplots/NutSaltEffectsSig.png", 
       arrangeGrob(grobs = LRNutSaltMain,
                   top = "Nutrient + Salt Main effects significant"))

###################################
###### NUT MAIN EFFECTS ONLY ######
###################################

SigModOverlap$LR_Nut_SigOnly # N=8
lapply(DE_Overlaps, function(x) 
  {intersect(x, SigModOverlap$LR_Nut_SigOnly)})
# Combo + nut only: 2 (plum2 & yellow4)
# Nut only: 6
lapply(SigDiffOverlap, function(x) {intersect(x, 
                                              SigModOverlap$LR_Nut_SigOnly)})
# Combo-Salt & Salt-Nut:4
# Combo-Nut & Salt-Nut: 2
# Salt-Nut: 1

# how many are different between nut + combo?
intersect(SigModOverlap$LR_Nut_SigOnly,
          lm_list$Combo_Nut) # lightpink4 and plum2

# sort by combo & nutrient p-values
LM_Results_ComboNutPsort <- LM_Results[order(LM_Results$Difference_p.DE_Combo,
                                                 LM_Results$Difference_p.DE_Nut,
                                                 decreasing = FALSE),]

LRNutMain <- lapply(LM_Results_ComboNutPsort[which(LM_Results_ComboNutPsort$Module %in%
                                                     SigModOverlap$LR_Nut_SigOnly),"Module"], 
                        function(x) {
                          Plot_ModDE(ModuleName=x,
                                     dataframe = LM_Results,
                                     color=gsub("ME", "", x), 
                                     title=gsub("ME", "", x))
                        }
)

grid.arrange(grobs=LRNutMain, 
             top = "Nutrient Main effect significant")

ggsave(file = "ResultsFiles/Coexpression/barplots/NutEffectSig.png", 
       arrangeGrob(grobs = LRNutMain,
                   top = "Nutrient Main effect significant"))

###################################
###### SALT MAIN EFFECTS ONLY #####
###################################

SigModOverlap$LR_Salt_SigOnly # N=2
lapply(DE_Overlaps, function(x) 
{intersect(x, SigModOverlap$LR_Salt_SigOnly)})
# Combo + salt only: 1; salt only: 1

# how many are different between salt + combo?
intersect(SigModOverlap$LR_Salt_SigOnly,
          lm_list$Combo_Salt) #0

# Used sorted Salt + combo p-values

LRSaltMain <- lapply(LM_Results_SaltNutComboPsort[which(LM_Results_SaltNutComboPsort$Module %in%
                                                     SigModOverlap$LR_Salt_SigOnly),"Module"], 
                    function(x) {
                      Plot_ModDE(ModuleName=x,
                                 dataframe = LM_Results,
                                 color=gsub("ME", "", x), 
                                 title=gsub("ME", "", x))
                    }
)

grid.arrange(grobs=LRSaltMain, 
             top = "Salt Main effect significant")
ggsave(file = "ResultsFiles/Coexpression/barplots/SaltEffectSig.png", 
       arrangeGrob(grobs = LRSaltMain,
                   top = "Salt Main effect significant"))

###################################
######## INTERACTION ONLY #########
###################################

SigModOverlap$LR_NutxSalt_SigOnly # N=1
lapply(DE_Overlaps, function(x) 
{intersect(x, SigModOverlap$LR_NutxSalt_SigOnly)}) # sig in Salt only


LR_interaction <- Plot_ModDE(SigModOverlap$LR_NutxSalt_SigOnly,
                        LM_Results,
                        color=gsub("ME", "", SigModOverlap$LR_NutxSalt_SigOnly), 
                        title=gsub("ME", "", SigModOverlap$LR_NutxSalt_SigOnly))
LR_interaction + ggtitle("Interaction only Significant")

ggsave(file = "ResultsFiles/Coexpression/barplots/NutxSaltSig.png")

###################################
###### DE BUT NOT LR MODEL? #######
###################################

# are there any modules that are significantly DE but not represented as significant?

setdiff()
AllSigTreatMods



##########
barplot(AllData_ordered$MEyellow4, names.arg = AllData_ordered$Treatment)


##########
########## SCRATCH
p <- ggplot(data = AllData, aes(x=Treatment, y=MEplum1, group=Accession))
p + geom_bar(stat = "identity", position = position_dodge(), fill="grey", color="black")
#aes(fill=Accession))

plum1_means <- as.data.frame(t(LM_Results[which(LM_Results$Module=="MEplum1"),c(1,18,21,24)]))
plum1_se <- as.data.frame(t(LM_Results[which(LM_Results$Module=="MEplum1"),c(1,19,22,25)]))
plum1_pvals <- as.data.frame(t(LM_Results[which(LM_Results$Module=="MEplum1"),c(1,20,23,26)]))

plum1 <- data.frame(Treatment=c("Nutrient", "Salt", "Combo"), 
                    DE=as.numeric(plum1_means[2:4,1]),
                    SE=as.numeric(plum1_se[2:4,1]),
                    pval=as.numeric(plum1_pvals[2:4,1]))
plum1$Treatment <- factor(plum1$Treatment, levels=c("Nutrient", "Salt", "Combo"))
plum1$plabs <- ifelse(plum1$pval < 0.0001, "***",
                      ifelse(plum1$pval < 0.001, "**",
                             ifelse(plum1$pval < 0.01, "*",
                                    paste0("p=",round(plum1$pval, digits=2)))))
plum1$plabPositions <- plum1$DE + sign(plum1$DE) * plum1$SE +
                                      sign(plum1$DE) * (0.25*plum1$SE) # add extra cushion

ymin <- min(plum1$DE) - 2*max(plum1$SE)
ymin_plot <- ifelse(ymin > 0, 0, ymin)
ymax <- max(plum1$DE) + 2*max(plum1$SE)
ymax_plot <- ifelse(ymax < 0, 0, ymax)


scale <- abs(max(plum1$DE) - min(plum1$DE))
0.2/0.02 #10
0.2/0.04 #5
0.2/10 # 0.02

p2 <- ggplot(data = plum1, aes(x=Treatment, y=DE))
p2 + geom_bar(stat="identity", fill="skyblue", alpha=0.7, color="black") +
  geom_errorbar(aes(ymin=DE-SE, ymax=DE+SE), width=0.05) +
  #geom_linerange(aes(ymin=DE-SE, ymax=DE+SE)) +
  theme_minimal() +
  ggtitle("Plum1") +
  geom_text(aes(label=plabs, y=plabPositions),
            vjust="outward") +
  ylim(ymin, ymax) 


labs(title="plum1_module", vjust=0.5)

# put p-values over bars

# default scale: -0.1 - 0.2

testplum <- Plot_ModDE("MEplum1", LM_Results, "skyblue", "Plum1")
