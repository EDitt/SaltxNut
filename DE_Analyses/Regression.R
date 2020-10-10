
#########################
######## SETUP ##########
#########################

library(dplyr)
library(ggplot2)
library(ggthemes)
source("Functions.R")


library(viridis)
library(RColorBrewer)

load("ResultsFiles/Overlap.RData")

lapply(SigOverlap, function(x) {length(x$Gene)})

#########################
###### DATAFRAME ########
#########################

# combine dataframes - log2FoldChange + Pvalues
my_columnsPL <- lapply(DEData, function(x) x[c(3,4,7)])
All_PLvalues <- do.call("cbind", my_columnsPL)
All_PLvalues <- cbind(DEData[[1]]["Gene"], All_PLvalues)

All_PLvalues[,c(1:10)]
length(All_PLvalues$Gene) #55,728

Sig_AllTrt <- All_PLvalues[which(All_PLvalues$condition_comboDE_results.padj < 0.05 |
                                   All_PLvalues$condition_nutDE_results.padj < 0.05 |
                                   All_PLvalues$condition_saltDE_results.padj < 0.05),]

length(Sig_AllTrt$Gene) #23,789

#########################
###### ALL SHARED #######
#########################

# Genes that are significant in all treatments and no pairwise differences among treatments

PairwiseDiffs <- union(my_dataSig$Nut_vs_Combo_results$Gene, 
                       union(my_dataSig$Nut_vs_Salt_results$Gene, my_dataSig$Salt_vs_Combo_results$Gene))

AllSame <- setdiff(SigOverlap$InCommonAll$Gene, PairwiseDiffs)
length(AllSame) #2627

# how many shared?
length(SigOverlap$InCommonAll$Gene) #5032

AllDiffs <- intersect(SigOverlap$InCommonAll$Gene, PairwiseDiffs)
length(AllDiffs) #2405

#########################
#### PLOT ALL SHARED ####
#########################


AllSameData <- Predictdf(All_PLvalues, AllSame, "condition_saltDE_results.log2FoldChange",
                              "condition_nutDE_results.log2FoldChange")

aggregate(AllSameData$Gene, by=list(AllSameData$Interval), length) #113 outside

AllDiffData <- Predictdf(All_PLvalues, AllDiffs, "condition_saltDE_results.log2FoldChange",
                         "condition_nutDE_results.log2FoldChange")

aggregate(AllDiffData$Gene, by=list(AllDiffData$Interval), length) #31 outside


pAll <- ggplot(data=AllSameData, aes(x=condition_nutDE_results.log2FoldChange, 
                                        y=condition_saltDE_results.log2FoldChange))
pAll + 
  geom_point(aes(fill=Interval), shape=21, stroke=0.0, size=2) +
  scale_fill_manual(values=c(alpha("blue", 0.5),
                             alpha("blue", 0.1))) +
  xlim(-5,5) + ylim(-5,5) 

pAll + geom_point(aes(shape=Interval, color=Interval, fill=Interval), size=2, stroke=0.5) +
  scale_shape_manual(values=c(21, 19)) +
  scale_color_manual(values=c(alpha("black", 0.7),
                              alpha("skyblue3", 0.1))) +
  scale_fill_manual(values=c(alpha("steelblue3", 0.5),
                             alpha("red", 0.1))) +
  xlim(-5,5) + ylim(-5,5) + 
  geom_smooth(method="lm", color="skyblue4", size=0.5) +
  #geom_smooth(method="lm", aes(y=lwr), linetype="dashed", color="skyblue3", size=0.5) + 
  #geom_smooth(method="lm", aes(y=upr), linetype="dashed", color="skyblue3", size=0.5) +
  geom_hline(yintercept=0, linetype="dotted", color="black", size=0.5) +
  geom_vline(xintercept=0, linetype="dotted", color="black", size=0.5) +
  geom_abline(slope=1, intercept=0, linetype="dotted", color="black", size=0.5) +
  
  geom_smooth(data=AllDiffData, aes(x=condition_nutDE_results.log2FoldChange,
                 y=condition_saltDE_results.log2FoldChange), method="lm")
  
pAllDiff <- ggplot(data=AllDiffData, aes(x=condition_nutDE_results.log2FoldChange,
                                         y=condition_saltDE_results.log2FoldChange)) 
pAllDiff+geom_point()+xlim(-10,10) + ylim(-10,10) 

### Plot both together?

pAll + geom_point(data=AllSameData[which(AllSameData$Interval=="Outside"),],
                  aes(x=condition_nutDE_results.log2FoldChange,
                      y=condition_saltDE_results.log2FoldChange), 
                  size=2, shape=21, stroke=0.5,
                  color=alpha("black", 0.5),
                  fill=alpha("skyblue3", 0.5)) +
  #xlim(-10,10) + ylim(-10,10) + 
  geom_smooth(method="lm", color="skyblue3", size=1) +
  geom_smooth(method="lm", aes(y=lwr), linetype="dashed", color="skyblue3", size=0.5) + 
  geom_smooth(method="lm", aes(y=upr), linetype="dashed", color="skyblue3", size=0.5) +
  geom_hline(yintercept=0, linetype="dotted", color="black", size=0.5) +
  geom_vline(xintercept=0, linetype="dotted", color="black", size=0.5) +
  #geom_abline(slope=1, intercept=0, linetype="dotted", color="black", size=0.5) +
  
  geom_smooth(data=AllDiffData, aes(x=condition_nutDE_results.log2FoldChange,
                                    y=condition_saltDE_results.log2FoldChange), 
              method="lm", size=1, color="goldenrod") +
  geom_point(data=AllDiffData[which(AllDiffData$Interval=="Outside"),], 
             aes(x=condition_nutDE_results.log2FoldChange,
             y=condition_saltDE_results.log2FoldChange),
             size=2, shape=21, stroke=0.5,
             color=alpha("black", 0.5),
             fill=alpha("goldenrod", 0.5)) +
  geom_smooth(data=AllDiffData, aes(x=condition_nutDE_results.log2FoldChange,
                                    y=lwr), 
              method="lm", linetype="dashed", color="goldenrod", size=0.5) +
  geom_smooth(data=AllDiffData, aes(x=condition_nutDE_results.log2FoldChange,
                                    y=upr), 
              method="lm", linetype="dashed", color="goldenrod", size=0.5) +
  theme_base(base_size = 12) +
  labs(x="Log2 fold-change Low Nutrient",
       y="Log2 fold-change High Salt", size=2)

#########################
####### ALL DIFFS #######
#########################

lapply(SigDiffOverlap, function(x) {length(x$Gene)})

SharedDiffs <- lapply(SigDiffOverlap, function(x) {intersect(x$Gene, AllDiffs)})
lapply(SharedDiffs, function(x) {length(x)})

# Total: 2405
# all pairs different, N=97 <- additive?
# Nut v. Combo, Nut v. Salt, N=1418 "salt dominant"
# Nut v. Salt, Salt v. Combo, N=73 "nut dominant"
# Nut v. Combo, Salt v. Combo, N=15 <- additive?
# Nut v. Combo, N=277
# Nut v. Salt, N=511 <- average?
# Salt v. Combo, N=14

# Salt > Nutrient
# *** how to illustrate these categories?


pAve <- ggplot(data=All_PLvalues[which(All_PLvalues$Gene %in% AllDiffs),], aes(x=NutSaltAve, 
                                    y=condition_comboDE_results.log2FoldChange))
pAve + geom_point() +
  #xlim(-10,10) + ylim(-10,10) +
  geom_abline(slope=1, intercept=0, linetype="dotted", color="black", size=0.5)

# diff between salt-nut
All_PLvalues$SaltminusNut <- All_PLvalues$condition_saltDE_results.log2FoldChange -
                            All_PLvalues$condition_nutDE_results.log2FoldChange

pDiff <- ggplot(data=All_PLvalues[which(All_PLvalues$Gene %in% AllDiffs),], aes(x=SaltminusNut, 
                                 y=condition_comboDE_results.log2FoldChange))
pDiff + geom_point() +
  xlim(-10,10) + ylim(-10,10)

pSN <- ggplot(data=All_PLvalues[which(All_PLvalues$Gene %in% AllDiffs),], 
              aes(x=condition_saltDE_results.log2FoldChange, 
                  y=condition_nutDE_results.log2FoldChange))
pSN + geom_point()

#########################
##### SALT DOMINANT #####
#########################

SaltDom <- Predictdf(All_PLvalues, SharedDiffs$`my_dataSig$Nut_vs_Combo_results[1]my_dataSig$Nut_vs_Salt_results[1]Only`, 
                     "condition_comboDE_results.log2FoldChange", "condition_saltDE_results.log2FoldChange")
aggregate(SaltDom$Gene, by=list(SaltDom$Interval), length) #28 outside
pSaltDom <- ggplot(data=SaltDom, aes(x=condition_saltDE_results.log2FoldChange, 
                                     y=condition_comboDE_results.log2FoldChange))
pSaltDom + geom_point() +
  #xlim(-10,10) + ylim(-10,10) +
  geom_abline(slope=1, intercept=0, linetype="dotted", color="black", size=0.5)

#########################
##### NUT DOMINANT ######
#########################

NutDom <- Predictdf(All_PLvalues, SharedDiffs$`my_dataSig$Nut_vs_Salt_results[1]my_dataSig$Salt_vs_Combo_results[1]Only`, 
                     "condition_comboDE_results.log2FoldChange", "condition_nutDE_results.log2FoldChange")
aggregate(NutDom$Gene, by=list(NutDom$Interval), length) #3 outside
pNutDom <- ggplot(data=NutDom, aes(x=condition_nutDE_results.log2FoldChange, 
                                     y=condition_comboDE_results.log2FoldChange))
pNutDom + geom_point() +
  #xlim(-10,10) + ylim(-10,10) +
  geom_abline(slope=1, intercept=0, linetype="dotted", color="black", size=0.5)

#########################
######## ADDITIVE #######
#########################

All_PLvalues$NutSaltAdd <- All_PLvalues$condition_nutDE_results.log2FoldChange +
  All_PLvalues$condition_saltDE_results.log2FoldChange

AddGenes <- union(SharedDiffs$InCommonAll, SharedDiffs$`my_dataSig$Nut_vs_Combo_results[1]my_dataSig$Salt_vs_Combo_results[1]Only`)

AllDiffAdd <- Predictdf(All_PLvalues, AddGenes, 
                        "condition_comboDE_results.log2FoldChange", "NutSaltAdd")
aggregate(AllDiffAdd$Gene, by=list(AllDiffAdd$Interval), length) # 5 outside

pAdd <- ggplot(data=AllDiffAdd, aes(x=NutSaltAdd, 
                                   y=condition_comboDE_results.log2FoldChange))
pAdd + geom_point() +
  xlim(-10,10) + ylim(-10,10) +
  geom_abline(slope=1, intercept=0, linetype="dotted", color="black", size=0.5)

#########################
## NUTRIENT CATEGORIES ##
#########################

# Genes that are significant in Nutrient + Combo only, no pairwise diffs
UnconditionNutGenes <- setdiff(SigOverlap$`my_dataSig$condition_comboDE_results[1]my_dataSig$condition_nutDE_results[1]Only`$Gene,
                               my_dataSig$Nut_vs_Combo_results$Gene)
length(UnconditionNutGenes) #2377

# DE in Nutrient + Combo only, pairwise differences OR DE in Nutrient only, pairwise differences
CondNutGenes <- union(SigOverlap$`my_dataSig$condition_nutDE_results[1]Only`$Gene,
                      intersect(my_dataSig$Nut_vs_Combo_results$Gene,
                      SigOverlap$`my_dataSig$condition_comboDE_results[1]my_dataSig$condition_nutDE_results[1]Only`$Gene))
length(CondNutGenes) #10256

length(SigOverlap$`my_dataSig$condition_nutDE_results[1]Only`$Gene) #8636
length(intersect(my_dataSig$Nut_vs_Combo_results$Gene,
      SigOverlap$`my_dataSig$condition_comboDE_results[1]my_dataSig$condition_nutDE_results[1]Only`$Gene)) #1620

#########################
### NUT UNCONDITIONAL ###
#########################

UncondNut <- Predictdf(All_PLvalues, UnconditionNutGenes, 
                  "condition_comboDE_results.log2FoldChange", "condition_nutDE_results.log2FoldChange")

pUNut <- ggplot(data=UncondNut, aes(x=condition_nutDE_results.log2FoldChange, 
                               y=condition_comboDE_results.log2FoldChange))

pUNut + geom_point(aes(shape=Interval, color=Interval, fill=Interval), size=2, stroke=0.5) +
  xlim(-10,10) + ylim(-10,10) +
  scale_shape_manual(values=c(21, 19)) +
  scale_color_manual(values=c(alpha("black", 0.7),
                              alpha("skyblue3", 0.1))) +
  scale_fill_manual(values=c(alpha("steelblue3", 0.5),
                             alpha("red", 0.1))) +
  geom_smooth(method="lm", size=0.5) +
  geom_smooth(aes(y=lwr), method="lm", linetype="dashed", color="skyblue3", size=0.5) + 
  geom_smooth(aes(y=upr), method="lm", linetype="dashed", color="skyblue3", size=0.5) +
  theme_base(base_size = 12) +
  labs(x="Log2 fold-change Low Nutrient",
       y="Log2 fold-change Low Nutrient + High Salt", size=2) +
  geom_hline(yintercept=0, linetype="dotted", color="black", size=0.5) +
  geom_vline(xintercept=0, linetype="dotted", color="black", size=0.5)
  

#########################
#### NUT CONDITIONAL ####
#########################

CondNut <- Predictdf(All_PLvalues, CondNutGenes, 
                       "condition_comboDE_results.log2FoldChange", "condition_nutDE_results.log2FoldChange")

pCNut <- ggplot(data=CondNut, aes(x=condition_nutDE_results.log2FoldChange, 
                                    y=condition_comboDE_results.log2FoldChange))

pCNut + geom_point(aes(shape=Interval, color=Interval, fill=Interval), size=2, stroke=0.5) +
  xlim(-30,20) + ylim(-30,20) +
  scale_shape_manual(values=c(21, 19)) +
  scale_color_manual(values=c(alpha("black", 0.7),
                              alpha("goldenrod", 0.1))) +
  scale_fill_manual(values=c(alpha("goldenrod", 0.5),
                             alpha("red", 0.1))) +
  geom_smooth(method="lm", color="goldenrod", size=0.5) +
  geom_smooth(aes(y=lwr), method="lm", linetype="dashed", color="goldenrod", size=0.5) + 
  geom_smooth(aes(y=upr), method="lm", linetype="dashed", color="goldenrod", size=0.5) +
  theme_base(base_size = 12) +
  labs(x="Log2 fold-change Low Nutrient",
       y="Log2 fold-change Low Nutrient + High Salt", size=2) +
  geom_hline(yintercept=0, linetype="dotted", color="black", size=0.5) +
  geom_vline(xintercept=0, linetype="dotted", color="black", size=0.5)


#########################
#### NUT UNCOND/COND ####
#########################


pUNut + geom_smooth(method="lm", size=0.5, color="skyblue3") +
  geom_smooth(aes(y=lwr), method="lm", linetype="dashed", color="skyblue3", size=0.5) + 
  geom_smooth(aes(y=upr), method="lm", linetype="dashed", color="skyblue3", size=0.5) +
  geom_point(data=UncondNut[which(UncondNut$Interval=="Outside"),],
             aes(x=condition_nutDE_results.log2FoldChange,
                 y=condition_comboDE_results.log2FoldChange),
             shape=21, color=alpha("black", 0.5), fill=alpha("steelblue3", 0.5)) +
  geom_point(data=CondNut[which(CondNut$Interval=="Outside"),],
             aes(x=condition_nutDE_results.log2FoldChange,
                 y=condition_comboDE_results.log2FoldChange),
             shape=21, color=alpha("black", 0.5), fill=alpha("goldenrod", 0.5)) +
  geom_smooth(data=CondNut, aes(x=condition_nutDE_results.log2FoldChange,
                 y=condition_comboDE_results.log2FoldChange),
              method="lm", size=0.5, color="goldenrod") +
  geom_smooth(data=CondNut, 
              aes(x=condition_nutDE_results.log2FoldChange,
                  y=lwr), method="lm", linetype="dashed", size=0.5, color="goldenrod") +
  geom_smooth(data=CondNut, 
              aes(x=condition_nutDE_results.log2FoldChange,
                  y=upr), method="lm", linetype="dashed", size=0.5, color="goldenrod") +
  xlim(-30,30) + ylim(-30,30) +
  theme_base(base_size = 12) +
  labs(x="Log2 fold-change Low Nutrient",
       y="Log2 fold-change Low Nutrient + High Salt", size=2) +
  geom_hline(yintercept=0, linetype="dotted", color="black", size=0.5) +
  geom_vline(xintercept=0, linetype="dotted", color="black", size=0.5) +
  geom_abline(slope=1, intercept=0, linetype="dotted", color="black", size=0.5)

#########################
###### NUT RESULTS ######
#########################

# Total: 12,633
# Unconditional: 2377 (18.8%)
# Conditional: 10256
  # Direction change: 1,276 (10.1%)
  # Magnitude higher in Combo: 52 (0.05%)
  # Magnitude lower in Combo: 8928 (87%)

NutDir <- DirectionDf(CondNut, "condition_comboDE_results.log2FoldChange",
                       "condition_nutDE_results.log2FoldChange")
lapply(NutDir, function(x) {length(x$Gene)})

# maybe I should exclude the ones that aren't significant in combo?
length(which(CondNut$condition_comboDE_results.padj >= 0.05)) #8487
length(which(CondNut$condition_comboDE_results.padj < 0.05)) #1620
length(which(is.na(CondNut$condition_comboDE_results.padj))) #149
length(CondNut$Gene) #10256

CondNutComboSig <- CondNut[which(CondNut$condition_comboDE_results.padj < 0.05),]

length(which(!CondNut$Gene %in% CondNutComboSig$Gene)) #8636

NutDir2 <- DirectionDf(CondNutComboSig, "condition_comboDE_results.log2FoldChange",
                      "condition_nutDE_results.log2FoldChange")
lapply(NutDir2, function(x) {length(x$Gene)})

# Total: 12,633
# Unconditional: 2377 (18.8%)
# Conditional: 10256
    # N.S. in Combo: 8,636 (84.2%)
    # Sig in Combo: 1,620 (15.8%)
      # Direction change: 38 (2.3%)
      # Magnitude higher in Combo: 6 (0.4%)
      # Magnitude lower in Combo: 1576 (97.3%)


#########################
#### SALT CATEGORIES ####
#########################

# Genes that are significant in Salt + Combo only, no pairwise diffs
UnconditionSaltGenes <- setdiff(SigOverlap$`my_dataSig$condition_comboDE_results[1]my_dataSig$condition_saltDE_results[1]Only`$Gene,
                               my_dataSig$Salt_vs_Combo_results$Gene)
length(UnconditionSaltGenes) #1142

# DE in Salt + Combo only, pairwise differences OR DE in Salt only
CondSaltGenes <- union(SigOverlap$`my_dataSig$condition_saltDE_results[1]Only`$Gene,
                      intersect(my_dataSig$Salt_vs_Combo_results$Gene,
                                SigOverlap$`my_dataSig$condition_comboDE_results[1]my_dataSig$condition_saltDE_results[1]Only`$Gene))
length(CondSaltGenes) #2088

# how many significant only in Salt treatment?
length(SigOverlap$`my_dataSig$condition_saltDE_results[1]Only`$Gene) #2018

# combo + salt only=1212
# Only 70 are significant in Combo+Salt overlap but different between them

CondSaltGenes2 <-intersect(my_dataSig$Salt_vs_Combo_results$Gene,
              SigOverlap$`my_dataSig$condition_comboDE_results[1]my_dataSig$condition_saltDE_results[1]Only`$Gene)

#########################
### SALT UNCONDITIONAL ###
#########################

UncondSalt <- Predictdf(All_PLvalues, UnconditionSaltGenes, 
                       "condition_comboDE_results.log2FoldChange", "condition_saltDE_results.log2FoldChange")

pUSalt <- ggplot(data=UncondSalt, aes(x=condition_saltDE_results.log2FoldChange, 
                                    y=condition_comboDE_results.log2FoldChange))

pUSalt + geom_point(aes(shape=Interval, color=Interval, fill=Interval), size=2, stroke=0.5) +
  xlim(-5,5) + ylim(-5,5) +
  scale_shape_manual(values=c(21, 19)) +
  scale_color_manual(values=c(alpha("black", 0.7),
                              alpha("skyblue3", 0.1))) +
  scale_fill_manual(values=c(alpha("steelblue3", 0.5),
                             alpha("red", 0.1))) +
  geom_smooth(method="lm", size=0.5) +
  geom_smooth(aes(y=lwr), method="lm", linetype="dashed", color="skyblue3", size=0.5) + 
  geom_smooth(aes(y=upr), method="lm", linetype="dashed", color="skyblue3", size=0.5) +
  theme_base(base_size = 12) +
  labs(x="Log2 fold-change High Salt",
       y="Log2 fold-change Low Nutrient + High Salt", size=2) +
  geom_hline(yintercept=0, linetype="dotted", color="black", size=0.5) +
  geom_vline(xintercept=0, linetype="dotted", color="black", size=0.5)


#########################
#### SALT CONDITIONAL ####
#########################

CondSalt <- Predictdf(All_PLvalues, CondSaltGenes, 
                     "condition_comboDE_results.log2FoldChange", "condition_saltDE_results.log2FoldChange")

pCSalt <- ggplot(data=CondSalt, aes(x=condition_saltDE_results.log2FoldChange, 
                                  y=condition_comboDE_results.log2FoldChange))

pCSalt2 + geom_point(aes(shape=Interval, color=Interval, fill=Interval), size=2, stroke=0.5) +
  #xlim(-30,20) + ylim(-30,20) +
  scale_shape_manual(values=c(21, 19)) +
  scale_color_manual(values=c(alpha("black", 0.7),
                              alpha("goldenrod", 0.1))) +
  scale_fill_manual(values=c(alpha("goldenrod", 0.5),
                             alpha("red", 0.1))) +
  geom_smooth(method="lm", color="goldenrod", size=0.5) +
  geom_smooth(aes(y=lwr), method="lm", linetype="dashed", color="goldenrod", size=0.5) + 
  geom_smooth(aes(y=upr), method="lm", linetype="dashed", color="goldenrod", size=0.5) +
  theme_base(base_size = 12) +
  labs(x="Log2 fold-change High Salt",
       y="Log2 fold-change Low Nutrient + High Salt", size=2) +
  geom_hline(yintercept=0, linetype="dotted", color="black", size=0.5) +
  geom_vline(xintercept=0, linetype="dotted", color="black", size=0.5)

# what about 70 that are different but significant in combo?

CondSalt2 <- Predictdf(All_PLvalues, CondSaltGenes2, 
                      "condition_comboDE_results.log2FoldChange", "condition_saltDE_results.log2FoldChange")

pCSalt2 <- ggplot(data=CondSalt2, aes(x=condition_saltDE_results.log2FoldChange, 
                                    y=condition_comboDE_results.log2FoldChange))

# 4 change direction - upregulated in salt + nutrient (though down in salt only)

#########################
###### SALT RESULTS ######
#########################

SaltDir <- DirectionDf(CondSalt2, "condition_comboDE_results.log2FoldChange",
                       "condition_saltDE_results.log2FoldChange")
lapply(SaltDir, function(x) {length(x$Gene)})

# Total: 3,230
  # Unconditional: 1,142 (35.3%)
  # Conditional: 2,088 (64.6%)
    # N.S. in Combo: 2,018 (96.6%%)
    # Sig in Combo: 70 (3.35%)
      # Direction change: 7 (10%) *all up in combo
      # Magnitude higher in Combo: 1 (1.4%)
      # Magnitude lower in Combo: 62 (88.6%)

#########################
##### ANTAGONISTIC ######
#########################

# significant in nutrient + salt & opposite directions

my_dataSigUp <- lapply(my_dataSig, MoreCritNum, column=3, critNum=0)
my_dataSigDown <- lapply(my_dataSig, LessCritNum, column=3, critNum=0)

UpSaltDownNut <- intersect(my_dataSigUp$condition_saltDE_results$Gene,
                           my_dataSigDown$condition_nutDE_results$Gene)
length(UpSaltDownNut) #366

DownSaltUpNut <- intersect(my_dataSigDown$condition_saltDE_results$Gene,
                           my_dataSigUp$condition_nutDE_results$Gene)
length(DownSaltUpNut) #389

AntagonisticGenes <- union(UpSaltDownNut, DownSaltUpNut)

# how many are also significant in combo?
AntagonisticGenes_ComboSig <-intersect(AntagonisticGenes, my_dataSig$condition_comboDE_results$Gene)
length(AntagonisticGenes_ComboSig) #149

#SharedDiffs <- 

lapply(SigDiffOverlap, function(x) {length(intersect(x$Gene, AntagonisticGenes_ComboSig))})
# all pairwise differences: 24 -> "Additive?
# Nut vs Combo, Nut vs. Salt: 121 -> "Salt dominant"
# Nut vs. Salt, Salt vs. Combo: 3 -> "Nutrient dominant"
# Nut vs. Combo, Salt vs. Combo: 0
# Nut vs. Combo only: 0
# Nut vs. Salt only: 1
# Salt vs. Combo only: 0


AntData <- All_PLvalues[which(All_PLvalues$Gene %in% AntagonisticGenes),]
AntDataComSig <- All_PLvalues[which(All_PLvalues$Gene %in% AntagonisticGenes_ComboSig),]
pAnt <- ggplot(data=AntData, aes(x=NutSaltAdd, 
                                  y=condition_comboDE_results.log2FoldChange))
pAnt + geom_point()

pAnt2 <- ggplot(data=AntDataComSig, aes(x=NutSaltAdd, 
                                 y=condition_comboDE_results.log2FoldChange))
pAnt2 + geom_point()

pAnt3 <- ggplot(data=AntData, aes(x=NutSaltAve, 
                                 y=condition_comboDE_results.log2FoldChange))
pAnt3 + geom_point()

##################################
# are any genes in combo a diff direction than both nut + salt


length(intersect(my_dataSigUp$condition_comboDE_results$Gene,
                 intersect(my_dataSigDown$condition_nutDE_results$Gene,
                           my_dataSigDown$condition_saltDE_results$Gene))) #0
ComboDown_otherUp <- (intersect(my_dataSigDown$condition_comboDE_results$Gene,
                 intersect(my_dataSigUp$condition_nutDE_results$Gene,
                           my_dataSigUp$condition_saltDE_results$Gene))) #1



