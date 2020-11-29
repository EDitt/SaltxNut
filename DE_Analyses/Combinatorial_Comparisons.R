#########################
######## SETUP ##########
#########################

library(dplyr)
library(ggplot2)
source("Functions.R")
library(UpSetR)


library(viridis)
library(RColorBrewer)

load("ResultsFiles/Overlap.RData")
# DEData: All data for all pairwise contrasts
# SigOverlap: Overlap between DE genes for 3 treatments
# SigDiffOverlap: Overlap between sig. pairwise differences: Combo-Nut, Combo-Salt, Nut-Salt

lapply(SigOverlap, function(x) {length(x$Gene)})
lapply(SigDiffOverlap, function(x) {length(x$Gene)})

my_dataSig <- lapply(DEData, SigDEdf, PvaluesCol=7, CritP=0.05) #list that contains only the genes that are significant

# combine dataframes - log2FoldChange + Pvalues
my_columnsPL <- lapply(DEData, function(x) x[c(3,4,7)])
All_PLvalues <- do.call("cbind", my_columnsPL)
All_PLvalues <- cbind(DEData[[1]]["Gene"], All_PLvalues)


#########################
##### NUTRIENT ONLY #####
#########################

# (Not significant in Combo or Salt)

# N=8636

# up vs. down
Nut_only <- my_dataSig$condition_nutDE_results[which(my_dataSig$condition_nutDE_results$Gene %in%
                                                       SigOverlap$`my_dataSig$condition_nutDE_results[1]Only`$Gene),]
length(which(Nut_only$log2FoldChange > 0)) #4628
length(which(Nut_only$log2FoldChange < 0)) #4008


#########################
##### SALT ONLY #####
#########################

# (Not significant in Combo or Nutrient)

# N=2018

# up vs. down
Salt_only <- my_dataSig$condition_saltDE_results[which(my_dataSig$condition_saltDE_results$Gene %in%
                                                       SigOverlap$`my_dataSig$condition_saltDE_results[1]Only`$Gene),]
length(which(Salt_only$log2FoldChange > 0)) #1110
length(which(Salt_only$log2FoldChange < 0)) #908

#########################
# NUTRIENT & COMBO ONLY #
#########################

Nut_Combo <- Compare_shared(SigOverlap$`my_dataSig$condition_comboDE_results[1]my_dataSig$condition_nutDE_results[1]Only`$Gene,
                            my_dataSig$Nut_vs_Combo_results$Gene)

lapply(Nut_Combo, function(x) {length(x)})
# 2377 shared, same; 1620 shared different


##### Of the ones that are different, which are: higher in combo, lower in combo, opposite in combo?

NutCombo_df <- merge(my_dataSig$condition_nutDE_results, my_dataSig$condition_comboDE_results, by="Gene")
length(NutCombo_df$Gene) #9029

NutComboDiff <- NutCombo_df[which(NutCombo_df$Gene %in% Nut_Combo$Shared_different),]
length(NutComboDiff$Gene) #to check


NutComboDiff_compare <- Treat_compare(NutComboDiff, "log2FoldChange.x", "log2FoldChange.y", "Nutrient", "Combo")
lapply(NutComboDiff_compare, function(x) {length(x$Gene)})
# up in Nutrient, reduced in combo: 738 (increased mag. in combo: 5)
# down in Nutrient, reduced in combo: 838 (increased mag. in combo: 1)
# Opposite directions- 14 (Up in nutrient) + 24 (Down in Nutrient)

#########################
### SALT & COMBO ONLY ###
#########################

Salt_Combo <- Compare_shared(SigOverlap$`my_dataSig$condition_comboDE_results[1]my_dataSig$condition_saltDE_results[1]Only`$Gene,
                             my_dataSig$Salt_vs_Combo_results$Gene)
lapply(Salt_Combo, function(x) {length(x)})
# 1142 shared, same; 70 shared different

##### Of the ones that are different, which are: higher in combo, lower in combo, opposite in combo?
SaltCombo_df <- merge(my_dataSig$condition_saltDE_results, my_dataSig$condition_comboDE_results, by="Gene")
length(SaltCombo_df$Gene) #6244

SaltComboDiff <- SaltCombo_df[which(SaltCombo_df$Gene %in% Salt_Combo$Shared_different),]
length(SaltComboDiff$Gene) #70

SaltComboDiff_compare <- Treat_compare(SaltComboDiff, "log2FoldChange.x", "log2FoldChange.y", "Salt", "Combo")
lapply(SaltComboDiff_compare, function(x) {length(x$Gene)})
# up in Salt, reduced in combo: 24 (increased mag. in combo: 1)
# down in Salt, reduced in combo: 38 (incrased mag. in combo: 0)
# Opposite directions: 0 (Up in salt), 7 (down in salt)


#########################
##### ALL 3 SHARED ######
#########################

lapply(SigOverlap, function(x) {length(x$Gene)})
# 5032 are in common among all 3
# 1855 are in nutrient + salt, but not combo

# how many are shared and neither nut or salt are different from combo?
Shared_all <- Compare_shared(SigOverlap$InCommonAll$Gene,
                             union(my_dataSig$Nut_vs_Combo_results$Gene,
                                   my_dataSig$Salt_vs_Combo_results$Gene))
lapply(Shared_all, function(x) {length(x)})
# 3138 shared same; 1894 shared, different


# how many in common among all 3 have different combinations of pairwise differences?
InCommonDiffs <- lapply(SigDiffOverlap, function(x) {intersect(x$Gene, SigOverlap$InCommonAll$Gene)})
lapply(InCommonDiffs, function(x) {length(x)})
# 97 - all 3 pairwise differences
# 1418 - nutrient different from other 2
# 73 - salt different from other 2
# 15 - combo different from other 2
# 277 - nut v. combo only
# 511 - nut v. salt only****
# 14 - salt v. combo only


#########################
####### NUTRIENT ########
#########################


# how many are different with combo-nutrient only?
Shared_compareNut <- Compare_shared(SigOverlap$InCommonAll$Gene,
                                    my_dataSig$Nut_vs_Combo_results$Gene)
lapply(Shared_compareNut, function(x) {length(x)})
#3225 shared same (only 87 more than when considering salt differences); 1807 shared different


##### Of the ones that are different, which are: higher in combo, lower in combo, opposite in combo?

NutCombo_df <- merge(my_dataSig$condition_nutDE_results, my_dataSig$condition_comboDE_results, by="Gene")
length(NutCombo_df$Gene) #9029

NutComboDiff2 <- NutCombo_df[which(NutCombo_df$Gene %in% Shared_compareNut$Shared_different),]
length(NutComboDiff2$Gene) #to check

NutComboDiff2_compare <- Treat_compare(NutComboDiff2, "log2FoldChange.x", "log2FoldChange.y", "Nutrient", "Combo")
lapply(NutComboDiff2_compare, function(x) {length(x$Gene)})
# Up in Nutrient, reduced in combo: 605 (increased in combo: 99)
# Down in Nutrient, reduced in combo: 870 (incrased in combo: 94)
# Opposite direction (Up in Nutrient: 68, Down in Nutrient 71)

