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
###### DIRECTIONS #######
#########################

my_dataSigUp <- lapply(my_dataSig, MoreCritNum, column=3, critNum=0)
my_dataSigDown <- lapply(my_dataSig, LessCritNum, column=3, critNum=0)

SigOverlapUp <- GeneSets(my_dataSigUp$condition_comboDE_results[1], my_dataSigUp$condition_nutDE_results[1], 
                         my_dataSigUp$condition_saltDE_results[1])
lapply(SigOverlapUp, function(x) {length(x$Gene)})

SigOverlapDown <- GeneSets(my_dataSigDown$condition_comboDE_results[1], my_dataSigDown$condition_nutDE_results[1], 
                           my_dataSigDown$condition_saltDE_results[1])
lapply(SigOverlapDown, function(x) {length(x$Gene)})

seq(1,7,by=1)

# how many in each category are in different directions
SigOverlap_SameDir <- lapply(seq(1,7,by=1), function(x) {union(SigOverlapUp[[x]]$Gene, SigOverlapDown[[x]]$Gene)})
names(SigOverlap_SameDir) <- names(SigOverlap)
lapply(SigOverlap_SameDir, function(x) {length(x)})

SigOverlap_DiffDir <- lapply(names(SigOverlap), function(x) {setdiff(SigOverlap[[x]]$Gene, SigOverlap_SameDir[[x]])})
names(SigOverlap_DiffDir) <- names(SigOverlap)
lapply(SigOverlap_DiffDir, function(x) {length(x)})

#########################
# NUTRIENT-COMBO OVERALL #
#########################

# not factoring in salt significance

NutCombo_Cats <- Treatment_compare(my_dataSig$condition_nutDE_results,
                                   my_dataSig$condition_comboDE_results,
                                   my_dataSig$Nut_vs_Combo_results,
                                   "Nutrient",
                                   "Combo")
lapply(NutCombo_Cats, function(x) {length(x)})
# Nutrient only: 10,491
# Combo-Nutrient same: 5602
# reduced in Combo: 3051
# increased in Combo: 199
# diff direction in Combo: 177

# factoring Salt significance
NutCombo_Cats_noSalt <- Treatment_compare(my_dataSig$condition_nutDE_results,
                                          my_dataSig$condition_comboDE_results[which(my_dataSig$condition_comboDE_results$Gene %in%
                                                                                       SigOverlap$`my_dataSig$condition_comboDE_results[1]my_dataSig$condition_nutDE_results[1]Only`$Gene),],
                                          my_dataSig$Nut_vs_Combo_results,
                                          "Nutrient",
                                          "Combo")
lapply(NutCombo_Cats_noSalt, function(x) {length(x)})
# Combo-Nutrient same: 2377
# Combo reduced: 1576
# Combo increased: 6
# Diff direction: 38


#########################
# SALT-COMBO OVERALL #
#########################

# not factoring in nutrient significance

SaltCombo_Cats <- Treatment_compare(my_dataSig$condition_saltDE_results,
                                   my_dataSig$condition_comboDE_results,
                                   my_dataSig$Salt_vs_Combo_results,
                                   "Salt",
                                   "Combo")
lapply(SaltCombo_Cats, function(x) {length(x)})
# Salt only: 3,873
# Combo-Salt same: 5975
# reduced in Combo: 126
# increased in Combo: 126
# diff direction in Combo: 17

# factoring Nut significance
SaltCombo_Cats_noNut <- Treatment_compare(my_dataSig$condition_saltDE_results,
                                          my_dataSig$condition_comboDE_results[which(my_dataSig$condition_comboDE_results$Gene %in%
                                                                                       SigOverlap$`my_dataSig$condition_comboDE_results[1]my_dataSig$condition_saltDE_results[1]Only`$Gene),],
                                          my_dataSig$Salt_vs_Combo_results,
                                          "Salt",
                                          "Combo")
lapply(SaltCombo_Cats_noNut, function(x) {length(x)})
# Combo-Salt same: 1142
# reduced in Combo: 62
# increased in Combo: 1
# diff direction in Combo: 7



#########################
##### NUT-SALT ONLY #####
#########################

# not significant in combo
Nut_Salt_cats <- Treatment_compare(my_dataSig$condition_nutDE_results,
                                   my_dataSig$condition_saltDE_results[which(my_dataSig$condition_saltDE_results$Gene %in%
                                                                               SigOverlap$`my_dataSig$condition_nutDE_results[1]my_dataSig$condition_saltDE_results[1]Only`$Gene),],
                                   my_dataSig$Nut_vs_Salt_results,
                                   "Nutrient",
                                   "Salt")
lapply(Nut_Salt_cats, function(x) {length(x)})
# 1115 shared same, "1-treatment-specific"
# reduced in salt: 107
# increased in salt: 27
# Diff directions: 606 - antagonistic

# 1115 plus the 134 (different but same direction)= 1,249 "One treatment-specific"



#########################
##### ALL 3 SHARED ######
#########################

lapply(SigOverlap, function(x) {length(x$Gene)})
# 5032 are in common among all 3
# 1855 are in nutrient + salt, but not combo

# how many of the shared responses are antagonistic between Nutrient and Salt?
All_sharedNutSalt <- NutSalt_df[which(NutSalt_df$Gene %in% SigOverlap$InCommonAll$Gene),]
length(All_sharedNutSalt$Gene) #to check

All_shared_NutSaltcompare <- Treat_compare(All_sharedNutSalt, "log2FoldChange.x", "log2FoldChange.y", "Nutrient", "Salt")
lapply(All_shared_NutSaltcompare, function(x) {length(x$Gene)})
# reduced in Salt (1605 Up, 1730 Down)
# increased in Salt (753 Up, 795 Down)
# opposite (76 Up in Nutrient, Down Salt; 73 Up in Salt, Down Nutrient) = 149


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

# 511 - nut v. salt only**** intermediate?

# 14 - salt v. combo only



AntagonisticDiffs <- lapply(InCommonDiffs, function(x) {intersect(x, union(All_shared_NutSaltcompare$Up_Nutrient_OppositeDir_Salt$Gene,
                                                                           All_shared_NutSaltcompare$Down_Nutrient_OppositeDir_Salt$Gene))})
lapply(AntagonisticDiffs, function(x) {length(x)})
# 24 (In common all)
# 121 (Nut v. Combo, Nut vs Salt)
# 3 (Nut v. Salt, Salt v. Combo)
# 1 (Nut v. Salt only)

#########################
# ALL 3 SHARED - NUTRIENT #
#########################

# Nut-Combo differences:
Nut_Combo_sharedall <- Compare_shared(SigOverlap$InCommonAll$Gene, my_dataSig$Nut_vs_Combo_results$Gene)
lapply(Nut_Combo_sharedall, function(x) {length(x)})
# shared same: 3225 (-3138 = 87), shared different 1807

All_shared_NutComboDiff <- NutCombo_df[which(NutCombo_df$Gene %in% Nut_Combo_sharedall$Shared_different),]
length(All_shared_NutComboDiff$Gene) #to check

All_shared_NutComboDiff_compare <- Treat_compare(All_shared_NutComboDiff, "log2FoldChange.x", "log2FoldChange.y", "Nutrient", "Combo")
lapply(All_shared_NutComboDiff_compare, function(x) {length(x$Gene)})

# reduced in Combo (605 Up in Nutrient, 870 Down in Nutrient) = 1,475
  # n.s. different between Salt- Combo (572 Up Nutrient, 824 Down in Nutrient) = 1,396
      # - remaining 79: 

# increased in Combo (99 Up in Nutrient, 94 Down in Nutrient) = 193
  # n.s. different between Salt-Combo (92 Up Nutrient, 85 Down in Nutrient) = 177
      # - remaining 16:

# Opposite Direction (68 Up in Nutrient, 71 Down in Nutrient) = 139
  # n.s. different between Salt-Combo: (58 Up Nutrient, 64 Down in Nutrient) = 122
    # - remaining 17:

lapply(All_shared_NutComboDiff_compare, function(x) {length(setdiff(x$Gene, my_dataSig$Salt_vs_Combo_results$Gene))})

# compare with similar results from salt

#########################
# ALL 3 SHARED - SALT #
#########################

# Salt-Combo differences:
Salt_Combo_sharedall <- Compare_shared(SigOverlap$InCommonAll$Gene, my_dataSig$Salt_vs_Combo_results$Gene)
lapply(Salt_Combo_sharedall, function(x) {length(x)})
# shared same: 4833 (-3138 = 1,695), shared different 199

Salt_Combo_sharedall_compare <- Treat_compare(SaltCombo_df[which(SaltCombo_df$Gene %in% Salt_Combo_sharedall$Shared_different),], 
                                             "log2FoldChange.x", "log2FoldChange.y", "Salt", "Combo")
lapply(Salt_Combo_sharedall_compare, function(x) {length(x$Gene)})

# reduced in Combo (25 Up in Salt, 39 Down in Salt) = 64
  # n.s. different between nut-combo: (13 Up in Salt, 24 Down in Salt): 37
    # remaining 27:
# increased in Combo (55 Up in Salt, 70 Down in Salt) = 125
  # n.s. different between nut-combo: (24 Up in Salt, 22 Down in Salt): 46
    # remaining 79:
# Opposite Direction (2 Up in Salt, 8 Down in Salt) = 10
  # n.s. different from nut-combo: (0 Up in Salt, 3 Down in Salt): 3
    # remaining 7:
lapply(Salt_Combo_sharedall_compare, function(x) {length(setdiff(x$Gene, my_dataSig$Nut_vs_Combo_results$Gene))})






