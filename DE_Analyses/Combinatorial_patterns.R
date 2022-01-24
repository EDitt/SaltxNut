#########################
######## SETUP ##########
#########################

library(dplyr)
library(ggplot2)
source("Functions.R")

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

my_dataSigUp <- lapply(my_dataSig, MoreCritNum, column=3, critNum=0)
my_dataSigDown <- lapply(my_dataSig, LessCritNum, column=3, critNum=0)

load("ResultsFiles/StressAdditionCats.RData")
lapply(NutCombo_Cats, function(x) {length(x)})
lapply(SaltCombo_Cats, function(x) {length(x)})

#########################
###### ANTAGONISTIC #####
#########################

# Significant in Salt & Nutrient (possibly Combo)

# how many shared by Salt/Nut are in different directions?
NutUpSaltDown <- intersect(my_dataSigUp$condition_nutDE_results$Gene, my_dataSigDown$condition_saltDE_results$Gene)
length(NutUpSaltDown) #389
NutDownSaltUp <- intersect(my_dataSigDown$condition_nutDE_results$Gene, my_dataSigUp$condition_saltDE_results$Gene)
length(NutDownSaltUp) #366
NutSaltOpposite <- union(NutUpSaltDown, NutDownSaltUp)
length(NutSaltOpposite) #755


# how many antagonistic cancel each other out in Combo?
lapply(SigOverlap, function(x) {length(intersect(x$Gene, NutSaltOpposite))})
#606 Nut & Salt only; #149 in common all

Cancelled <- intersect(SigOverlap$`my_dataSig$condition_nutDE_results[1]my_dataSig$condition_saltDE_results[1]Only`$Gene,
                       NutSaltOpposite)
length(Cancelled) #606

# what do the remaining ones look like?
Antagonistic_cats <- lapply(SigDiffOverlap, function(x) {intersect(x$Gene, intersect(NutSaltOpposite, SigOverlap$InCommonAll$Gene))})
lapply(Antagonistic_cats, function(x) {length(x)})
#24 in common all; 121 Nut different from other 2; 3 Salt different from other 2; 1 Nut vs. Salt only

# the 1 that's only diff Nut-Salt??
All_PLvalues[which(All_PLvalues$Gene %in% Antagonistic_cats$`my_dataSig$Nut_vs_Salt_results[1]Only`),]
# up in Combo, up in nutrient, down in salt <- call this Nut dominant

# opposite directions, no differences between combo-salt
Antagonistic_SaltDom <- Antagonistic_cats$`my_dataSig$Nut_vs_Combo_results[1]my_dataSig$Nut_vs_Salt_results[1]Only`
length(Antagonistic_SaltDom) #121

# opposite directions, no differences between combo-nut
Antagonistic_NutDom <- union(Antagonistic_cats$`my_dataSig$Nut_vs_Salt_results[1]my_dataSig$Salt_vs_Combo_results[1]Only`,
                             Antagonistic_cats$`my_dataSig$Nut_vs_Salt_results[1]Only`)
length(Antagonistic_NutDom) # 3 + 1 (Nut-Salt diff only)

# the remaining could be intermediate?
lapply(NutCombo_Cats, function(x) {length(intersect(x, Antagonistic_cats$InCommonAll))})
# 7 reduced in combo (relative to nutrient); 17 different direction (relative to nutrient)

lapply(SaltCombo_Cats, function(x) {length(intersect(x, Antagonistic_cats$InCommonAll))})
# 17 reduced in combo (relative to salt): 7 different direction (relative to salt)

Ant_AllDiffCats <- lapply(NutCombo_Cats, function(x) {intersect(x, Antagonistic_cats$InCommonAll)})

# reduced magnitude, same direction as nutrient:
Antagonistic_NutDirDom <- Ant_AllDiffCats$Combo_reduced

# reduced magnitude, same direction as salt:
Antagonistic_SaltDirDom <- Ant_AllDiffCats$Combo_Diff_direction

Antagonistic_List <- list(Cancelled, Antagonistic_SaltDom, Antagonistic_NutDom, Antagonistic_NutDirDom, Antagonistic_SaltDirDom)
names(Antagonistic_List) <- c("Cancelled", "Antagonistic_SaltDominant", "Antagonistic_NutDominant", "Antagonistic_NutrientDir", "Antagnoistic_SaltDir")
lapply(Antagonistic_List, function(x) {length(x)})

lapply(Antagonistic_List, function(x) {length(intersect(x, intersect(SigOverlap$InCommonAll$Gene, SigDiffOverlap$InCommonAll$Gene)))})

#########################
### COMBO INTERMEDIATE ##
#########################

# Significant in all 3 stresses

# Signifiant in all 3 & not in opposite directions between Salt/Nut
Combo_int_Salthigher <- subset(All_PLvalues[which(!All_PLvalues$Gene %in% NutSaltOpposite &
                                                    All_PLvalues$Gene %in% SigOverlap$InCommonAll$Gene),], 
                               condition_comboDE_results.log2FoldChange > condition_nutDE_results.log2FoldChange &
                                condition_comboDE_results.log2FoldChange < condition_saltDE_results.log2FoldChange)
length(Combo_int_Salthigher$Gene) #1274

Combo_int_Nuthigher <- subset(All_PLvalues[which(!All_PLvalues$Gene %in% NutSaltOpposite &
                                                   All_PLvalues$Gene %in% SigOverlap$InCommonAll$Gene),], 
                              condition_comboDE_results.log2FoldChange < condition_nutDE_results.log2FoldChange &
                                condition_comboDE_results.log2FoldChange > condition_saltDE_results.log2FoldChange)
length(Combo_int_Nuthigher$Gene) #1290

Combo_Intermed <- union(Combo_int_Salthigher$Gene, Combo_int_Nuthigher$Gene)
length(Combo_Intermed) #2564

lapply(SigDiffOverlap, function(x) {length(intersect(x$Gene, Combo_Intermed))})
# 70- In common all
# 924- Nut vs. Combo & Nut vs. Salt
# 42- Nut vs. Salt & Salt vs. Combo
# 500- Nut vs. Salt only

# to be significant, all pairwise diffs
Combo_Intermed_Sig <- intersect(Combo_Intermed, SigDiffOverlap$InCommonAll$Gene)
length(Combo_Intermed_Sig) #N=70

lapply(NutCombo_Cats, function(x) {length(intersect(x, Combo_Intermed_Sig))})
# 66 reduced in Combo (relative to Nutrient)
# 4 increased in Combo (relative to Nutrient)

lapply(SaltCombo_Cats, function(x) {length(intersect(x, Combo_Intermed))})
# 4 reduced in Combo (relative to Salt)
# 66 increased in Combo (relative to Salt)

Intermediate_list <- list(intersect(NutCombo_Cats$Combo_reduced, Combo_Intermed), intersect(NutCombo_Cats$Combo_increased, Combo_Intermed))
names(Intermediate_list) <- c("Intermediate_Nut_higher", "Intermediate_Salt_higher")


#########################
#### COMBO ADDITIVE #####
#########################

# Significant in all, Nut v. Combo & Salt v. Combo significant

Sig_allCombodiff <- intersect(SigOverlap$InCommonAll$Gene, intersect(my_dataSig$Nut_vs_Combo_results$Gene, my_dataSig$Salt_vs_Combo_results$Gene))
length(Sig_allCombodiff) #N=112 (2 pairwise diffs only + all pairwise diffs)
length(setdiff(Sig_allCombodiff, NutSaltOpposite)) #88

# DE in all three stressors and Combo different with Nut and Salt, and magnitude is larger than both
Sig_AllCombodiff_values <- All_PLvalues[which(All_PLvalues$Gene %in% Sig_allCombodiff &
                                                !All_PLvalues$Gene %in% NutSaltOpposite),] #N=88

Combo_additive <- subset(Sig_AllCombodiff_values, abs(condition_comboDE_results.log2FoldChange) > abs(condition_nutDE_results.log2FoldChange) &
                           abs(condition_comboDE_results.log2FoldChange) > abs(condition_saltDE_results.log2FoldChange))
length(Combo_additive$Gene) #12

# are these the same between Salt-Nut?
lapply(SigDiffOverlap, function(x) {length(intersect(x$Gene, Combo_additive$Gene))})
# In common all: N=3; Nut/Salt same: N=9

lapply(NutCombo_Cats, function(x) {length(intersect(x, Combo_additive$Gene))})
# 12 increased (relative to combo)
lapply(SaltCombo_Cats, function(x) {length(intersect(x, Combo_additive$Gene))})
# 12 increased (relative to combo)

#########################
# COMBO LOWER THAN BOTH #
#########################

# opposite of additive
Combo_subtract <- subset(Sig_AllCombodiff_values, abs(condition_comboDE_results.log2FoldChange) < abs(condition_nutDE_results.log2FoldChange) &
                           abs(condition_comboDE_results.log2FoldChange) < abs(condition_saltDE_results.log2FoldChange))
length(Combo_subtract$Gene) #6

lapply(SigDiffOverlap, function(x) {length(intersect(x$Gene, Combo_subtract$Gene))})

lapply(NutCombo_Cats, function(x) {length(intersect(x, Combo_subtract$Gene))})
# 6 reduced (relative to nutrient)

lapply(SaltCombo_Cats, function(x) {length(intersect(x, Combo_subtract$Gene))})
# 6 reduced (relative to salt)

Add_Subtract_list <- list(Combo_additive$Gene, Combo_subtract$Gene)
names(Add_Subtract_list) <- c("Additive_combo_higher", "Subtractive_combo_lower")

#########################
### COMBO SYNERGISTIC ###
#########################

# any for which combo is a diff. direction than other 2?
Combo_diffDirUp <- intersect(my_dataSigUp$condition_comboDE_results$Gene, intersect(my_dataSigDown$condition_nutDE_results$Gene,
                                                                             my_dataSigDown$condition_saltDE_results$Gene))
length(intersect(Combo_diffDirUp, Sig_allCombodiff)) #0

Combo_diffDirDown <- intersect(my_dataSigDown$condition_comboDE_results$Gene, intersect(my_dataSigUp$condition_nutDE_results$Gene,
                                                                                    my_dataSigUp$condition_saltDE_results$Gene))
length(Combo_diffDirDown) #1
length(intersect(Combo_diffDirDown, Sig_allCombodiff)) #0 (when looking at ones that are sig different with combo)
                           
lapply(NutCombo_Cats, function(x) {length(intersect(x, Combo_diffDirDown))})
lapply(SaltCombo_Cats, function(x) {length(intersect(x, Combo_diffDirDown))}) # combo-salt same???                      
All_PLvalues[which(All_PLvalues$Gene %in% Combo_diffDirDown),] # down in combo, up in nut/salt; but not sig. different between nut/salt or salt/combo

### NONE

#########################
### 1-STRESS SPECIFIC ###
#########################

# Significant in Nut/Salt only, not different between them

One_stress <- setdiff(SigOverlap$`my_dataSig$condition_nutDE_results[1]my_dataSig$condition_saltDE_results[1]Only`$Gene,
                      my_dataSig$Nut_vs_Salt_results$Gene)
length(One_stress) #1115 (134 are different but not "cancelled")

# Will call "One-stress-specific" anything that is significant in Salt & Nut but not in opposite directions
One_stress_specific <- setdiff(SigOverlap$`my_dataSig$condition_nutDE_results[1]my_dataSig$condition_saltDE_results[1]Only`$Gene,
                               Cancelled)
length(One_stress_specific) #1249

#########################
###### All SHARED #######
#########################

# Genes that are significant in all treatments and no pairwise differences
PairwiseDiffs <- union(my_dataSig$Nut_vs_Combo_results$Gene, union(my_dataSig$Nut_vs_Salt_results$Gene, my_dataSig$Salt_vs_Combo_results$Gene))
All_shared <- setdiff(SigOverlap$InCommonAll$Gene, PairwiseDiffs)

length(All_shared) #2627

Cats <- list(One_stress_specific, All_shared)
names(Cats) <- c("One_stress_specific", "All_shared")

### to compare to other datasets (e.g. Heliaphen):
write.table(SigOverlap$InCommonAll$Gene, file="/Volumes/GoogleDrive/My Drive/Active Projects/Transcriptomics_Exp/Results/ToCompare/AllStress.txt", quote=FALSE, row.names = FALSE) # N=5032
write.table(All_shared, file= "/Volumes/GoogleDrive/My Drive/Active Projects/Transcriptomics_Exp/Results/ToCompare/AllStress_Same.txt", quote=FALSE, row.names = FALSE) # N=2627

#########################
#### SAVE CATEGORIES ####
#########################

Combo_patterns <- c(Antagonistic_List, Intermediate_list, Add_Subtract_list, Cats)

lapply(Combo_patterns, function(x) {length(x)})


