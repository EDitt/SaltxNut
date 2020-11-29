#########################
######## SETUP ##########
#########################

library(dplyr)
library(ggplot2)
source("Functions.R")

library(viridis)
library(RColorBrewer)
library(ggsci)
library(scales)

load("ResultsFiles/Overlap.RData")
# DEData: All data for all pairwise contrasts
# SigOverlap: Overlap between DE genes for 3 treatments
# SigDiffOverlap: Overlap between sig. pairwise differences: Combo-Nut, Combo-Salt, Nut-Salt

lapply(SigOverlap, function(x) {length(x$Gene)})
lapply(SigDiffOverlap, function(x) {length(x$Gene)})

my_dataSig <- lapply(DEData, SigDEdf, PvaluesCol=7, CritP=0.05) #list that contains only the genes that are significant

# total number of DE genes
TotNumDEgenes <- Reduce("+", lapply(SigOverlap, function(x) {length(x$Gene)})) #23,789

#########################
###### ANTAGONISTIC #####
#########################

# how many shared by Salt/Nut are in different directions?
my_dataSigUp <- lapply(my_dataSig, MoreCritNum, column=3, critNum=0)
my_dataSigDown <- lapply(my_dataSig, LessCritNum, column=3, critNum=0)
NutUpSaltDown <- intersect(my_dataSigUp$condition_nutDE_results$Gene, my_dataSigDown$condition_saltDE_results$Gene)
length(NutUpSaltDown) #389
NutDownSaltUp <- intersect(my_dataSigDown$condition_nutDE_results$Gene, my_dataSigUp$condition_saltDE_results$Gene)
length(NutDownSaltUp) #366
NutSaltOpposite <- union(NutUpSaltDown, NutDownSaltUp)
length(NutSaltOpposite) #755

#########################
### COMBO INTERMEDIATE ##
#########################

# DE and different among all 3 stressors, Combo is intermediate

# combine dataframes - log2FoldChange + Pvalues
my_columnsPL <- lapply(DEData, function(x) x[c(3,4,7)])
All_PLvalues <- do.call("cbind", my_columnsPL)
All_PLvalues <- cbind(DEData[[1]]["Gene"], All_PLvalues)

# different among all 3 stressors:
Sig_Alldiff <- intersect(SigOverlap$InCommonAll$Gene, SigDiffOverlap$InCommonAll$Gene)
length(Sig_Alldiff) #97

Sig_Alldiff_values <- All_PLvalues[which(All_PLvalues$Gene %in% Sig_Alldiff &
                                           !All_PLvalues$Gene %in% NutSaltOpposite),] #N=73

Combo_int1 <- subset(Sig_Alldiff_values, condition_comboDE_results.log2FoldChange > condition_nutDE_results.log2FoldChange &
                       condition_comboDE_results.log2FoldChange < condition_saltDE_results.log2FoldChange)
length(Combo_int1$Gene) #43

Combo_int2 <- subset(Sig_Alldiff_values, condition_comboDE_results.log2FoldChange < condition_nutDE_results.log2FoldChange &
                       condition_comboDE_results.log2FoldChange > condition_saltDE_results.log2FoldChange)
length(Combo_int2$Gene) #27

Combo_Intermed <- union(Combo_int1$Gene, Combo_int2$Gene)


#########################
#### COMBO ADDITIVE #####
#########################

Sig_allCombodiff <- intersect(SigOverlap$InCommonAll$Gene, intersect(my_dataSig$Nut_vs_Combo_results$Gene, my_dataSig$Salt_vs_Combo_results$Gene))
length(Sig_allCombodiff) #N=112
length(setdiff(Sig_allCombodiff, NutSaltOpposite)) #88

# DE in all three stressors and Combo different with Nut and Salt, and magnitude is larger than both
Sig_AllCombodiff_values <- All_PLvalues[which(All_PLvalues$Gene %in% Sig_allCombodiff &
                                           !All_PLvalues$Gene %in% NutSaltOpposite),] #N=88

Combo_additive <- subset(Sig_AllCombodiff_values, abs(condition_comboDE_results.log2FoldChange) > abs(condition_nutDE_results.log2FoldChange) &
                           abs(condition_comboDE_results.log2FoldChange) > abs(condition_saltDE_results.log2FoldChange))
length(Combo_additive$Gene) #12

# are these the same between Salt-Nut?
length(setdiff(Combo_additive$Gene, my_dataSig$Nut_vs_Salt_results$Gene)) #9

#########################
# NUTRIENT-COMBO OVERALL #
#########################

# not factoring in salt significance

NutCombo_Cats <- Treatment_compare(my_dataSig$condition_nutDE_results,
                                   my_dataSig$condition_comboDE_results,
                                   my_dataSig$Nut_vs_Combo_results,
                                   "Nutrient",
                                   "Combo")
NutCatNums <- lapply(NutCombo_Cats, function(x) {length(x)})
# Nutrient only: 10,491
# Combo-Nutrient same: 5602
# reduced in Combo: 3051
# increased in Combo: 199
# diff direction in Combo: 177

# how many of these are not significant in Salt? (e.g. Nut-specific)
NutSpecific_CatNums <- lapply(NutCombo_Cats, function(x) {setdiff(x, my_dataSig$condition_saltDE_results$Gene)})
lapply(NutSpecific_CatNums, function(x) {length(x)})
# Nutrient only: 8,636 - "Nutrient-specific, conditional"
# Combo-Nutrient same: 2,377 - "Nutrient-specific, ***Unconditional"
# reduced in Combo: 1576 - "Nutrient-specific, conditional"
# increased in Combo: 6
# diff direction in Combo: 38

NutUnSpecific_CatNums <- lapply(NutCombo_Cats, function(x) {intersect(x, my_dataSig$condition_saltDE_results$Gene)})
lapply(NutUnSpecific_CatNums, function(x) {length(x)})
# Nutrient only: 1855: 
  # "1-stress Specific" OR antagonistic
  # need to see how many are in opposite directions in salt/nut
# Combo-Nutrient Same: 3225: "General Stress Response"
# Combo reduced: 1475
# Combo increased: 193
# Combo diff direction: 139

# For General Stress Response:
lapply(SigDiffOverlap, function(x) {length(intersect(x$Gene, NutUnSpecific_CatNums$Combo_Nutrient_Same))})
# N=73: Nut vs. Salt, Salt vs. Combo; N=511: Nut vs. Salt; N=14: Salt vs Combo
# 3225 - (73 + 14) = 3,138 for "General Stress Response"

# Combo reduced, Not Nutrient-specific:
# how many are the same between Salt-Combo?
length(setdiff(NutUnSpecific_CatNums$Combo_reduced, my_dataSig$Salt_vs_Combo_results$Gene)) #1396
# how many are also different between Salt-Combo?
length(intersect(NutUnSpecific_CatNums$Combo_reduced, my_dataSig$Salt_vs_Combo_results$Gene)) #79
  # how many are antagonistic?
Combo_red_Antag <- intersect(NutUnSpecific_CatNums$Combo_reduced, NutSaltOpposite) #7
  # how many are intermediate between Salt-Nut?
Combo_red_Inter <- intersect(NutUnSpecific_CatNums$Combo_reduced, Combo_Intermed) #66
  # remaining
length(setdiff(intersect(NutUnSpecific_CatNums$Combo_reduced, my_dataSig$Salt_vs_Combo_results$Gene),
               union(Combo_red_Antag, Combo_red_Inter)))



#### Unchanged:
Nut_Unchanged <- list(NutSpecific_CatNums$Combo_Nutrient_Same, setdiff(NutUnSpecific_CatNums$Combo_Nutrient_Same, my_dataSig$Salt_vs_Combo_results$Gene),
                      intersect(NutUnSpecific_CatNums$Combo_Nutrient_Same, my_dataSig$Salt_vs_Combo_results$Gene))
names(Nut_Unchanged) <- c("Unconditional_StressSpecific", "Unconditional_GeneralStress", "Other") ###come up with diff. name for "other"
lapply(Nut_Unchanged, function(x) {length(x)})

#### No Longer Significant:
Nut_ComboNS <-list(NutSpecific_CatNums$Nutrient_only, intersect(NutUnSpecific_CatNums$Nutrient_only, NutSaltOpposite), 
                   setdiff(NutUnSpecific_CatNums$Nutrient_only, NutSaltOpposite))
names(Nut_ComboNS) <- c("Conditional_StressSpecific", "Antagonistic", "One-stress_specific")
lapply(Nut_ComboNS, function(x) {length(x)})
#134 of 1-stress specific (N=1249) are diff between Nut-Salt

#### Decreased Magnitude:
Nut_Combo_dec <- list(NutSpecific_CatNums$Combo_reduced, setdiff(NutUnSpecific_CatNums$Combo_reduced, my_dataSig$Salt_vs_Combo_results$Gene),
                      Combo_red_Inter, Combo_red_Antag, 
                      setdiff(intersect(NutUnSpecific_CatNums$Combo_reduced, my_dataSig$Salt_vs_Combo_results$Gene),
                              union(Combo_red_Antag, Combo_red_Inter)))
names(Nut_Combo_dec) <- c("Conditional_StressSpecific", "Salt_dominant", "Intermediate", "Antagonistic", "other")
lapply(Nut_Combo_dec, function(x) {length(x)})

#########################
### SALT-COMBO OVERALL ###
#########################

# not factoring in nutrient significance

SaltCombo_Cats <- Treatment_compare(my_dataSig$condition_saltDE_results,
                                   my_dataSig$condition_comboDE_results,
                                   my_dataSig$Salt_vs_Combo_results,
                                   "Salt",
                                   "Combo")
SaltCatNums <- lapply(SaltCombo_Cats, function(x) {length(x)})
# Salt only: 3,873
# Combo-Salt same: 5975
# reduced in Combo: 126
# increased in Combo: 126
# diff direction in Combo: 17

# how many of these are not significant in Nut?
SaltSpecific_CatNums <- lapply(SaltCombo_Cats, function(x) {setdiff(x, my_dataSig$condition_nutDE_results$Gene)})
lapply(SaltSpecific_CatNums, function(x) {length(x)})
# Salt only: 2,018
# Combo-Salt same: 1,142 - "Salt-specific, Unconditional"
# reduced in Combo: 62
# increased in Combo: 1
# diff direction in Combo: 7

SaltUnSpecific_CatNums <- lapply(SaltCombo_Cats, function(x) {intersect(x, my_dataSig$condition_nutDE_results$Gene)})
lapply(SaltUnSpecific_CatNums, function(x) {length(x)})

Salt_Unchanged <- list(SaltSpecific_CatNums$Combo_Salt_Same, setdiff(SaltUnSpecific_CatNums$Combo_Salt_Same, my_dataSig$Nut_vs_Combo_results$Gene),
                      intersect(SaltUnSpecific_CatNums$Combo_Salt_Same, my_dataSig$Nut_vs_Combo_results$Gene))
names(Salt_Unchanged) <- c("Unconditional_StressSpecific", "Unconditional_GeneralStress", "Other") ###come up with diff. name for "other"
lapply(Salt_Unchanged, function(x) {length(x)})

#########################
#### CREATE DATAFRAME ###
#########################

#### Save:
#save(NutCombo_Cats, SaltCombo_Cats, file="ResultsFiles/StressAdditionCats.RData")

In_Combo <- c("Not_Significant", "Un_changed", "Lower_Magnitude", "Higher_Magnitude", "Opposite_direction")

Nut_df <- data.frame(unlist(NutCatNums, recursive = TRUE, use.names = TRUE))
colnames(Nut_df) <- c("Number")
Nut_df$Stress <- c("Nutrient")
Nut_df_labels <- cbind(Nut_df, In_Combo)
Nut_total <- sum(Nut_df_labels$Number)
Nut_df_labels$Prop2 <- Nut_df_labels$Number / Nut_total

Salt_df <- data.frame(unlist(SaltCatNums, recursive = TRUE, use.names = TRUE))
colnames(Salt_df) <- c("Number")
Salt_df$Stress <- c("Salt")
Salt_df_labels <- cbind(Salt_df, In_Combo)
Salt_total <- sum(Salt_df_labels$Number)
Salt_df_labels$Prop2 <- Salt_df_labels$Number / Salt_total

All_df <- rbind(Nut_df_labels, Salt_df_labels, make.row.names = FALSE)
str(All_df)
All_df$Stress <- factor(All_df$Stress)
All_df$In_Combo <- factor(All_df$In_Combo, 
                          levels = c("Opposite_direction", "Higher_Magnitude",
                                     "Lower_Magnitude", "Not_Significant",
                                     "Un_changed"))

# add Proportion
### should I make it a proportion of all genes in each single stress (so bars sum to 1 for both)?
All_df$Proportion <- All_df$Number / TotNumDEgenes

# change order of levels

levels(All_df$In_Combo)

#########################
####### BAR CHART #######
#########################

viridis(5)
viridis(8)
brewer.pal(8, "Dark2")
show_col(pal_jco("default")(10))

# proportion out of all DE genes
p <- ggplot(All_df, aes(fill=In_Combo, y=Proportion, x=Stress))

# proportion out of all DE genes for single stress
p2 <- ggplot(All_df, aes(fill=In_Combo, y=Prop2, x=Stress))

# Use the pal_jco palette, except changed the blue color to one from viridis palette
# make red the darker shade
# - I like this one*****
# Switch order of pink and yellow colors??? 
p2 + geom_bar(position="stack", stat="identity", color="black", size=0.4) + ylim(0,1) +
  ylab("Proportion of DE Genes") + xlab("Single Stress") +
  scale_fill_manual(values=c("#3B528BFF", "#868686FF", "#EFC000FF", "#A73030FF",  "#8F7700FF"),
                    breaks=c("Un_changed", "Not_Significant", "Lower_Magnitude",
                             "Higher_Magnitude", "Opposite_direction"),
                    labels=c("Un-changed", "No Longer Differentially Expressed",
                             "Decreased Magnitude", "Increased Magnitude",
                             "Opposite Direction"),
                    name="Addition of Second Stress") +
  theme_bw(base_size = 14)


######################################
# experimenting with colors-
p + geom_bar(position="stack", stat="identity") + ylim(0,1) +
  ylab("Proportion of DE Genes") + xlab("Single Stress") +
  scale_fill_manual(values=c(viridis(5)),
                    breaks=c("Un_changed", "Not_Significant", "Lower_Magnitude",
                             "Higher_Magnitude", "Opposite_direction"),
                    labels=c("Un-changed", "No Longer Differentially Expressed",
                             "Decreased Magnitude", "Increased Magnitude",
                             "Opposite Direction"),
                    name="Effect in Combination") +
  theme_bw()


p + geom_bar(position="stack", stat="identity") + ylim(0,1) +
  ylab("Proportion of DE Genes") + xlab("Single Stress") +
  scale_fill_manual(values=c("#0073C2FF", "#868686FF", "#EFC000FF", "#CD534CFF", "#8F7700FF"),
                    breaks=c("Un_changed", "Not_Significant", "Lower_Magnitude",
                             "Higher_Magnitude", "Opposite_direction"),
                    labels=c("Un-changed", "No Longer Differentially Expressed",
                             "Decreased Magnitude", "Increased Magnitude",
                             "Opposite Direction"),
                    name="Effect in Combination") +
  theme_bw()



p + geom_bar(position="stack", stat="identity") + ylim(0,1) +
  ylab("Proportion of DE Genes") + xlab("Single Stress") +
  scale_fill_jco(
                    breaks=c("Un_changed", "Not_Significant", "Lower_Magnitude",
                             "Higher_Magnitude", "Opposite_direction"),
                    labels=c("Un-changed", "No Longer Differentially Expressed",
                             "Decreased Magnitude", "Increased Magnitude",
                             "Opposite Direction"),
                    name="Effect in Combination") +
  theme_bw()

