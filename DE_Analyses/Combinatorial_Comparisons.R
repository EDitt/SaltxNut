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


#########################
#### CREATE DATAFRAME ###
#########################

#### Save:
save(NutCombo_Cats, SaltCombo_Cats, file="ResultsFiles/StressAdditionCats.RData")

In_Combo <- c("Not_Significant", "Un_changed", "Lower_Magnitude", "Higher_Magnitude", "Opposite_direction")

Nut_df <- df_from_List(NutCatNums, In_Combo, "Nutrient")
Salt_df <- df_from_List(SaltCatNums, In_Combo, "Salt")

All_df <- rbind(Nut_df, Salt_df, make.row.names = FALSE)
str(All_df)
All_df$Stress <- factor(All_df$Stress)
All_df$labels <- factor(All_df$labels, 
                          levels = c("Opposite_direction", "Higher_Magnitude",
                                     "Lower_Magnitude", "Not_Significant",
                                     "Un_changed"))

# add Proportion
### should I make it a proportion of all genes in each single stress (so bars sum to 1 for both)?
All_df$Proportion <- All_df$Number / TotNumDEgenes

#########################
####### BAR CHART #######
#########################

viridis(5)
viridis(8)
brewer.pal(8, "Dark2")
show_col(pal_jco("default")(10))

# proportion out of all DE genes
#p <- ggplot(All_df, aes(fill=labels, y=Proportion, x=Stress))

# proportion out of all DE genes for single stress
p2 <- ggplot(All_df, aes(fill=labels, y=Prop, x=Stress))

# Use the pal_jco palette, except changed the blue color to one from viridis palette
# make red the darker shade
# - I like this one*****
# Switch order of pink and yellow colors??? 
p2 + geom_bar(position="stack", stat="identity", color="black", size=0.4) + ylim(0,1) +
  ylab("Proportion of DE Genes") + xlab("Single Stress") +
  scale_fill_manual(values=c( "#8F7700FF", "#A73030FF", "#EFC000FF", "#868686FF", "#3B528BFF"),
                    #values=c("#3B528BFF", "#868686FF", "#EFC000FF", "#A73030FF",  "#8F7700FF"),
                    breaks=c("Un_changed", "Not_Significant", "Lower_Magnitude",
                             "Higher_Magnitude", "Opposite_direction"),
                    labels=c("Un-changed", "No Longer Differentially Expressed",
                             "Decreased Magnitude", "Increased Magnitude",
                             "Opposite Direction"),
                    name="Addition of Second Stress") +
  theme_bw(base_size = 14)

#########################
### NUT/SALT SPECIFIC ###
#########################

# how many of the Nutrient categories are not significant in Salt? (e.g. Nut-specific)
NutSpecific_Cats <- lapply(NutCombo_Cats, function(x) {setdiff(x, my_dataSig$condition_saltDE_results$Gene)})
NutSpecific_CatNums <- lapply(NutSpecific_Cats, function(x) {length(x)})

# how many of these are not significant in Nut?
SaltSpecific_Cats <- lapply(SaltCombo_Cats, function(x) {setdiff(x, my_dataSig$condition_nutDE_results$Gene)})
SaltSpecific_CatNums <- lapply(SaltSpecific_Cats, function(x) {length(x)})

labels <- c("Not_DE", "Unchanged", "Decreased_mag", "Increased_mag", "Opposite_dir")
NutSpecific_df <- df_from_List(NutSpecific_CatNums, labels, "Nutrient_specific")
SaltSpecific_df <- df_from_List(SaltSpecific_CatNums, labels, "Salt_specific")

#########################
###### UN-SPECIFIC ######
#########################

lapply(Combo_patterns, function(x) {length(x)})
# shared with both: 606 cancelled, 1249 1-stress-specific, 
# 2627 all shared, 12 additive-combo higher, 6 "subtractive" combo lower


NutUnSpecific_Cats <- lapply(NutCombo_Cats, function(x) {intersect(x, my_dataSig$condition_saltDE_results$Gene)})
NutUnSpecific_CatNums <- lapply(NutUnSpecific_Cats, function(x) {length(x)}) 
# Nutrient only: 1855, Combo-Nutrient_same: 3225, Combo_reduced: 1475, combo_increased: 193, diff direction: 139

SaltUnSpecific_Cats <- lapply(SaltCombo_Cats, function(x) {intersect(x, my_dataSig$condition_nutDE_results$Gene)})
SaltUnSpecific_CatNums <- lapply(SaltUnSpecific_Cats, function(x) {length(x)}) # Salt only: 1855
# Salt only: 1855, Combo-Salt_same: 4833, Combo_reduced: 64, combo_increased: 125, diff direction:10


### Overlap
### Same:
lapply(Combo_patterns, function(x) {length(intersect(x, NutUnSpecific_Cats$Combo_Nutrient_Same))})
# 2627-all shared, 4-antagonistic_nut-dominant
lapply(SaltUnSpecific_Cats, function(x) {length(intersect(x, NutUnSpecific_Cats$Combo_Nutrient_Same))})
# 3138- combo-salt same: 2627 all shared + 511 (Nut v. Salt only)
# Combo reduced (N=37), Combo increased (N=47), Diff direction (N=3)



lapply(Combo_patterns, function(x) {length(intersect(x, NutUnSpecific_Cats$Combo_reduced))})
# Out of 1475: antagonistic_NutDir (N=7), Intermediate_Nut-higher (N=66), subtractive (6)


lapply(Combo_patterns, function(x) {length(intersect(x, SaltUnSpecific_Cats$Combo_Salt_Same))})
# 2627-all shared, 121-antagonistic_salt-dominant



lapply(SaltUnSpecific_Cats, function(x) {length(intersect(x, NutUnSpecific_Cats$Combo_reduced))})
# 1396- combo-salt same, 6 combo reduced (subtractive), 66 combo increased (Intermediate), 7 diff direction (antagonistic)

lapply(SaltUnSpecific_Cats, function(x) {length(intersect(x, NutUnSpecific_Cats$Nutrient_only))})
# 1855- salt-only

#########################
####### DATAFRAME #######
#########################

Categories <- c("Un_changed", "Not_DE", "Mag_decreased")
Stress <- c("Nutrient", "Salt")


# for NS in Combo:
Categories <- c("Cancelled", "One_stress_specific", "Stress_specific_conditional")
Stress <- c("Nutrient", "Salt")

# for Un-changed:
Categories <- c("Shared_response", "Stress_specific_Unconditional")

# for Mag. decreased
Categories <- c("Stress_specific_conditional", "Antagonistic", "")

######### SCRATCH

#########################
### NUT-COMBO SEPARATE ##
#########################



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
## SALT-COMBO SEPARATE ##
#########################



SaltUnSpecific_CatNums <- lapply(SaltCombo_Cats, function(x) {intersect(x, my_dataSig$condition_nutDE_results$Gene)})
lapply(SaltUnSpecific_CatNums, function(x) {length(x)})

Salt_Unchanged <- list(SaltSpecific_CatNums$Combo_Salt_Same, setdiff(SaltUnSpecific_CatNums$Combo_Salt_Same, my_dataSig$Nut_vs_Combo_results$Gene),
                       intersect(SaltUnSpecific_CatNums$Combo_Salt_Same, my_dataSig$Nut_vs_Combo_results$Gene))
names(Salt_Unchanged) <- c("Unconditional_StressSpecific", "Unconditional_GeneralStress", "Other") ###come up with diff. name for "other"
lapply(Salt_Unchanged, function(x) {length(x)})

###################################### SCRATCH
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



p2 + geom_bar(position="stack", stat="identity") + ylim(0,1) +
  ylab("Proportion of DE Genes") + xlab("Single Stress") +
  scale_fill_jco(
                    breaks=c("Un_changed", "Not_Significant", "Lower_Magnitude",
                             "Higher_Magnitude", "Opposite_direction"),
                    labels=c("Un-changed", "No Longer Differentially Expressed",
                             "Decreased Magnitude", "Increased Magnitude",
                             "Opposite Direction"),
                    name="Effect in Combination") +
  theme_bw()

