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
my_dataSigUp <- lapply(my_dataSig, MoreCritNum, column=3, critNum=0)
my_dataSigDown <- lapply(my_dataSig, LessCritNum, column=3, critNum=0)

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

lapply(Combo_patterns, function(x) {length(x)}) ### where is the code for this list???
# shared with both: 606 cancelled, 1249 1-stress-specific, 
# 2627 all shared, 12 additive-combo higher, 6 "subtractive" combo lower


NutUnSpecific_Cats <- lapply(NutCombo_Cats, function(x) {intersect(x, my_dataSig$condition_saltDE_results$Gene)})
NutUnSpecific_CatNums <- lapply(NutUnSpecific_Cats, function(x) {length(x)}) 
# Nutrient only: 1855, Combo-Nutrient_same: 3225, Combo_reduced: 1475, combo_increased: 193, diff direction: 139

SaltUnSpecific_Cats <- lapply(SaltCombo_Cats, function(x) {intersect(x, my_dataSig$condition_nutDE_results$Gene)})
SaltUnSpecific_CatNums <- lapply(SaltUnSpecific_Cats, function(x) {length(x)}) # Salt only: 1855
# Salt only: 1855, Combo-Salt_same: 4833, Combo_reduced: 64, combo_increased: 125, diff direction:10

labels <- c("Not_DE", "Unchanged", "Decreased_mag", "Increased_mag", "Opposite_dir")
NutUnSpecific_df <- df_from_List(NutUnSpecific_CatNums, labels, "Nutrient_Unspecific")
SaltUnSpecific_df <- df_from_List(SaltUnSpecific_CatNums, labels, "Salt_Unspecific")

## which genes are in the same categories?
lapply(SaltUnSpecific_Cats, function(x) {length(intersect(x, NutUnSpecific_Cats$Combo_Diff_direction))})
# all 1855 of the "not DE" category
# 3138 of the "Unchanged"
# 6 of the reduced, 12 of the increased, 0 diff direction

#########################
#### SAVE GENE LISTS ####
#########################

save(NutSpecific_Cats, SaltSpecific_Cats, NutUnSpecific_Cats, SaltUnSpecific_Cats, 
     file="ResultsFiles/GeneSets/MultiStressCompare.RData")

#########################
### UN-SPECIFIC MATCH ###
#########################

lists <- 1:length(SaltUnSpecific_Cats)

# overlap across categories:
overlap_combinations <- lapply(SaltUnSpecific_Cats[lists], function(x) {lapply(NutUnSpecific_Cats[lists], 
                                                                               function (y) {length(intersect(x,y))})})
overlap_comb_df <- data.frame(unlist(overlap_combinations, recursive = TRUE, use.names = TRUE))
colnames(overlap_comb_df) <- c("Number")
overlap_comb_df$labels <- as.character(rownames(overlap_comb_df))
overlap_comb_df$Catinteraction_string <- strsplit(overlap_comb_df$labels, "[.]")
overlap_comb_df$Salt_cat <- as.factor(sapply(overlap_comb_df$Catinteraction_string, 
                                             "[", 1))
overlap_comb_df$Nut_cat <- as.factor(sapply(overlap_comb_df$Catinteraction_string, 
                                            "[", 2))

# to alter nutrient categories to match salt
overlap_comb_df$Stress <- "Nut_Unspecific"
overlap_comb_df_Nut <- overlap_comb_df[order(overlap_comb_df$Salt_cat),]
df_to_bind_Nut <- overlap_comb_df_Nut[,c("Number","Stress","labels")]

# to alter salt categories to match nutrient
#overlap_comb_df$Stress <- "Salt_Unspecific"
#overlap_comb_df_Salt <- overlap_comb_df[order(overlap_comb_df$Nut_cat),]
#df_to_bind_Salt <- overlap_comb_df_Salt[,c(1,6,2)]

### Nut Unspecific with Salt cats
#Unspecific_df <- rbind(NutUnSpecific_df[,c(1:3)], df_to_bind_Salt,  make.row.names=FALSE)
#Unspecific_df$Stress <- factor(Unspecific_df$Stress)
#Unspecific_df2 <- subset(Unspecific_df, Unspecific_df$Number > 0)
#Unspecific_df2$labels <- droplevels(Unspecific_df2$labels)
#Unspecific_df2$labels <- factor(Unspecific_df2$labels, 
#                                levels = c("Opposite_dir", "Increased_mag",
#                                           "Decreased_mag", "Not_DE",
#                                           "Unchanged",
#                                           "Combo_reduced.Combo_Diff_direction", 
#                                           "Combo_Salt_Same.Combo_Diff_direction",
#                                           "Combo_increased.Combo_increased",
#                                           "Combo_reduced.Combo_increased", 
#                                           "Combo_Salt_Same.Combo_increased", 
#                                           "Combo_Diff_direction.Combo_reduced",
#                                          "Combo_increased.Combo_reduced",
#                                           "Combo_reduced.Combo_reduced",
#                                           "Combo_Salt_Same.Combo_reduced",
#                                           "Salt_only.Nutrient_only",
#                                           "Combo_Diff_direction.Combo_Nutrient_Same",
#                                           "Combo_increased.Combo_Nutrient_Same",
#                                           "Combo_reduced.Combo_Nutrient_Same",
#                                           "Combo_Salt_Same.Combo_Nutrient_Same"))

### Salt Unspecific with Nutrient
df_all_full <- rbind(NutSpecific_df[,c(1:3)], df_to_bind_Nut,
                SaltUnSpecific_df[,c(1:3)], SaltSpecific_df[,c(1:3)],
                make.row.names=FALSE)
df_all_full$Stress <- factor(df_all_full$Stress,
                        levels = c("Nutrient_specific", "Nut_Unspecific",
                                   "Salt_Unspecific", "Salt_specific"))

df_all <- subset(df_all_full, df_all_full$Number > 0)
df_all$labels <- factor(df_all$labels)
df_all$labels <- droplevels(df_all$labels)
levels(df_all$labels)

# I changed the levels of labels so that the "shared" responses between salt and nutrient matched up
df_all$labels <- factor(df_all$labels, 
                        levels = c("Opposite_dir", "Increased_mag", "Decreased_mag", "Not_DE", "Unchanged",
                                    "Combo_Diff_direction.Combo_reduced", "Combo_Diff_direction.Combo_Nutrient_Same",
                                    "Combo_increased.Combo_increased", "Combo_increased.Combo_reduced",
                                    "Combo_increased.Combo_Nutrient_Same",
                                    "Combo_reduced.Combo_Diff_direction", "Combo_reduced.Combo_increased",
                                    "Combo_reduced.Combo_reduced", "Combo_reduced.Combo_Nutrient_Same",
                                    "Salt_only.Nutrient_only",
                                    "Combo_Salt_Same.Combo_Diff_direction", "Combo_Salt_Same.Combo_increased",
                                    "Combo_Salt_Same.Combo_reduced", "Combo_Salt_Same.Combo_Nutrient_Same"))

#########################
######### GRAPH #########
#########################

Opposite_col <- c("olivedrab4")
Unchanged_col <- c("#4A6990FF") #periwinkle blue
Increased_col <- c("#A73030FF")
#Increased_col <- c("#CD534CFF") #lighter red
Decreased_col <- c("#EFC000FF")
#NotDE_col <- c("#868686FF")
NotDE_col <- c("grey55")

str(df_all)

all_p <- ggplot(df_all, aes(fill=labels, y=Number, x=Stress))
all_p + geom_bar(position="stack", stat="identity", color="black", size=0.2) +
  ylab("Number of DE Genes") +
  scale_fill_manual(values=c(Opposite_col, Increased_col, Decreased_col, NotDE_col, Unchanged_col,
                             Decreased_col, Unchanged_col,
                             Increased_col, Decreased_col, Unchanged_col,
                             Opposite_col, Increased_col, Decreased_col, Unchanged_col,
                             NotDE_col,
                             Opposite_col, Increased_col, Decreased_col, Unchanged_col),
                    name="Addition of Second Stress",
                    breaks=c("Unchanged", "Decreased_mag",
                             "Increased_mag", "Opposite_dir", "Not_DE"),
                    labels=c("Un-changed",
                             "Decreased Magnitude", "Increased Magnitude",
                             "Opposite Direction", "No Longer Differentially Expressed")) +
  theme_bw(base_size = 14)
