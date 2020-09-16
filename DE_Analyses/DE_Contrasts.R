#########################
######## SETUP ##########
#########################

library(dplyr)
library(ggplot2)
source("Functions.R")
library(UpSetR)


#library(viridis)
#library(RColorBrewer)



#########################
###### DE OVERLAP #######
#########################

# read in dataframes
DEData <- ImportCSVs("ResultsFiles/GeneSets", 0.05) 
# function that imports all csv's in directory as a list, sorts by Gene name, changes first column to "Gene"

my_dataSig <- lapply(DEData, SigDEdf, PvaluesCol=7, CritP=0.05) #list that contains only the genes that are significant

SigOverlap <- GeneSets(my_dataSig$condition_comboDE_results[1], my_dataSig$condition_nutDE_results[1], 
                       my_dataSig$condition_saltDE_results[1]) #different sets of overlap for DE genes
SigDiffOverlap <- GeneSets(my_dataSig$Nut_vs_Combo_results[1], my_dataSig$Nut_vs_Salt_results[1], 
                           my_dataSig$Salt_vs_Combo_results[1]) #different sets of overlap for pairwise treatment contrasts
names(SigDiffOverlap)

lapply(SigOverlap, function(x) {length(x$Gene)}) # numbers for Venn Diagram
#5032 in common all (21%), 3997 Combo + Nut (17%), 1855 Nut + Salt (8%), 1212 Combo + Salt (5%),
# 1039 Combo only (4%), 8636 Nut only (36%), 2018 salt only (8%)
#23,791 total

#########################
###### UPSET PLOT #######
#########################

SigOverlap_Graph <- lapply(my_dataSig, function(x) {x$Gene})

upset(fromList(SigOverlap_Graph),
      #keep.order = TRUE,
      order.by = "freq",
      #group.by = "sets", 
      nsets = 13,
      #empty.intersections = "on",
      nintersects = 20)

#########################
###### UP vs DOWN #######
#########################

my_dataSigUp <- lapply(my_dataSig, MoreCritNum, column=3, critNum=0)
my_dataSigDown <- lapply(my_dataSig, LessCritNum, column=3, critNum=0)

SigOverlapUp_Graph <- lapply(my_dataSigUp[1:3], function(x) {x$Gene})
SigOverlapDown_Graph <- lapply(my_dataSigDown[1:3], function(x) {x$Gene})

upset(fromList(c(SigOverlapUp_Graph, SigOverlapDown_Graph)),
      #keep.order = TRUE,
      order.by = "freq",
      #group.by = "sets", 
      nsets = 13,
      #empty.intersections = "on",
      nintersects = 20)


#########################
###### NUT x SALT #######
#########################

# how many genes shared by nutrient + salt (but not combo) are in the opposite direction
# 1855 genes found in nutrient + salt only

length(intersect(my_dataSigUp$condition_nutDE_results$Gene, 
                 my_dataSigDown$condition_saltDE_results$Gene)) #389 (313 ns in combo, 67 combo is down)
length(intersect(my_dataSigUp$condition_saltDE_results$Gene, 
                 my_dataSigDown$condition_nutDE_results$Gene)) #366 (293 ns in combo, 71 combo is up)

SaltNutUniq <-
  intersect(SigOverlap$`my_dataSig$condition_nutDE_results[1]my_dataSig$condition_saltDE_results[1]Only`$Gene,
            SigDiffOverlap$`my_dataSig$Nut_vs_Combo_results[1]my_dataSig$Salt_vs_Combo_results[1]Only`$Gene)
# N=61


