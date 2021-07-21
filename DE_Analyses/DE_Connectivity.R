# Connectivity and genes that are in different categories

#########################
######## SETUP ##########
#########################

library(dplyr)
library(ggplot2)
source("Functions.R")
library(car)
library(emmeans)

load("ResultsFiles/Overlap.RData")
# DEData: All data for all pairwise contrasts
# SigOverlap: Overlap between DE genes for 3 treatments
# SigDiffOverlap: Overlap between sig. pairwise differences: Combo-Nut, Combo-Salt, Nut-Salt

lapply(SigOverlap, function(x) {length(x$Gene)})
lapply(SigDiffOverlap, function(x) {length(x$Gene)})


#########################
###### DATAFRAME ########
#########################

# combine dataframes - base mean, log2FoldChange + Pvalues
my_columnsPL <- lapply(DEData, function(x) x[c(2,3,4,7)])
All_PLvalues <- do.call("cbind", my_columnsPL)
All_PLvalues <- cbind(DEData[[1]]["Gene"], All_PLvalues)

All_PLvalues[,c(1:10)]
length(All_PLvalues$Gene) #55,728

# these are the categories
load("ResultsFiles/GeneSets/MultiStressCompare.RData")
# NutSpecific_Cats, SaltSpecific_Cats, NutUnSpecific_Cats, SaltUnSpecific_Cats


load("ResultsFiles/Overlap.RData")
# DEData, SigOverlap, SigDiffOverlap

my_dataSig <- lapply(DEData, SigDEdf, PvaluesCol=7, CritP=0.05)

#########################
### ASSIGN CATEGORIES ###
#########################

# nutrient but not salt
Nutrient <- setdiff(my_dataSig$condition_nutDE_results$Gene, my_dataSig$condition_saltDE_results$Gene)
length(Nutrient) #12633
Nut_Unconditional <- setdiff(SigOverlap$`my_dataSig$condition_comboDE_results[1]my_dataSig$condition_nutDE_results[1]Only`$Gene,
                             my_dataSig$Nut_vs_Combo_results$Gene)
length(Nut_Unconditional) #2377
Nut_conditional <- union(intersect(SigOverlap$`my_dataSig$condition_comboDE_results[1]my_dataSig$condition_nutDE_results[1]Only`$Gene, 
                         my_dataSig$Nut_vs_Combo_results$Gene),
                          SigOverlap$`my_dataSig$condition_nutDE_results[1]Only`$Gene)
length(Nut_conditional) #10256

lapply(NutSpecific_Cats, length)

Salt_Unconditional <- setdiff(SigOverlap$`my_dataSig$condition_comboDE_results[1]my_dataSig$condition_saltDE_results[1]Only`$Gene,
                              my_dataSig$Salt_vs_Combo_results$Gene)
length(Salt_Unconditional) #1142

Salt_conditional <- union(intersect(SigOverlap$`my_dataSig$condition_comboDE_results[1]my_dataSig$condition_saltDE_results[1]Only`$Gene,
                              my_dataSig$Salt_vs_Combo_results$Gene),
                          SigOverlap$`my_dataSig$condition_saltDE_results[1]Only`$Gene)

length(Salt_conditional) #2088

lapply(SaltSpecific_Cats, length)

SigBoth <- SigOverlap$`my_dataSig$condition_nutDE_results[1]my_dataSig$condition_saltDE_results[1]Only`$Gene
length(SigBoth) #1855
SigAll <- SigOverlap$InCommonAll$Gene
length(SigAll) # 5032

ComboDiffs <- union(my_dataSig$Nut_vs_Combo_results$Gene,
                    my_dataSig$Salt_vs_Combo_results$Gene)

Shared_conditional <- union(intersect(SigAll, ComboDiffs),
                            SigOverlap$`my_dataSig$condition_nutDE_results[1]my_dataSig$condition_saltDE_results[1]Only`$Gene)
length(Shared_conditional) #3749

Shared_Unconditional <- setdiff(SigAll, ComboDiffs)
length(Shared_Unconditional) #3138

#########################
##### CONNECTIVITY ######
#########################

# connectivity using TOM

Connectivity <- read.csv("/Users/emilydittmar/Google Drive/Active Projects/Transcriptomics_Exp/Analyses/WGCNA/ConsensusNetwork_June2021/GeneTO.csv", header=T)
Connectivity <- read.csv("/Users/eld72413/Google Drive/Active Projects/Transcriptomics_Exp/Analyses/WGCNA/ConsensusNetwork_June2021/GeneTO.csv", header=T)
length(Connectivity$genes) # 34408
#length(Connectivity[which(is.na(Connectivity$HA_kWithin)),"probes"]) # 17160

Connectivity$genes <- gsub("\\.", ":", Connectivity$genes)

Connectivity$Category <- ifelse(Connectivity$genes %in% Nut_conditional, "Nut_conditional",
                                ifelse(Connectivity$genes %in% Nut_Unconditional, "Nut_Unconditional",
                                       ifelse(Connectivity$genes %in% Salt_Unconditional, "Salt_Unconditional",
                                              ifelse(Connectivity$genes %in% Salt_conditional, "Salt_conditional",
                                                     ifelse(Connectivity$genes %in% Shared_conditional, "Shared_conditional",
                                                            ifelse(Connectivity$genes %in% Shared_Unconditional, "Shared_Unconditional",
                                                                   "NS"))))))
#check
aggregate(Connectivity$genes, by=list(Connectivity$Category), length)
# doesn't contain values for all the genes

# add average expression
All_PLvalues$AveExp <- (All_PLvalues$condition_nutDE_results.baseMean +
  All_PLvalues$condition_saltDE_results.baseMean +
  All_PLvalues$condition_comboDE_results.baseMean) / 3


#########################
######### PLOT ##########
#########################

boxplot(Connectivity$Median ~ Connectivity$Category)
boxplot(log(Connectivity$Median) ~ Connectivity$Category)
boxplot(Connectivity$TO_sums ~ Connectivity$Category)
boxplot(log(Connectivity$TO_sums) ~ Connectivity$Category)

CombinedData <- merge(All_PLvalues[,c(1,26)], Connectivity, 
                      by.x = "Gene", by.y = "genes")
length(CombinedData$Gene)
boxplot(log(CombinedData$AveExp) ~ CombinedData$Category, notch = TRUE)
boxplot(CombinedData$TO_sums ~ CombinedData$Category, notch = TRUE)
boxplot(log(CombinedData$TO_sums) ~ CombinedData$Category, notch = TRUE)
boxplot(log(CombinedData$Median) ~ CombinedData$Category, notch = TRUE)

hist(CombinedData$TO_sums)
hist(CombinedData$Median)
hist(log(CombinedData$Median)) # normally distributed
hist(log(CombinedData$TO_sums)) 
hist(CombinedData$AveExp) # very non-normal
hist(log(CombinedData$AveExp))

CombinedData$Category <- factor(CombinedData$Category)
# kruskall-wallis test
kruskal.test(CombinedData$TO_sums, 
             CombinedData$Category)
kruskal.test(TO_sums ~ Category,
             data = CombinedData[which(CombinedData$Category != "NS"),]) # p<0.0001
kruskal.test(log(AveExp) ~ Category,
             data = CombinedData[which(CombinedData$Category != "NS"),]) # p<0.0001

## relativize by average expression?
CombinedData$TO.exp <- CombinedData$TO_sums / CombinedData$AveExp
boxplot(log(CombinedData$TO.exp) ~ CombinedData$Category, notch = TRUE)

#########################
# NUT DIFFS WITH COMBO  #
#########################

# relationship with connectivity and differences in expression from single stress-combo?

CombinedData <- merge(All_PLvalues[,c(1,15,23,26)], Connectivity, 
                      by.x = "Gene", by.y = "genes")

p <- ggplot(data=CombinedData[which(CombinedData$Gene %in% my_dataSig$condition_nutDE_results$Gene),], 
            aes(x=log(TO_sums), y=abs(Nut_vs_Combo_results.log2FoldChange)))
p_graph <- p + 
  #geom_point(aes(color=Category), alpha=0.9) +
  geom_smooth(method=lm, aes(color=Category)) +
  theme_minimal() +
  #ylim(0,2) +
  ylab("|Log fold change between Nutrient and Combo|") +
  xlab("Log Topological Overlap Sum per gene")

p2 <- ggplot(data=CombinedData[which(CombinedData$Gene %in% my_dataSig$condition_nutDE_results$Gene),], 
            aes(x=log(AveExp), y=abs(Nut_vs_Combo_results.log2FoldChange)))
p2 + 
  #geom_point(aes(color=Category))
  geom_smooth(method=lm, aes(color=Category)) +
  theme_minimal()

p3 <- ggplot(data=CombinedData[which(CombinedData$Gene %in% my_dataSig$condition_nutDE_results$Gene),], 
             aes(x=log(TO_sums), y=log(AveExp)))
p3 + 
  #geom_point(aes(color=Category))
  geom_smooth(method=lm, aes(color=Category)) +
  theme_minimal()

p4 <- ggplot(data=CombinedData[which(CombinedData$Gene %in% my_dataSig$condition_nutDE_results$Gene),], 
            aes(x=log(Median), y=abs(Nut_vs_Combo_results.log2FoldChange)))
p4 + 
  #geom_point(aes(color=Category))
  geom_smooth(method=lm, aes(color=Category)) +
  theme_minimal()

### add module info (taken from other script)
p_module <- ggplot(data=LM_wLog2FoldChange, aes(x=log(TO_Sum_Mean), 
                   y=abs(Log2Fold_Nut.Combo)
                   ))
p_module+geom_point()
                   
p_graph + geom_point(data=LM_wLog2FoldChange, 
                     mapping=aes(x=log(TO_Sum_Mean), 
                         y=abs(Log2Fold_Nut.Combo)))

#########################
## SIGNIFICANT SLOPES? ##
#########################
hist(abs(CombinedData$Nut_vs_Combo_results.log2FoldChange))

Nut_mod <- lm(abs(Nut_vs_Combo_results.log2FoldChange) ~ log(TO_sums) +
                    log(AveExp) + Category +
                Category:log(TO_sums), 
              data=CombinedData[which(CombinedData$Gene %in% my_dataSig$condition_nutDE_results$Gene),])
hist(resid(Nut_mod))
summary(Nut_mod)
drop1(Nut_mod, test = "Chi")
Anova(Nut_mod)


#########################
# SALT DIFFS WITH COMBO  #
#########################

q <- ggplot(data=CombinedData[which(CombinedData$Gene %in% my_dataSig$condition_saltDE_results$Gene),], 
            aes(x=log(TO_sums), y=abs(Salt_vs_Combo_results.log2FoldChange)))
q + 
  #geom_point(aes(color=Category)) +
  geom_smooth(method=lm, aes(color=Category)) +
  theme_minimal() +
  ylab("Log fold change between Salt and Combo") +
  xlab("Log Topological Overlap Sum per gene")

q2 <- ggplot(data=CombinedData[which(CombinedData$Gene %in% my_dataSig$condition_saltDE_results$Gene),], 
            aes(x=log(AveExp), y=abs(Salt_vs_Combo_results.log2FoldChange)))
q2 + 
  #geom_point(aes(color=Category)) +
  geom_smooth(method=lm, aes(color=Category)) +
  theme_minimal()

q3 <- ggplot(data=CombinedData[which(CombinedData$Gene %in% my_dataSig$condition_saltDE_results$Gene),], 
             aes(x=log(TO_sums), y=log(AveExp)))
q3 + 
  #geom_point(aes(color=Category)) +
  geom_smooth(method=lm, aes(color=Category)) +
  theme_minimal()

#########################
# MODULE ASSIGNMENTS  #
#########################

Mod_Assignment <- read.csv("/Users/emilydittmar/Google Drive/Active Projects/Transcriptomics_Exp/Analyses/WGCNA/ConsensusNetwork_June2021/Gene_ModuleMembership.csv", header=T)
head(Mod_Assignment[,c(1:2)])
length(Mod_Assignment$Probe)

Mod_Assignment$Probe <- gsub("\\.", ":", Mod_Assignment$Probe)

Combined2 <- merge(CombinedData, Mod_Assignment[,c(1:2)], by.x = "Gene", by.y = "Probe")


Log2Fold_Nut.Combo <- ConnectivitySummary(Combined2$Nut_vs_Combo_results.log2FoldChange,
                            Combined2$ModuleLabel, "Module")

Log2Fold_Salt.Combo <- ConnectivitySummary(Combined2$Salt_vs_Combo_results.log2FoldChange,
                                          Combined2$ModuleLabel, "Module")


Log2FoldChange_combo <- cbind.data.frame(Log2Fold_Nut.Combo[,c(1:3,5:6)],
                                         Log2Fold_Salt.Combo[,c(3,5)])
colnames(Log2FoldChange_combo)[c(3,4,6,7)] <- c("Log2Fold_Nut.Combo", "SE_Log2Fold_Nut.Combo",
                                                "Log2Fold_Salt.Combo", "SE_Log2Fold_Salt.Combo")


############# SCRATCH below
Unconditional <- union(
              setdiff(SigOverlap$`my_dataSig$condition_comboDE_results[1]my_dataSig$condition_nutDE_results[1]Only`$Gene, 
                      my_dataSig$Nut_vs_Combo_results$Gene), #2377
                  union(
              setdiff(SigOverlap$`my_dataSig$condition_comboDE_results[1]my_dataSig$condition_saltDE_results[1]Only`$Gene, 
                      my_dataSig$Salt_vs_Combo_results$Gene), #1142
              setdiff(Shared, my_dataSig$Nut_vs_Salt_results$Gene) #4048
                ))

Conditional <- union(
            setdiff(Nutrient, Unconditional), #10256
                union(
            setdiff(Salt, Unconditional),
            setdiff(Shared, Unconditional)
            ))

All_PLvalues$Stress <- as.factor(ifelse(All_PLvalues$Gene %in% Nutrient,
                              "Nutrient", ifelse(All_PLvalues$Gene %in% Salt,
                                                 "Salt", ifelse(All_PLvalues$Gene %in% Shared,
                                                                "Shared", "NS"))))
All_PLvalues$Dependence <- as.factor(ifelse(All_PLvalues$Gene %in% Conditional,
                                  "Conditional", ifelse(All_PLvalues$Gene %in% Unconditional,
                                                        "Unconditional", "NS")))

# check
aggregate(All_PLvalues$Gene, 
          by=list(All_PLvalues$Stress, All_PLvalues$Dependence), length)
# check
lapply(NutSpecific_Cats, function(x) {length(x)})
#nutrient conditional: 8636 + 1576 +6 +38= 10,256
# nutrient unconditional: 2377 

lapply(SaltSpecific_Cats, function(x) {length(x)})
# salt conditional: 2018 + 62 + 1 + 7 = 2088
# salt unconditional: 1142

#########################
##### CONNECTIVITY ######
#########################

Connectivity <- read.csv("/Users/emilydittmar/Google Drive/Active Projects/Transcriptomics_Exp/Analyses/WGCNA/ConsensusNetwork_June2021/Gene_Connectivity.csv", header=T)
Connectivity <- read.csv("/Users/eld72413/Google Drive/Active Projects/Transcriptomics_Exp/Analyses/WGCNA/ConsensusNetwork_June2021/Gene_Connectivity.csv", header=T)
length(Connectivity$probes) # 34408
length(Connectivity[which(is.na(Connectivity$HA_kWithin)),"probes"]) # 17160

Connectivity$probes <- gsub("\\.", ":", Connectivity$probes)

CombinedData <- merge(All_PLvalues, Connectivity, by.x = "Gene", by.y = "probes")
length(CombinedData$Gene) #34,408

CombinedData$kTotal_perGene <- CombinedData$kTotal / 34408
hist(CombinedData$kTotal_perGene)
min(CombinedData$kTotal_perGene) # 0.0156
max(CombinedData$kTotal_perGene) # 0.06
# does this range make sense?

#############################
# CONNECTIVITY RELATIONSHIPS #
#############################

Mod_Connect <- function(Yvar, dataset) {
  mod <- lm(Yvar ~ Stress +
              Dependence +
              Stress:Dependence, data=dataset)
  return(mod)
}

#############################
# CONNECTIVITY PER CATEGORY #
#############################

# relationship between expression level and connectivity?
plot(log(CombinedData$condition_saltDE_results.baseMean) ~
       CombinedData$kTotal) # no obvious relationship

boxplot(CombinedData$kTotal ~ CombinedData$Dependence)
boxplot(CombinedData$kTotal ~ CombinedData$Stress)

# significant relationship with kTotal
CombinedDataSig <- subset(CombinedData, Stress!="NS")
Cmod <- Mod_Connect(CombinedDataSig$kTotal, CombinedDataSig)
summary(Cmod)
Anova(Cmod, type = "II", test = "F") # all factors significant
emmeans(Cmod, ~Dependence|Stress, type="response")
pairs(emmeans(Cmod, ~Dependence|Stress, type="response")) # all pairwise differences are significant (nutrient p-value=0.0001, salt and shared are <0.0001)
pairs(emmeans(Cmod, ~Dependence*Stress, type="response"))
# only pairs that are not different: 
# unconditional nutrient and conditional salt, conditional salt and unconditional shared (p=0.0295)

CombinedData$Category <- paste0(CombinedData$Stress, "_", CombinedData$Dependence)

par(mfrow=c(2,1))
boxplot(CombinedData$kTotal ~ CombinedData$Category)
boxplot(log(CombinedData$condition_comboDE_results.baseMean) ~ CombinedData$Category)

# look at intra-modular versus inter-modular connectivity?
par(mfrow=c(2,1))
boxplot(CombinedData$kWithin ~ CombinedData$Category)
boxplot(CombinedData$kOut ~ CombinedData$Category)

# proportion within modular connectivity?
Cmod2 <- Mod_Connect(CombinedDataSig$kWithin/CombinedDataSig$kTotal, CombinedDataSig)
Anova(Cmod2, type = "II", test = "F")
emmeans(Cmod2, ~Dependence|Stress, type="response") # unconditional lower, small difference
pairs(emmeans(Cmod2, ~Dependence|Stress, type="response")) #nutrient (p=0.0439), shared (p=0.0572), salt (p<0.0001)
pairs(emmeans(Cmod2, ~Dependence*Stress, type="response"))

# inter modular connectivity?
Cmod3 <- Mod_Connect(CombinedDataSig$kOut/CombinedDataSig$kTotal, CombinedDataSig)
Anova(Cmod3, type = "II", test = "F")
emmeans(Cmod3, ~Dependence|Stress, type="response") # conditional lower, small difference
pairs(emmeans(Cmod3, ~Dependence|Stress, type="response"))
pairs(emmeans(Cmod3, ~Dependence*Stress, type="response"))

Cmod4 <- Mod_Connect(CombinedDataSig$kWithin/CombinedDataSig$kOut, CombinedDataSig)
Anova(Cmod4, type = "II", test = "F")
emmeans(Cmod4, ~Dependence|Stress, type="response") # unconditional lower, small difference
pairs(emmeans(Cmod4, ~Dependence|Stress, type="response"))
pairs(emmeans(Cmod4, ~Dependence*Stress, type="response"))
pairs(emmeans(Cmod4, ~Stress, type="response"))

par(mfrow=c(2,1))
boxplot(CombinedData$kWithin/CombinedData$kTotal ~ CombinedData$Category,
        main="Intramodular connectivity")
boxplot(CombinedData$kOut/CombinedData$kTotal ~ CombinedData$Category, 
        main="Intermodular connectivity")

#############################
## EXPRESSION PER CATEGORY ##
#############################

# differences in expression?
hist(log(CombinedDataSig$condition_comboDE_results.baseMean))
boxplot(log(CombinedDataSig$condition_comboDE_results.baseMean) ~ CombinedDataSig$Category)
boxplot(log(CombinedData$condition_nutDE_results.baseMean) ~ CombinedDataSig$Category)
boxplot(log(CombinedData$condition_saltDE_results.baseMean) ~ CombinedDataSig$Category)

plot(log(CombinedDataSig$condition_comboDE_results.baseMean), CombinedDataSig$kTotal)

Expmod <- lm(log(condition_comboDE_results.baseMean) ~ Stress +
             Dependence +
             Stress:Dependence, data=CombinedDataSig)
summary(Expmod)
Anova(Expmod, type = "II", test = "F") # all factors significant

emmeans(Expmod, ~Dependence|Stress, type="response")
pairs(emmeans(Expmod, ~Dependence|Stress, type="response"))
# no difference with salt

pairs(emmeans(Expmod, ~Dependence*Stress, type="response"))
# not significant pairs: 
# conditional nutrient/unconditional salt, unconditional nutrient with conditional&unconditional salt & unconditional shared
# 

#############################
## Diff Single vs. Multi  ##
#############################

Data_stress <- split(CombinedData, CombinedData$Stress)

par(mfrow=c(1,2))
boxplot(abs(Data_stress$Nutrient$Nut_vs_Combo_results.log2FoldChange) ~
          Data_stress$Nutrient$Category)
boxplot(abs(log(Data_stress$Nutrient$condition_nutDE_results.baseMean)) ~
          Data_stress$Nutrient$Category)

plot(abs(Data_stress$Nutrient$Nut_vs_Combo_results.log2FoldChange) ~
          Data_stress$Nutrient$kTotal)


# look a the magnitude of differences in the nut vs. combo log2foldchange as a function of connectivity

p1 <- ggplot(data=Data_stress$Nutrient, 
             aes(x=kTotal, y=abs(Nut_vs_Combo_results.log2FoldChange)))
p1 + geom_point(aes(color=Category)) +
  geom_smooth()


#############################
## ACCOUNT FOR GENE NUMBER##
#############################

ModMembership <- read.csv("/Users/emilydittmar/Google Drive/Active Projects/Transcriptomics_Exp/Analyses/WGCNA/ConsensusNetwork_June2021/Gene_ModuleMembership.csv", header=T)
length(ModMembership$Probe) # 34408
length(ModMembership[which(ModMembership$ModuleLabel == "grey"), "Probe"]) # 17160
ModMembership$probes <- gsub("\\.", ":", ModMembership$Probe)

# merge dataframes
GeneInfo <- merge(ModMembership[,c(2,537)], Connectivity, by = "probes")
length(GeneInfo$probes) # 34408

# take out genes that aren't assigned to modules (i.e. grey)
GeneModInfo <- subset(GeneInfo, ModuleLabel != "grey")
length(GeneModInfo$probes) # 17248

# does number of genes in module affect within-module connectivity?
NumGenes <- aggregate(GeneModInfo$probes, by=list(GeneModInfo$ModuleLabel), length)
sum(NumGenes$x) #17248


### need to use 
######## SCRATCH ############


# not DE in salt, not diff between combo-nut
NutSpecific_Unconditional <- setdiff(SigOverlap$`my_dataSig$condition_comboDE_results[1]my_dataSig$condition_nutDE_results[1]Only`$Gene,
                                     my_dataSig$Nut_vs_Combo_results$Gene)
length(NutSpecific_Unconditional) #2377

# not DE in salt, diff between combo-nut OR not expressed in combo
NutSpecific_Conditional <- union(
  intersect(SigOverlap$`my_dataSig$condition_comboDE_results[1]my_dataSig$condition_nutDE_results[1]Only`$Gene,
            my_dataSig$Nut_vs_Combo_results$Gene),
  SigOverlap$`my_dataSig$condition_nutDE_results[1]Only`$Gene
)
length(NutSpecific_Conditional) #10256

# check
lapply(SigOverlap, function(x) {length(x$Gene)})
# combo + nut: 3997
# nut only: 8636
# total: 12,633
lapply(NutSpecific_Cats, function(x) {length(x)})

SaltSpecific_Unconditional <- setdiff(SigOverlap$`my_dataSig$condition_comboDE_results[1]my_dataSig$condition_saltDE_results[1]Only`$Gene,
                                      my_dataSig$Salt_vs_Combo_results$Gene)
length(SaltSpecific_Unconditional) #1142

SaltSpecific_Conditional <- union(
  intersect(SigOverlap$`my_dataSig$condition_comboDE_results[1]my_dataSig$condition_saltDE_results[1]Only`$Gene,
            my_dataSig$Salt_vs_Combo_results$Gene),
  SigOverlap$`my_dataSig$condition_saltDE_results[1]Only`$Gene
)
length(SaltSpecific_Conditional) # 70 + 2088
# total: 3230
# combo+salt: 1212
# salt only: 2018
lapply(SaltSpecific_Cats, function(x) {length(x)})

# shared genes:
# DE in salt & nut, nut diff between them
Shared_same <- setdiff(
  intersect(my_dataSig$condition_nutDE_results$Gene, my_dataSig$condition_saltDE_results$Gene),
  my_dataSig$Nut_vs_Salt_results$Gene
)
length(Shared_same) #4048

Shared_diff <- intersect(
  intersect(my_dataSig$condition_nutDE_results$Gene, my_dataSig$condition_saltDE_results$Gene),
  my_dataSig$Nut_vs_Salt_results$Gene
)
length(Shared_diff) #2839

# 5032+1855=6,887

##########################
# FACTORS IN SPREADSHEET #
##########################

All_PLvalues$Category <- ifelse(All_PLvalues$Gene %in% NutSpecific_Unconditional,
                                "Nutrient_unconditional", ifelse(All_PLvalues$Gene %in% NutSpecific_Conditional,
                                                                 "Nutrient_conditional", ifelse(All_PLvalues$Gene %in% SaltSpecific_Unconditional,
                                                                                                "Salt_unconditonal", ifelse(All_PLvalues$Gene %in% SaltSpecific_Conditional,
                                                                                                                            "Salt_conditional", ifelse(All_PLvalues$Gene %in% Shared_same,
                                                                                                                                                       "Shared_NotDiff", ifelse(All_PLvalues$Gene %in% Shared_diff,
                                                                                                                                                                                "Shared_Diff", "NS"))))))

aggregate(All_PLvalues$Gene, by=list(All_PLvalues$Category), length)
All_PLvalues$Category <- factor(All_PLvalues$Category)


CombinedData$StressType <- "NS"
CombinedData[grep("Nutrient", CombinedData$Category),"StressType"] <- "Nutrient"
CombinedData[grep("Salt", CombinedData$Category),"StressType"] <- "Salt"
CombinedData[grep("Shared", CombinedData$Category),"StressType"] <- "Shared"




