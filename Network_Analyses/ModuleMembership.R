###### modules with a lot of genes that are correlated (positive or negative?) 
### could be further "upstream"?

#### Using module membership-
### module membership is the correlation with each gene and the module eigengene

### module eigengenes also correlate, but explain diff. proportions of variance of the genes in the module

#########################
######## SETUP ##########
#########################



#########################
###### DATAFRAME ########
#########################

GenekME <- read.csv("/Users/emilydittmar/Google Drive/Active Projects/Transcriptomics_Exp/Analyses/WGCNA/ConsensusNetwork_June2021/Module_kMEs.csv", header=T)
GenekME <- read.csv("/Users/eld72413/Google Drive/Active Projects/Transcriptomics_Exp/Analyses/WGCNA/ConsensusNetwork_June2021/Module_kMEs.csv", header=T)
length(GenekME$ID) # 34408

ModMemPvals <- read.csv("/Users/emilydittmar/Google Drive/Active Projects/Transcriptomics_Exp/Analyses/WGCNA/ConsensusNetwork_June2021/Module_kME_Pvals.csv", header=T)

# get module assignments:
Mod_Assignment <- read.csv("/Users/emilydittmar/Google Drive/Active Projects/Transcriptomics_Exp/Analyses/WGCNA/ConsensusNetwork_June2021/Gene_ModuleMembership.csv", header=T)
head(Mod_Assignment[,c(1:2)])
length(Mod_Assignment$Probe) #34408

ModMembership <- merge(Mod_Assignment[,c(1:2)], GenekME, by.x = "Probe", by.y="ID")
length(ModMembership$Probe)

length(ModMembership[which(ModMembership$ModuleLabel == "grey"), "Probe"]) # 17160


#########################
# ASSIGN N.S. CORR ZERO #
#########################

Mods <- levels(ModMembership$ModuleLabel)
Mods2 <- Mods[Mods!="grey"] # grey module not included in df

# change non-significant correlations to zeros

for (Mod in Mods2) {
  PvalName <- paste0("meta.p.equalWeights.kME", Mod)
  kME_Name <- paste0("consensus.kME", Mod)
  ModMembership[which(ModMemPvals[,PvalName] > 0.05),
                kME_Name] <- 0
}

# check
ModMembership[which(ModMemPvals[,"meta.p.equalWeights.kMEyellowgreen"] > 0.05),
              "consensus.kMEyellowgreen"]
head(ModMembership[,c(1:4)])
head(ModMemPvals[,c(1:4)])


#########################
##### STATS FUNCTION ####
#########################

#### find the mean/median correlation for each gene across all modules
#### find the number of significant correlations for each gene

AveModMember <- function(Module, df, ColNamesVector) {
  Mod_df <- df[which(df$ModuleLabel==Module),]
  ModMeans <- apply(Mod_df[,c(ColNamesVector)], MARGIN=2,
                    FUN=mean)
  AssignedModCol <- paste0("consensus.kME", Module)
  OtherMods <- ModMeans[which(names(ModMeans) != AssignedModCol)]
  Mod_median <- median(OtherMods)
  Mod_mean <- mean(OtherMods)
  NumCor <- length(ModMeans[which(OtherMods > 0)])
  AssignedMod_mean <- ModMeans[[which(names(ModMeans) == AssignedModCol)]]
  result <- data.frame(ModuleName = c(Module),
                       Median_ModCor = c(Mod_median),
                       Mean_ModCor = c(Mod_mean),
                       Num_Cor = c(NumCor),
                       Assigned_Mod_mean = c(AssignedMod_mean))
  return (result)
}

#########################
### APPLY ACROSS DATA ###
#########################

colnames <- paste0("consensus.kME", Mods2)

ModInfo <- lapply(Mods2, function(x) {AveModMember(x, ModMembership, colnames)})
ModInfoDF <- do.call("rbind", ModInfo)

hist(ModInfoDF$Median_ModCor)
hist(log(ModInfoDF$Median_ModCor)) # log transform to make normally distributed
hist(ModInfoDF$Mean_ModCor) # more normally distributed than median
hist(ModInfoDF$Num_Cor) 
hist(ModInfoDF$Assigned_Mod_mean)

#########################
##### ADD INFO TO DF ####
#########################

# Mean relatived by within-module correlations
ModInfoDF$RelMean <- ModInfoDF$Mean_ModCor / ModInfoDF$Assigned_Mod_mean
hist(ModInfoDF$RelMean)

plot(ModInfoDF$Mean_ModCor ~ ModInfoDF$Assigned_Mod_mean) # positive relationship
plot(ModInfoDF$Mean_ModCor ~ ModInfoDF$Median_ModCor) # possibly negative relationship?

ModInfoDF$RelMed <- ModInfoDF$Median_ModCor / ModInfoDF$Assigned_Mod_mean
hist(ModInfoDF$RelMed)

# add number of genes

NumGenePerMod <- aggregate(ModMembership$Probe, by=list(ModMembership$ModuleLabel), length)

ModInfoDF <- merge(ModInfoDF, NumGenePerMod, by.x = "ModuleName", by.y = "Group.1")
colnames(ModInfoDF)[8] <- c("NumGenes")

#########################
# EXPLORE RELATIONSHIPS #
#########################

# does number of genes correlate with these stats
plot(ModInfoDF$Mean_ModCor ~ ModInfoDF$NumGenes) # doesn't appear correlated (1 outlier)
plot(ModInfoDF$Median_ModCor ~ ModInfoDF$NumGenes) # positive relationship
plot(log(ModInfoDF$Median_ModCor) ~ ModInfoDF$NumGenes, xlim=c(0,600)) 
plot(log(ModInfoDF$Mean_ModCor) ~ ModInfoDF$NumGenes)

plot(ModInfoDF$Assigned_Mod_mean ~ ModInfoDF$NumGenes)

plot(ModInfoDF$Mean_ModCor ~ ModInfoDF$Assigned_Mod_mean)
     #, ylim=c(0,0.9))
abline(a=-0.6, b=1)

plot(ModInfoDF$RelMean ~ ModInfoDF$Assigned_Mod_mean) # relativizing removes relationship

plot(ModInfoDF$Num_Cor ~ ModInfoDF$NumGenes, xlim=c(0,600))

write.csv(ModInfoDF, file = "ResultsFiles/Coexpression/ModuleMembershipStats.csv")

ModInfoDF[order(ModInfoDF$RelMean),]



######### SCRATCH TESTING

#test5 <- AveModMember("plum3", ModMembership, colnames)

# test on one module
levels(ModMembership$ModuleLabel)
plum3_test <- ModMembership[which(ModMembership$ModuleLabel=="plum3"),]
hist(plum3_test$consensus.kMEplum3)
mean(plum3_test$consensus.kMEplum3) # 0.66
min(plum3_test$consensus.kMEplum3) # 0.344
max(plum3_test$consensus.kMEplum3) # 0.978
median(plum3_test$consensus.kMEplum3) # 0.666

## how does this compare to the separate networks?
mean(Mod_Assignment[which(Mod_Assignment$ModuleLabel=="plum3"),"kME.set1.MEplum3"]) # 0.8277
mean(Mod_Assignment[which(Mod_Assignment$ModuleLabel=="plum3"),"kME.set2.MEplum3"]) # 0.69

mean(plum3_test$consensus.kMEdarkslateblue) # 0.0757 (set1: 0.41, set2: 0.03)



plum3_Means <- apply(plum3_test[,c(colnames)], MARGIN=2,
                     FUN=mean)
#plum3_Means <- lapply(plum3_test[,c(colnames)], MARGIN=2,
#                      FUN=mean)

plum3_Means[which(plum3_Means>0.6)] #plum3

median(plum3_Means[which(names(plum3_Means) != "consensus.kMEplum3")]) # 0.0056 (was 0.322 before converting ns values to 0)
mean(plum3_Means[which(names(plum3_Means) != "consensus.kMEplum3")]) # 0.042
length(plum3_Means[which(names(plum3_Means) != "consensus.kMEplum3" &
                           plum3_Means > 0)]) #57

hist(plum3_Means[which(names(plum3_Means) != "consensus.kMEplum3")])
max(plum3_Means[which(names(plum3_Means) != "kME.set1.MEplum3")]) #0.66 (was 0.826 for set 1 without removing zeros)

OtherMods_plum3 <- plum3_Means[which(names(plum3_Means) != "consensus.kMEplum3")]
median(OtherMods_plum3)

num <- plum3_Means[[which(names(plum3_Means) == "consensus.kMEplum3")]]

