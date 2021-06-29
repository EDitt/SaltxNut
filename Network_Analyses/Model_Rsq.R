
library(ggplot2)

# assessing relative importance of variables-
#### modules for which model with accession alone explains more than 90% of variance:

##############################
####### MODEL FUNCTION #######
##############################

Rsq_Mods <- function(Yvar, dataset) {
  mod1 <- lm(Yvar ~ Accession,
               data = dataset,
              contrast=list(Accession=contr.sum))
  AccRsq <- summary(mod1)$adj.r.squared
  AccP <- anova(mod1)[1,5]
  mod2 <- lm(Yvar ~ Treatment,
            data = dataset)
  TreatRsq <- summary(mod2)$adj.r.squared
  TreatP <- anova(mod2)[1,5]
  return(data.frame(row.names = c("Accession", "Treatment"),
                    Adj.r.squared = c(AccRsq, TreatRsq),
                    Pvalue = c(AccP, TreatP)))
}

# test
Rsq_Mods(AllData$MEbisque4, AllData)
Rsq_Mods(AllData$MEwhite, AllData)["Accession", "Adj.r.squared"]

################################
# HIGH R-SQUARED FOR ACCESSION #
################################

# run on all modules
Module_Rsq_results <- lapply (AllData[,c(17:104)], function(x) {Rsq_Mods(x, AllData) })
str(Module_Rsq_results)

# which modules have Accession explaining more than 90% of variation?

HighR_accession <- lapply (Module_Rsq_results, function(x) {which(x["Accession", "Adj.r.squared"] > 0.9) })
length(which(HighR_accession==1)) #34
HighR_accessionMods <- which(HighR_accession==1)

critp <- 0.05/length(AllData[,c(17:104)])

# are these significant?
HighR_accessionP <- lapply (Module_Rsq_results[names(HighR_accessionMods)], function(x) 
  {which(x["Accession", "Pvalue"] < critp) })
which(HighR_accessionP==0) # none

# for these modules, is the treatment model significant? what is the R^2 for treatment?
HighR_accessionTreatP <- lapply (Module_Rsq_results[names(HighR_accessionMods)], function(x) 
  {which(x["Treatment", "Pvalue"] < critp) })
which(HighR_accessionTreatP==1) # none

# what is the R^2 of these?
HighR_accessionTreatRsq <- lapply (Module_Rsq_results[names(HighR_accessionMods)], function(x) 
{x["Treatment", "Adj.r.squared"]})


HighR_accessionTreatRsq_vals <- do.call("rbind", HighR_accessionTreatRsq)
hist(HighR_accessionTreatRsq_vals)

### none of these modules have significant or positive R^2 values

####################################
# FUNCTION TO PLOT MODULE RAW DATA #
####################################

base_barplot <- function(var1, module_name){
  barplot(var1, main = module_name)
}

# put data in order of treatment and accession for visualization
AllData_ordered <- AllData[order(AllData$Treatment, AllData$Accession, decreasing = FALSE),]

################################
# GRAPH HIGH R^2 ACCESSION MODS #
################################

# use modules of interest to make a list for var1
highAccMods_yvals <- lapply(names(HighR_accessionMods), function(x) {AllData_ordered[,x]})

length(HighR_accessionMods) #34

# plot by 16 to a page
par(mfrow=c(4,4))
par(mar=c(1,2,1,1))
mapply(baseplot, c(highAccMods_yvals[1:16]), c(names(HighR_accessionMods[1:16]))) 

par(mfrow=c(5,4))
par(mar=c(1,2,1,1))
mapply(baseplot, c(highAccMods_yvals[17:34]), c(names(HighR_accessionMods[17:34]))) 


################################
###### COMBINE ALL DATA ########
################################

# add column for module identity
Module_Rsq_resultswLabels <-
  lapply(names(Module_Rsq_results), function(x) 
    { Module_Rsq_results[[x]]$Module <- x;return(Module_Rsq_results[[x]])})

All_accession <- lapply(Module_Rsq_resultswLabels, function(x) {x[1,]})
All_accessionDf <- do.call("rbind", All_accession)
colnames(All_accessionDf) <- c("Accession_Rsquared", "Accession_Pvalue", "Module")

All_treat <- lapply(Module_Rsq_resultswLabels, function(x) {x[2,]})
All_treatDf <- do.call("rbind", All_treat)
colnames(All_treatDf) <- c("Treatment_Rsquared", "Treatment_Pvalue", "Module")

All_RsqData <- merge(All_accessionDf, All_treatDf, by="Module", rownames=NULL)

All_RsqData$AccessionSig <- ifelse(All_RsqData$Accession_Pvalue < 0.05, "Significant", "NS")
All_RsqData$TreatmentSig <- ifelse(All_RsqData$Treatment_Pvalue < 0.05, "Significant", "NS")

aggregate(All_RsqData$Module, by=list(All_RsqData$AccessionSig, All_RsqData$TreatmentSig), length)
# 40 significant for accession only, 22 significant for treatment only, 24 significant for both, 2 n.s. for both

################################
###### COEFFICIENT GRAPHS ######
################################

dev.off()

#### what is the relationship between R^2 for accession and R^2 for treatment?
plot(All_RsqData$Accession_Rsquared ~ All_RsqData$Treatment_Rsquared,
     xlab="Treatment R^2", ylab = "Accession R^2")
abline(h=0.5, lty=2)
abline(h=0.9, lty=2)

# looks like a few very high R^2 values for accession - only ones over 0.9 are above 0.5

hist(All_RsqData$Accession_Rsquared, xlab="Accession R^2", main="")

### scatterplot with significance color-coded?
r_plot <- ggplot(data=All_RsqData, aes(x=Treatment_Rsquared, y=Accession_Rsquared))
r_plot + geom_point(aes(color=TreatmentSig,
                        shape=AccessionSig), size=2) +
  theme_minimal() +
  scale_shape_manual(values=c(1,19)) +
  geom_hline(yintercept = c(0.9, 0.5), linetype="dotted")
ggsave("ResultsFiles/Coexpression/QC/AccessionTreat_R2_color.pdf")

# what is the lowest value over 0.9?
min(All_RsqData[which(All_RsqData$Accession_Rsquared > 0.9),"Accession_Rsquared"]) # 0.963 !

#####################################
# MODERATE R-SQUARED FOR ACCESSION? #
#####################################

# a few R^2 values for accession are close to 0.5. I want to graph the raw values

# sort by highest R^2 values
All_RsqData_sorted <- All_RsqData[order(All_RsqData$Accession_Rsquared, decreasing = TRUE),]

Mod_Acc_Mods <- All_RsqData_sorted[which(All_RsqData_sorted$Accession_Rsquared > 0.25 &
                                           All_RsqData_sorted$Accession_Rsquared < 0.9 &
                                           All_RsqData_sorted$TreatmentSig=="NS"),"Module"] 
# R^2 ranges from 0.33-0.49 (N=6)

ModAccMods_yvals <- lapply(Mod_Acc_Mods, function(x) {AllData_ordered[,x]})

# put R^2 value with module name on graph
R2_vals <- All_RsqData_sorted[which(All_RsqData_sorted$Module %in% Mod_Acc_Mods),"Accession_Rsquared"]
round(R2_vals, digits = 2)
Names <- paste0(Mod_Acc_Mods, " Accession R^2=", round(R2_vals, digits = 2))

par(mfrow=c(3,2))
par(mar=c(1,2,1,1))
mapply(baseplot, c(ModAccMods_yvals), c(Names)) 

####################################
# RAW GRAPHS FOR SIG R^2 TREATMENT #
####################################

# same graphs (R^2 for accession between 0.25-0.5, but also significant for treatment)

Mod_Acc_Mods_Treatsig <- All_RsqData_sorted[which(All_RsqData_sorted$Accession_Rsquared > 0.25 &
                                           All_RsqData_sorted$Accession_Rsquared < 0.9 &
                                           All_RsqData_sorted$TreatmentSig=="Significant"),"Module"] 

ModAccMods_yvals_Treatsig <- lapply(Mod_Acc_Mods_Treatsig, function(x) {AllData_ordered[,x]})

# put R^2 value with module name on graph
R2_vals <- All_RsqData_sorted[which(All_RsqData_sorted$Module %in% Mod_Acc_Mods_Treatsig),
                              "Accession_Rsquared"]
# accession R^2 = 0.32-0.46, N=20

Names <- paste0(Mod_Acc_Mods_Treatsig, " Acc. R^2=", round(R2_vals, digits = 2))

par(mfrow=c(5,4))
par(mar=c(1,2,1,1))
mapply(baseplot, c(ModAccMods_yvals_Treatsig), c(Names)) 

###################################
###### REMAINING RAW GRAPHS #######
###################################

# graph the remaining modules with significant R2 for treatment that are also significant for accession (N=4)
LowAccR2_TreatSig <- All_RsqData_sorted[which(All_RsqData_sorted$Accession_Rsquared < 0.25 &
                                                    All_RsqData_sorted$TreatmentSig=="Significant" &
                                                    All_RsqData_sorted$AccessionSig=="Significant"),"Module"] 

LowAccR2_TreatSig_yvals <- lapply(LowAccR2_TreatSig, function(x) {AllData_ordered[,x]})
R2_vals <- All_RsqData_sorted[which(All_RsqData_sorted$Module %in% LowAccR2_TreatSig),
                              "Accession_Rsquared"] # accession R^2 = 0.06-0.22, N=4
Names <- paste0(LowAccR2_TreatSig, " Acc. R^2=", round(R2_vals, digits = 2))
par(mfrow=c(2,2))
par(mar=c(1,2,1,1))
mapply(baseplot, c(LowAccR2_TreatSig_yvals), c(Names)) 

# modules that are significant for treatment only (N=22)
TreatSigAccNS <- All_RsqData_sorted[which(All_RsqData_sorted$TreatmentSig=="Significant" &
                                          All_RsqData_sorted$AccessionSig=="NS"),"Module"]
TreatSigAccNS_yvals <- lapply(TreatSigAccNS, function(x) {AllData_ordered[,x]})
Names <- paste0(TreatSigAccNS, " Acc. R^2 n.s")
par(mfrow=c(4,3))
par(mar=c(1,2,1,1))
mapply(baseplot, c(TreatSigAccNS_yvals[1:11]), c(Names[1:11]))
mapply(baseplot, c(TreatSigAccNS_yvals[12:22]), c(Names[12:22]))
 
# modules that are n.s. in both (N=4)
TreatNSAccNS <- All_RsqData_sorted[which(All_RsqData_sorted$TreatmentSig=="NS" &
                                            All_RsqData_sorted$AccessionSig=="NS"),"Module"]

TreatNSAccNS_yvals <- lapply(TreatNSAccNS, function(x) {AllData_ordered[,x]})
Names <- paste0(TreatNSAccNS, " Acc. R^2 n.s")
par(mfrow=c(2,1))
par(mar=c(1,2,1,1))
mapply(baseplot, c(TreatNSAccNS_yvals), c(Names))

###########################
#### SAVE MODULE LISTS ####
###########################

aggregate(All_RsqData$Module, by=list(All_RsqData$AccessionSig, All_RsqData$TreatmentSig), length)

# >90% R^2 for accession (all NS for treatment)
HighR_accession <- names(HighR_accessionMods) # N=34

# 0.25-0.5% R^2 for accession, N.S. R^2 for treatment (N=6)
Mod_Acc_Mods

# 0.25-0.5% R^2 for accession, sig. R^2 for treatment (N=20)
Mod_Acc_Mods_Treatsig

# <0.25% R^2 for accession, sig R^2 for treatment (N=4)
LowAccR2_TreatSig

# R^2 significant for treatment only (N=22)
TreatSigAccNS

# N.S. R^2 for both (N=2)
TreatNSAccNS


ModuleRsquaredList <- list(HighR_accession, Mod_Acc_Mods, Mod_Acc_Mods_Treatsig, LowAccR2_TreatSig, TreatSigAccNS, TreatNSAccNS)
names(ModuleRsquaredList) <- c("AccessionR2_more95", "AccessionR2_25-50_nsTreat", "AccessionR2_25-50_sigTreat",
                               "AccessionR2_less25_sigTreat", "AccessionR2_NS_sigTreat", "AccessionR2_NS_nsTreat")

save(ModuleRsquaredList, file="ResultsFiles/Coexpression/QC/Module_RsquaredCats.RData")


####################################
# OVERLAP W FULL MODEL SIGNIFICANE #
####################################
# how many of these overlap with modules identified as significant for a treatment factor in the full model?


lapply(SigModOverlap, function(x) {intersect(x, HighR_accession)})
# 12 overlap the Nutrient only category

HighAccR2_SigLRNut <- intersect(SigModOverlap$LR_Nut_SigOnly, HighR_accession)
# darkmagenta, yellowgreen, brown2, lightgreen, midnightblue, royalblue, violet,
# mediumorchid, lightyellow, paleturquoise, coral1, plum3

lapply(SigModOverlap, function(x) {intersect(x, Mod_Acc_Mods)})
# no overlaps

lapply(SigModOverlap, function(x) {intersect(x, Mod_Acc_Mods_Treatsig)})
# 4 in common all, 4 Nut+Salt, 9 Nut+NutxSalt, 3 Nut only

lapply(SigModOverlap, function(x) {intersect(x, LowAccR2_TreatSig)})
# 1 Nut+NutxSalt, 1 Nut only, 1 Salt only

lapply(SigModOverlap, function(x) {intersect(x, TreatSigAccNS)})
# 2 in common all, 3 Nut+Salt, 9 Nut+NutxSalt, 4 Nut only, 1 Salt only, 1 NutxSalt only

lapply(SigModOverlap, function(x) {intersect(x, TreatNSAccNS)})
# no overlap



####################################
####################################
####################################
### scratch below

# assessing relative importance of variables?

LR_mod_results$MEbisque4

summary(LR_mod_results$MEbisque4)

test <- lm(MEbisque4 ~
            Accession, data = AllData,
          contrast=list(Accession=contr.sum))

summary(test) # R^2=0.9821***

test_full <- lm(MEbisque4 ~ Treatment +
                  Accession, data = AllData,
                contrast=list(Accession=contr.sum))
summary(test_full) #R^2=0.9819***

anova(test_full, test) # ns
drop1(test_full, test = "Chi") # treatment is ns

# modules that have treatment & accession as significant
intersect(LR_Accession_Sig, LR_Nut_Sig)

# MEdarkmagenta- significant for Nut main effect & accession, but pattern looks driven by accession
# MEwhite - sig. for Nut and accession but more treatment effect

Mod_darkmagenta1 <- lm(MEdarkmagenta ~
             Accession, data = AllData,
           contrast=list(Accession=contr.sum))
summary(Mod_darkmagenta1) #R^2=0.9859

Mod_darkmagenta2 <- lm(MEdarkmagenta ~
                         Accession + Treatment, data = AllData,
                       contrast=list(Accession=contr.sum))
summary(Mod_darkmagenta2) #R^2=0.9878

anova(Mod_darkmagenta2, Mod_darkmagenta1) #p=0.0009256
# not less than critp
Anova(Mod_darkmagenta2, test = "F", type = "II") # same p-value for treatment
drop1(Mod_darkmagenta2, test = "F")

Mod_darkmagenta3 <- lm(MEdarkmagenta ~
                         Treatment, data = AllData)
summary(Mod_darkmagenta3) # R^2=-0.02729, p=0.96

anova(Mod_darkmagenta2, Mod_darkmagenta3) ***

Mod_white <- lm(MEwhite ~
                         Accession, data = AllData,
                       contrast=list(Accession=contr.sum))
summary(Mod_white) # R^2=0.4, p<0.0001

Mod_white2 <- lm(MEwhite ~
                  Accession + Treatment, data = AllData,
                contrast=list(Accession=contr.sum))
summary(Mod_white2) #R^2=0.7005***

Mod_white3 <- lm(MEwhite ~
                   Treatment, data = AllData)
summary(Mod_white3) 
# R ^2=0.2721***

# what coefficient is the R^2 value?
summ <- summary(Mod_white3)
summ$adj.r.squared
summary(Mod_white3)$adj.r.squared #this works too
summary(Mod_white3)$coefficients
summary(Mod_white3)$coefficients[,4]
anova(Mod_white3)[,5]
str(anova(Mod_white3))
anova(Mod_white3)[1,5] # this has the p-value

#### try rms package - did not work
lmod = ols(MEwhite ~ Accession + Treatment, data = AllData, x=TRUE, y=TRUE)
lmod

library(rms)
w <- plot(anova(Mod_white2), what='proportion R2', pl=FALSE)
