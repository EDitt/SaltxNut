#qsub -I -q s_interq -l walltime=4:00:00 -l nodes=1:ppn=8 -l mem=50gb
#module load R/3.6.1-foss-2018a-X11-20180131-GACRC
#R

#########################
######## SETUP ##########
#########################

library(DESeq2)
library(tximport)

setwd("/scratch/eld72413/Salty_Nut/CultivatedOnly/DE_Analyses_Inbred")

### To parallelize:
library("BiocParallel")
register(MulticoreParam(8)) #register cores so can specify parallel=TRUE when need to parallelize. Need to increase memory if increasing cores.
#here I used 8 cores and 50 gb memory/core

#############################
### LOAD SAMPLE METADATA ####
#############################

### Experimental Data
#open up grouping file
experiment_InbrednoOut <- read.csv("StudyDesign_Inbred_noOut.csv", header=T)
head(experiment_InbrednoOut)

#############################
# IMPORT DATA WITH TXIMPORT #
#############################

# create vector pointing to quantification files - vector of filenames that contains sample IDs and combine with
files_cult <- file.path("/scratch/eld72413/Salty_Nut/CultivatedOnly/DE_Analyses/RSEMOut_Files", paste0("RSEMOut_", experiment_InbrednoOut$Plant, ".merged_.genes.results"))
names(files_cult) <- experiment_InbrednoOut$Plant #assign sampleIDs

all(file.exists(files_cult)) #make sure files actually exist (=TRUE)

txi.rsem.cult <- tximport(files_cult, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem.cult$counts)
dim(txi.rsem.cult$counts) #78260 x 104

#Save this R object
save(txi.rsem.cult, file = "txi.rsem.cult_InbrednoOut.RData")


###############################
# CONSTRUCT DESEQ DATA OBJECT #
###############################

### Variables
experiment_InbrednoOut$SampleDay<-factor(experiment_InbrednoOut$SampleDay) #day sampled 1-2
experiment_InbrednoOut$Cross<-factor(experiment_InbrednoOut$Cross) #1,2,3

### Variables of interest:
levels(experiment_InbrednoOut$Treatment)
experiment_InbrednoOut$Treatment <- relevel(experiment_InbrednoOut$Treatment, ref="Control") #relevel to make Control reference
experiment_InbrednoOut$Accession<-as.factor(make.names(experiment_InbrednoOut$Accession))  #make syntactically valid names (no characters R doesn't like)
levels(experiment_InbrednoOut$Accession)
levels(experiment_InbrednoOut$Group) #HA, RHA

### Other variables to account for:
levels(experiment_InbrednoOut$Reproductive) #whether plant was reproductive when sampled


## make DESeq2 data object:
ddsTxi_cult <- DESeqDataSetFromTximport(txi.rsem.cult,
                                        colData = experiment_InbrednoOut,
                                        design = ~ 0 + Group + Group:Cross + SampleDay + Reproductive + Treatment)
ddsTxi_cult #dim is 78260 x 104

### Pre-filtering to reduce size of data object:
keep <- rowSums(counts(ddsTxi_cult) >= 1) >= 3 ### at least 3 samples have a count of 1 or higher
length(which(keep==1)) #55,728
ddsTxi_cult_red <- ddsTxi_cult[keep,]
dim(ddsTxi_cult_red) #55728 x 104

save(ddsTxi_cult_red, file = "ddsTxi_cult_red.RData")


###############################
######## DE ANALYSES ##########
###############################

### Need to parallelize due to complex design and many samples (use parallel = TRUE)

# run DESEQ2 analysis:
dds_cult <- DESeq(ddsTxi_cult_red, parallel = TRUE) 
# this does all steps: 
# estimate size factors (to control for library size diffs) - here it's using info from avgTxLength
# estimate dispersion for each gene
# fit glm
# returns DESeqDataSet which contains all the fitted info within it

# save as R object
save(dds_cult, file = "dds_cult.RData")

###############################
###### EXTRACT RESULTS ########
###############################

#load("dds_cult.RData")
resultsNames(dds_cult) ## the elements of results to be extracted

results_Combo <- results(dds_cult, alpha = 0.05, name="TreatmentCombo")
results_Salt <- results(dds_cult, alpha = 0.05, name="TreatmentHighSalt")
results_Nut <- results(dds_cult, alpha = 0.05, name="TreatmentLowNut")

#significant differences between each treatment and combo-
# these take longer and should parallelize to run
results_Salt_vs_Combo <- results(dds_cult, contrast=c("Treatment","HighSalt","Combo"), alpha = 0.05, parallel=TRUE)
results_Nut_vs_Combo <- results(dds_cult, contrast=c("Treatment","LowNut","Combo"), alpha = 0.05, parallel=TRUE)
results_Nut_vs_Salt <- results(dds_cult, contrast=c("Treatment","LowNut","HighSalt"), alpha = 0.05, parallel=TRUE)

### summarize basic tallies
summary(results_Salt_vs_Combo) #2.6 % up, 3.2% down
summary(results_Nut_vs_Combo) #9.6% up, 11% down

write.csv(as.data.frame(results_Combo), 
          file="condition_comboDE_results.csv")

write.csv(as.data.frame(results_Salt), 
          file="condition_saltDE_results.csv")

write.csv(as.data.frame(results_Nut), 
          file="condition_nutDE_results.csv")

write.csv(as.data.frame(results_Salt_vs_Combo), 
          file="Salt_vs_Combo_results.csv")

write.csv(as.data.frame(results_Nut_vs_Combo), 
          file="Nut_vs_Combo_results.csv")

write.csv(as.data.frame(results_Nut_vs_Salt), 
          file="Nut_vs_Salt_results.csv")
