### GO TERM ENRICHMENT

### Using GFF3 file filtered for mRNA

#bash code:
#awk '{if ($3 == "mRNA") {print $9}}' Ha412HOv2.0-20181130.gff3 | awk -F "[;=]" '{$1 = $1; print $2,$12}' > ./Ha412HOv2_GOTerms.txt

#########################
##### DEPENDENCIES ######
#########################

library(tidyr)
library(goseq)
library(tximport)
library(matrixStats)
source("Functions.R")

############################
##### FORMAT FOR GOSEQ #####
############################

HA412_Genes <- read.table("DataFiles/Ha412HOv2_GOTerms.txt", fill = T)
colnames(HA412_Genes) <- c("Parent", "Ontology_term")
head(HA412_Genes)
str(HA412_Genes)

### Need to make a GO Term .txt file with only 1 line per gene
HA412_Genes$Ontology_term <- as.character(HA412_Genes$Ontology_term)
goMapping <- separate_rows(HA412_Genes, Ontology_term, sep = ",")

#save
save(goMapping, file = "DataFiles/goMapping.RData")
#GO_Terms must have GeneID ("Parent") as first column, and GO Term ("Ontology_term") in 2nd column

############################
## TRANSCRIPT LENGTH INFO ##
############################

load("txi.rsem.cult.RData")

GeneLength <- txi.rsem.cult$length
Medians <- rowMedians(GeneLength)
Median_Lengths <- cbind(row.names(GeneLength), Medians)
colnames(Median_Lengths) <- c("Gene","Median")
write.table(Median_Lengths, file = "GeneLengthMedians")

############################
######## READ FILES ########
############################

### load in files: transcript length (median value), list of DE genes, parsed GO Terms

Length_table <- read.table("DataFiles/GeneLengthMedians", header=T)
dim(Length_table) #78260
head(Length_table)
str(Length_table)
#LengthTable must have Gene ID as first column + Length as 2nd column


# need all genes being tested, read in one of the results files
result.file <- read.csv("ResultsFiles/GeneSets/condition_comboDE_results.csv", header=T)
AllGenes <- result.file$X
length(AllGenes) # 55728


############################
## GO ENRICHMENT ANALYSES ##
############################

### Using function "GO_Enrichment" from Functions.R script (relies on package goseq)
load("DataFiles/goMapping.RData") #goMapping object

load("ResultsFiles/GeneSets/MultiStressCompare.RData")
# This has 4 lists: NutSpecific_Cats, SaltSpecific_Cats, NutUnSpecific_Cats, SaltUnSpecific_Cats
# Within each list are genes separated by categories based on how they compare to combination treatment


############################
#### NUTRIENT-SPECIFIC #####
############################

lapply(NutSpecific_Cats, function(x) {length(x)})

GO_NutSpec <- lapply(NutSpecific_Cats, function(x) 
  {GO_Enrichment(AllGenes, x, Length_table, goMapping)})
lapply(GO_NutSpec, function(x) {length(x$category)}) 
#21 for nutrient only, 5 for un-changed, 29 for reduced in combo, 0 for combo increased, 0 for combo in opposite direction

###### Take direction into account?
nutDE <- read.csv("ResultsFiles/GeneSets/condition_nutDE_results.csv")
SigUpNut <- subset(nutDE, padj < 0.05 & log2FoldChange > 0) # 10092
SigDownNut <- subset(nutDE, padj < 0.05 & log2FoldChange < 0) # 9428

NutSpec_Up <- lapply(NutSpecific_Cats, function(x) {
  intersect(SigUpNut$X, x)
})
lapply(NutSpec_Up, function(x) {length(x)}) # 4628 nut. only, 1475 unchanged, 738 reduced in combo, 5 increased in combo, 14 diff direction

GO_NutSpecUP <- lapply(NutSpec_Up, function(x) 
{GO_Enrichment(AllGenes, x, Length_table, goMapping)})

lapply(GO_NutSpecUP, function(x) {length(x$category)}) # 21 nut. only, 4 unchanged, 6 reduced in combo

NutSpec_Down <- lapply(NutSpecific_Cats, function(x) {
  intersect(SigDownNut$X, x)
})
lapply(NutSpec_Down, function(x) {length(x)})

GO_NutSpecDOWN <- lapply(NutSpec_Down, function(x) 
{GO_Enrichment(AllGenes, x, Length_table, goMapping)})

lapply(GO_NutSpecDOWN, function(x) {length(x$category)}) # 31 nut. only, 25 unchanged, 80 reduced

############################
####### SALT-SPECIFIC ######
############################

lapply(SaltSpecific_Cats, function(x) {length(x)})

GO_SaltSpec <- lapply(SaltSpecific_Cats, function(x) 
{GO_Enrichment(AllGenes, x, Length_table, goMapping)})
lapply(GO_SaltSpec, function(x) {length(x$category)}) 
# 5 for salt only, 1 for combo reduced, 0 for everything else

###### Do I want to take into account whether they are up or down-regulated?
saltDE <- read.csv("ResultsFiles/GeneSets/condition_saltDE_results.csv")

SigUpSalt <- subset(saltDE, padj < 0.05 & log2FoldChange > 0) #4974
SigDownSalt <- subset(saltDE, padj < 0.05 & log2FoldChange < 0) #5143

SaltSpec_Up <- lapply(SaltSpecific_Cats, function(x) {
  intersect(SigUpSalt$X, x)
})
lapply(SaltSpec_Up, function(x) {length(x)}) 
#1110 salt only, 630 unchanged, 24 reduced in combo, 1 increased in combo

GO_SaltSpecUP <- lapply(SaltSpec_Up, function(x) 
{GO_Enrichment(AllGenes, x, Length_table, goMapping)})

lapply(GO_SaltSpecUP, function(x) {length(x$category)}) 
# 16 for salt only (5 of which were found for the analysis without regard to direction)

### down-regulated
SaltSpec_Down <- lapply(SaltSpecific_Cats, function(x) {
  intersect(SigDownSalt$X, x)
})
lapply(SaltSpec_Down, function(x) {length(x)}) #908 salt only, 512 unchanged 38 reduced in combo, 7 diff direction

GO_SaltSpec_DOWN <- lapply(SaltSpec_Down, function(x) 
{GO_Enrichment(AllGenes, x, Length_table, goMapping)})

lapply(GO_SaltSpec_DOWN, function(x) {length(x$category)}) # 1 in salt only (not found in other analysis), 1 in combo reduced (same as found in other category)


############################
######## UNSPECIFIC ########
############################

lapply(NutUnSpecific_Cats, function(x) {length(x)})

SingleStressOnly <- GO_Enrichment(AllGenes, NutUnSpecific_Cats$Nutrient_only, 
                                  Length_table, goMapping)
length(SingleStressOnly$term) #N=3

Unchanged <- GO_Enrichment(AllGenes, NutUnSpecific_Cats$Combo_Nutrient_Same, 
                                  Length_table, goMapping)
length(Unchanged$term) 


