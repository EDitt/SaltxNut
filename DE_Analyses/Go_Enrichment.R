### GO TERM ENRICHMENT

### Using GFF3 file filtered for mRNA

#bash code:
#awk '{if ($3 == "mRNA") {print $9}}' Ha412HOv2.0-20181130.gff3 | awk -F "[;=]" '{$1 = $1; print $2,$12}' > ./Ha412HOv2_GOTerms.txt

#########################
######## SETUP ##########
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

load("DataFiles/goMapping.RData") #goMapping object

# need all genes being tested, read in one of the results files
result.file <- read.csv("ResultsFiles/GeneSets/condition_comboDE_results.csv", header=T)
AllGenes <- result.file$X
length(AllGenes) # 55728


############################
## GO ENRICHMENT ANALYSES ##
############################

### Using function "GO_Enrichment" from Functions.R script (relies on package goseq)

load("MultiStressCompare.RData")
# This has 4 lists: NutSpecific_Cats, SaltSpecific_Cats, NutUnSpecific_Cats, SaltUnSpecific_Cats
# Within each list are genes separated by categories based on how they compare to combination treatment


