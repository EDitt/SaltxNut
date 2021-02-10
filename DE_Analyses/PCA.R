#srun --pty  -p inter_p  --mem=22G --nodes=1 --ntasks-per-node=1 --time=12:00:00 --job-name=qlogin /bin/bash -l
#module load R/3.6.2-foss-2019b
#R

#########################
######## SETUP ##########
#########################

library(DESeq2)
library(tximport)

setwd("/scratch/eld72413/Salty_Nut/CultivatedOnly/DE_Analyses_Inbred")

library( "genefilter" )
library("gplots")
library("RColorBrewer")

#########################
######### DATA ##########
#########################

#load("ddsTxi_cult_red.RData")

load("dds_cult.RData")

# df "Sig_AllTrt" of genes significant in at least 1 treatment
load("SigGenes.RData")

# transform
vsd <- vst(dds_cult, blind = FALSE)

#########################
######## ORDER? #########
#########################

# top variable genes
topVarGenes <- head( order( rowVars( assay(vsd) ), decreasing=TRUE ), 35 )

# all significant genes
SigGenes <- assay(vsd[which(rownames(vsd) %in% Sig_AllTrt$Gene),])
topSigGenes <- head( order( rowVars( SigGenes ), decreasing=TRUE ), 35 )


#########################
####### HEATMAP #########
#########################

pdf(file = "Heatmap.pdf", width = 12, height = 12);
heatmap.2( assay(vsd)[ topVarGenes, ], scale="row",
trace="none", dendrogram="column",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()

pdf(file = "HeatmapSig.pdf", width = 12, height = 12);
heatmap.2( head(SigGenes,35), scale="row",
trace="none", dendrogram="column",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()

pdf(file = "HeatmapSig2.pdf", width = 12, height = 12);
heatmap.2( SigGenes[ topSigGenes, ], scale="row",
trace="none", dendrogram="column",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()

pdf(file = "HeatmapSigAll.pdf", width = 12, height = 12);
heatmap.2( SigGenes, scale="row",
trace="none", dendrogram="column",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()
