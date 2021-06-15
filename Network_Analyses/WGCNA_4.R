
# To obtain information about each gene

load("Consensus-NetworkConstruction-man.RData")
# objects: consMEs, moduleLabels, consTree

load("multiExpr.RData")
nSets = checkSets(multiExpr)$nSets

#############################
# GET GENE MODULE ASSIGNMENTS 
#############################

length(colnames(multiExpr[[1]]$data)) #34,408
length(moduleLabels) #34,408
modGenes <- as.data.frame(cbind(moduleLabels, colnames(multiExpr[[1]]$data)))
colnames(modGenes)[2] <- c("Gene")
#check
aggregate(modGenes$Gene, by=list(modGenes$moduleLabels), length)
table(moduleLabels)

################################
###### MODULE MEMBERSHIP #######
################################
### For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile.
# this allows us to quantify the similarity of all genes to every module

# gene names
probes <- colnames(multiExpr[[1]]$data)

# recalculate module eigengenes in the "alphabetic" order and calculate gene significances and module memberships in each dataset
# I don't understand why ^ they do this, so I will take original module assignments

kME = list();
for (set in 1:nSets)
{
kME[[set]] = corAndPvalue(multiExpr[[set]]$data, consMEs[[set]]$data);
}

# combine the Z scores of correlations from each set to form a "meta-Z" score and corresponding p-value
kME.metaZ = (kME[[1]]$Z + kME[[2]]$Z)/sqrt(2)
kME.metaP = 2*pnorm(abs(kME.metaZ), lower.tail = FALSE)

# form matrices to hold the kME. Re-shaping trick to put the values and associated p-values and meta-analysis results next to each other
kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[1]]$p, kME[[2]]$p, kME.metaZ, kME.metaP);
MEnames = colnames(consMEs[[1]]$data);
nMEs = checkSets(consMEs)$nGenes
dim(kMEmat) = c(nGenes, 6*nMEs)
rownames(kMEmat) = probes;
colnames(kMEmat) = spaste(
c("kME.set1.", "kME.set2.", "p.kME.set1.", "p.kME.set2.", "Z.kME.meta.", "p.kME.meta"),
rep(MEnames, rep(6, nMEs)))

# put together full information and write into a plain text CSV file (probes not sorted in any particular way)

info = data.frame(Probe = probes,
ModuleLabel = moduleLabels,
kMEmat)

head(info[,1:3])
#check
aggregate(info$Probe, by=list(info$ModuleLabel), length)
table(moduleLabels)

write.csv(info, file = "Gene_ModuleMembership.csv", row.names = FALSE, quote = FALSE)

#############################
##### GENE CONNECTIVITY #####
#############################

#The function intramodularConnectivity computes:

# ktotal = whole network connectivity
# kWithin = within module connectivity
# kOut = kTotal - kWithin
# kDiff = kIn - kOut = 2*kIn-kTotal

Alldegrees1 = list()
for (set in 1:nSets)
{
Alldegrees1[[set]] = intramodularConnectivity.fromExpr(multiExpr[[set]]$data, moduleLabels,
	networkType = "signed", ignoreColors = "grey",
	getWholeNetworkConnectivity = TRUE)
}


# softConnectivity: FYI: connecitivty of genes with less than 17 valid samples will be returned as NA.
# ..calculating connectivities....100% 
# softConnectivity: FYI: connecitivty of genes with less than 18 valid samples will be returned as NA.
# ..calculating connectivities....100% 

# how to summarize across the 2 datasets?

colnames(Alldegrees1[[1]]) <- paste0("HA_", colnames(Alldegrees1[[1]]))
colnames(Alldegrees1[[2]]) <- paste0("RHA_", colnames(Alldegrees1[[2]]))
All_Connectivity <- cbind(probes, Alldegrees1[[1]], Alldegrees1[[2]])

write.csv(All_Connectivity, 'Gene_Connectivity.csv', row.names = FALSE, quote = FALSE)

#############################
####### TOP HUB GENES ####### 
#############################

softPower = 12

top_hubs = list()
for (set in 1:nSets)
{
top_hubs[[set]] <- chooseTopHubInEachModule(
multiExpr[[set]]$data,
moduleLabels,
omitColors = "grey",
power = softPower,
type = "signed",
corFnc = bicor)
}

names(top_hubs[[1]])
top_hub_1 <- data.frame(matrix(unlist(top_hubs[[1]]), nrow=length(top_hubs[[1]]), byrow=T), stringsAsFactors=FALSE)
colnames(top_hub_1) <- "TopHub_HA"
top_hub_1$Modules <- names(top_hubs[[1]])

top_hub_2 <- data.frame(matrix(unlist(top_hubs[[2]]), nrow=length(top_hubs[[2]]), byrow=T), stringsAsFactors=FALSE)
colnames(top_hub_2) <- "TopHub_RHA"
top_hub_2$Modules <- names(top_hubs[[2]])

top_hubs_all <- merge(top_hub_1, top_hub_2, by="Modules")

write.csv(top_hubs_all , file = "topHubGenes.csv", row.names = FALSE, quote = FALSE)

