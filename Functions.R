#########################
##### DE GENE SETS ######
#########################

NumSigDE <- function(vector) {
  sig <- which(vector < 0.05)
  return (length(sig))
}

SigDE <- function(vector) {
  sig <- which(vector < 0.05)
  return (sig)
}

LessCritNum <- function(dataset, column, critNum) {
  less <- which(dataset[,column] < critNum)
  data_less <- dataset[less,]
  return (data_less)
}

MoreCritNum <- function(dataset, column, critNum) {
  more <- which(dataset[,column] > critNum)
  data_more <- dataset[more,]
  return (data_more)
}


VennNumbers <- function (dataset1, dataset2, dataset3) {
  name1 <- deparse(substitute(dataset1))
  name2 <- deparse(substitute(dataset2))
  name3 <- deparse(substitute(dataset3))
  sig1_2 <- intersect(dataset1, dataset2)
  sig2_3 <- intersect(dataset2, dataset3)
  sig1_3 <- intersect(dataset1, dataset3)
  sig1_2only <- setdiff(sig1_2, dataset3)
  sig2_3only <- setdiff(sig2_3, dataset1)
  sig1_3only <- setdiff(sig1_3, dataset2)
  all <- intersect(sig1_2, dataset3)
  sig1only <- setdiff(dataset1, union(dataset2, dataset3))
  sig2only <- setdiff(dataset2, union(dataset1, dataset3))
  sig3only <- setdiff(dataset3, union(dataset1, dataset2))
  #print(paste0("Common in all three = ", length(all)))
  #print(paste0("In ", name1, " and ", name2, " only = ", length(sig1_2only)))
  #print(paste0("In ", name2, " and ", name2, " only = ", length(sig2_3only)))
  #print(paste0("In ", name1, " and ", name3, " only = ", length(sig1_3only)))
  #print(paste0(name1, " only = ", length(sig1only)))
  #print(paste0(name2, " only = ", length(sig2only)))
  #print(paste0(name3, " only = ", length(sig3only)))
  my_list <- data.frame("Category" = c("InCommonAll", paste0(name1, name2, "Only"), paste0(name2, name3, "Only"), paste0(name1, name3, "Only"), paste0(name1, "Only"), paste0(name2, "Only"), paste0(name3, "Only")), 
                        "Number" = c(length(all), length(sig1_2only), length(sig2_3only), length(sig1_3only), length(sig1only), length(sig2only), length(sig3only)))
  #my_list <- data.frame("Category" = c("InCommonAll", paste0(names[1], names[2], "Only"), paste0(names[2], names[3], "Only"), paste0(names[1], names[3], "Only"), paste0(names[1], "Only"), paste0(names[2], "Only"), paste0(names[3], "Only")), 
  #"Number" = c(length(all), length(sig1_2only), length(sig2_3only), length(sig1_3only), length(sig1only), length(sig2only), length(sig3only)))
  return (my_list)
}

SigDEdf <- function(dataset, PvaluesCol, CritP) {
  sig <- which(dataset[,PvaluesCol] < CritP)
  Sigdf <- dataset[sig,]
  return (Sigdf)
}

GeneSets <- function (dataset1, dataset2, dataset3) {
  name1 <- deparse(substitute(dataset1))
  name2 <- deparse(substitute(dataset2))
  name3 <- deparse(substitute(dataset3))
  sig1_2 <- intersect(dataset1, dataset2)
  sig2_3 <- intersect(dataset2, dataset3)
  sig1_3 <- intersect(dataset1, dataset3)
  sig1_2only <- setdiff(sig1_2, dataset3)
  #name(sig1_2only) <- paste0(name1, name2, "Only")
  sig2_3only <- setdiff(sig2_3, dataset1)
  #name(sig2_3only) <- paste0(name2, name3, "Only")
  sig1_3only <- setdiff(sig1_3, dataset2)
  #name(sig1_3only) <- paste0(name1, name3, "Only") 
  InCommonAll <- intersect(sig1_2, dataset3)
  sig1only <- setdiff(dataset1, union(dataset2, dataset3))
  #name(sig1only) <- paste0(name1, "Only")
  sig2only <- setdiff(dataset2, union(dataset1, dataset3))
  #name(sig2only) <- paste0(name2, "Only")
  sig3only <- setdiff(dataset3, union(dataset1, dataset2))
  #name(sig3only) <- paste0(name3, "Only")
  Indices <- list(InCommonAll, sig1_2only, sig2_3only, sig1_3only, sig1only, sig2only, sig3only)
  names(Indices) <- c("InCommonAll", paste0(name1, name2, "Only"), paste0(name2, name3, "Only"), paste0(name1, name3, "Only"), paste0(name1, "Only"), paste0(name2, "Only"), paste0(name3, "Only"))
  return (Indices)
}

## Import files

ImportCSVs <- function (DirPath, critP) {
  my_files <- list.files(path = DirPath, pattern = "*.csv", full.names = TRUE)
  my_data <- lapply(my_files, read.csv)
  names(my_data) <- gsub("\\.csv$", "", my_files)
  names(my_data) <- gsub(DirPath, "", names(my_data))
  names(my_data) <- gsub("/", "", names(my_data))
  for (i in seq_along(my_data)) {
    name <- names(my_data[i])
    colnames(my_data[[i]])[1] <- "Gene"
    my_data[[i]] <- my_data[[i]][order(my_data[[i]]$Gene),]
  }
  return(my_data)
}

PairwiseCats <- function (Overlap, Diff, BothUp, BothDown) {
  Unconditional <- setdiff(Overlap, Diff)
  Up_Unconditional <- intersect(Unconditional, BothUp)
  Down_Unconditional <- intersect(Unconditional, BothDown)
  Conditional <- intersect(Overlap, Diff)
  Up_Conditional <- intersect(Conditional, BothUp)
  Down_Conditional <- intersect(Conditional, BothDown)
  DiffDir_Conditional <- setdiff(Conditional, union(Up_Conditional, Down_Conditional))
  my_list <- list(Unconditional_Up = Up_Unconditional,
                  Unconditional_Down = Down_Unconditional,
                  Conditional_Up = Up_Conditional,
                  Conditional_Down = Down_Conditional,
                  Conditional_Diff_Dir = DiffDir_Conditional)
  return(my_list)
}

Prioritized <- function (Diff, OverlapUp, BothUp_Col, OverlapDown, BothDown_Col) {
  Prioritized_Up <- lapply(OverlapUp, function(x) {intersect(x$Gene, Diff)})
  Prioritized_Down <- lapply(OverlapDown, function(x) {intersect(x$Gene, Diff)})
  my_list <- list(All_Up = Prioritized_Up[[1]],
                  Both_Up = Prioritized_Up[[BothUp_Col]],
                  All_Down = Prioritized_Down[[1]],
                  Both_Down = Prioritized_Down[[BothDown_Col]])
  return(my_list)
}

DirectionDf <- function(dataset, Yvar, Xvar) {
  Dir1 <- dataset[which(dataset[,Xvar] > 0 &
                          dataset[,Yvar] < 0),]
  Dir2 <- dataset[which(dataset[,Xvar] < 0 &
                          dataset[,Yvar] > 0),]
  DirGenes <- union(Dir1$Gene, Dir2$Gene)
  ReducedGenes <- setdiff(dataset[which(abs(dataset[,Xvar]) >
                                          abs(dataset[,Yvar])),"Gene"], 
                          DirGenes)
  ReducedData <- dataset[which(dataset$Gene %in% ReducedGenes),]
  IncreasedGenes <- setdiff(dataset[which(abs(dataset[,Xvar]) <
                                            abs(dataset[,Yvar])),"Gene"], 
                            DirGenes)
  IncreasedData <- dataset[which(dataset$Gene %in% IncreasedGenes),]
  Yname <- deparse(substitute(Yvar))
  Xname <- deparse(substitute(Xvar))
  MyList <- list("Down" = Dir1,
                 "Up" = Dir2,
                 "Reduced" = ReducedData,
                 "Increased" = IncreasedData)
  return (MyList)
}

#########################
##### REGRESSION #####
#########################

Predictdf <- function(dataset, geneset, Yvar, Xvar) {
  dfSubset <- dataset[which(dataset$Gene %in% geneset),]
  Mod <- lm(dfSubset[,Yvar] ~ dfSubset[,Xvar])
  sum<-summary(Mod)
  residSE <- sum$sigma
  dfPredict <- as.data.frame(predict(Mod, level=0.95, interval= 'prediction'))
  dfConfid <- as.data.frame(predict(Mod, level=0.95, interval= 'confidence'))
  colnames(dfConfid) <- c("CI_fit", "CI_lower", "CI_upper")
  resid <- resid(Mod)
  dfSubsetPredict <- cbind(dfSubset, resid, dfPredict, dfConfid)
  dfSubsetPredict$Interval <- ifelse(dfSubsetPredict[,Yvar] < 
                                       dfSubsetPredict$lwr |
                                       dfSubsetPredict[,Yvar] > 
                                       dfSubsetPredict$upr, "Outside", "Within")
  dfSubsetPredict$ResidSE_Interval <- ifelse(abs(dfSubsetPredict$resid) > 
                                               residSE, "Outside", "Within")
  dfSubsetPredict$CI <- ifelse(dfSubsetPredict[,Yvar] < 
                                       dfSubsetPredict$CI_lower |
                                       dfSubsetPredict[,Yvar] > 
                                       dfSubsetPredict$CI_upper, "Outside", "Within")
  return (dfSubsetPredict)
}

#########################
##### GO ENRICHMENT #####
#########################

### Function for GO Term Enrichment

#LengthTable must have Gene ID as first column + Length as 2nd column
#GO_Terms must have GeneID ("Parent") as first column, and GO Term ("Ontology_term") in 2nd column

GO_Enrichment <- function (AllGenes, DE_Genes, LengthTable, GO_Terms) {
  common <- intersect(AllGenes, LengthTable[,1])
  LengthTable <- as.data.frame(LengthTable[LengthTable[,1] %in% common, ]) #"all" genes
  LengthTable <- LengthTable[order(LengthTable[,1]),] #order by gene ID
  gene.vector=as.integer(AllGenes %in% DE_Genes)
  names(gene.vector)=AllGenes
  pwf <- nullp(gene.vector, bias.data = LengthTable[,2])
  GO.wall <- goseq(pwf, gene2cat = GO_Terms)
  enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.05]
  enrichedGO_table <- GO.wall[is.element(GO.wall$category, enriched.GO),]
  return(enrichedGO_table)
}

#########################
##### COEXPRESSION ######
#########################

Module_Cluster <- function(MEDissThres, datExpr, dynamicColors) {
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  mergedMEs = merge$newMEs
  datVA <- moduleEigengenes(datExpr, merge$colors)$varExplained
  df <- data.frame(matrix(unlist(datVA ), nrow=length(datVA), byrow=T))
  ModGeneNums <- as.data.frame(table(merge$colors))
  Metric <- c("NumMods", "MeanGeneNum", "MinGeneNum", "MaxGeneNum",
              "MeanVarExp", "MinVarExp", "MaxVarExp")
  Value <- c(length(mergedMEs), mean(ModGeneNums$Freq), min(ModGeneNums$Freq), max(ModGeneNums$Freq),
             mean(df[,1]), min(df[,1]), max(df[,1]))
  Result <- data.frame(Metric, Value)
  return(Result)
}