## figures for Module analyses

library(ggfortify)

pca_res <- prcomp(AllData[,lmTreatmentSigMods], center = TRUE, scale. = TRUE)
summary(pca_res) # PC1=35%, PC2=13%
autoplot(pca_res, data=AllData, 
         colour='Treatment', shape='Accession', frame = TRUE) +
  theme_minimal()
# looks remarkably similar to the one of all DE genes


