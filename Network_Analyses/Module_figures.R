## figures for Module analyses

library(ggplot2)
library(ggfortify)

pca_res <- prcomp(AllData[,lmTreatmentSigMods], center = TRUE, scale. = TRUE)
summary(pca_res) # PC1=35%, PC2=13%
autoplot(pca_res, data=AllData, 
         colour='Treatment', shape='Accession', frame = TRUE) +
  theme_minimal()
# looks remarkably similar to the one of all DE genes



#########################
#### PLOT CATEGORIES ####
#########################


# load ANOVA results
LM_Results <- read.csv("ResultsFiles/Coexpression/Module_Anova.csv", header=T) # 88 modules

## R^2 categories for simple models:
load("ResultsFiles/Coexpression/QC/Module_RsquaredCats.RData")

# categories:
load("ResultsFiles/Coexpression/LR_SigModOverlap.RData")
DE_Overlaps

# AllSigMods, SigModOverlap, DE_Overlaps, SigDiffOverlap, 
# LR_list, lm_list, DE_Combo_sig, DE_Salt_sig, DE_Nut_sig

### assign categories

LM_Results$Category <- ifelse(LM_Results$Module %in%ModuleRsquaredList$AccessionR2_more95,
                              "VariesByAccession", ifelse(LM_Results$Module %in% DE_Overlaps$InCommonAll,
                                                          "All", ifelse(LM_Results$Module %in% DE_Overlaps$DE_Combo_sigDE_Salt_sigOnly,
                                                                           "Salt_Combo", ifelse(LM_Results$Module %in% DE_Overlaps$DE_Salt_sigDE_Nut_sigOnly,
                                                                                                "Nut_Salt", ifelse(LM_Results$Module %in% DE_Overlaps$DE_Combo_sigDE_Nut_sigOnly,
                                                                                                                   "Nut_Combo", ifelse(LM_Results$Module %in% DE_Overlaps$DE_Salt_sigOnly,
                                                                                                                                       "Salt", ifelse(LM_Results$Module %in% DE_Overlaps$DE_Nut_sigOnly,
                                                                                                                                                      "Nutrient", 
                                                                                                                                                      "Other")))))))
#sig diff between Nut-Salt
LM_Results$Nut_SaltDiff <- ifelse(LM_Results$Difference_p.Nut.Salt < 0.05,
                                  "Sig. Different",
                                  "Not significantly different")

#sig diff between Nut-Combo
LM_Results$Nut_ComboDiff <- ifelse(LM_Results$Difference_p.Nut.Combo < 0.05,
                                  "Sig. Different",
                                  "Not significantly different")

LM_Results$Salt_ComboDiff <- ifelse(LM_Results$Difference_p.Salt.Combo < 0.05,
                                   "Sig. Different",
                                   "Not significantly different")


# check
aggregate(LM_Results$Module, by=list(LM_Results$Category), length)

plot(LM_Results$emmean.LowNut, LM_Results$emmean.Salt) # means
plot(LM_Results$Difference.DE_Nut, LM_Results$Difference.DE_Salt) # somewhat positive relationship when you look at differences from control

pmod <- ggplot(LM_Results, aes(emmean.LowNut, emmean.Salt))
pmod + geom_point()
pmod + geom_point(aes(fill=Category, color = Category, shape = Nut_SaltDiff), size = 3, stroke = 1)

# negative relationship with control or combo?
plot(LM_Results$emmean.LowNut, LM_Results$emmean.Control) # negative
plot(LM_Results$emmean.LowNut, LM_Results$emmean.Combo) # doesn't look negative
plot(LM_Results$emmean.Salt, LM_Results$emmean.Control) # looks positive
plot(LM_Results$emmean.Salt, LM_Results$emmean.Combo) # slightly negative?

pmod2 <- ggplot(LM_Results, aes(Difference.DE_Nut, Difference.DE_Salt))
pmod2 + geom_point(aes(fill=Category, color = Category, shape = Nut_SaltDiff), size = 3, stroke = 1) +
  geom_abline(intercept=0, slope= 1, linetype="dotted")

pmod3 <- ggplot(LM_Results, aes(Difference.DE_Nut, Difference.DE_Combo))
pmod3 + geom_point(aes(fill=Category, color = Category, shape = Nut_ComboDiff), size = 3, stroke = 1) +
  geom_abline(intercept=0, slope= 1, linetype="dotted")

pmod3 + geom_point(aes(fill=Difference_p.Nut.Combo), size = 3, stroke = 1) +
  geom_abline(intercept=0, slope= 1, linetype="dotted")

pmod4 <- ggplot(LM_Results, aes(Difference.DE_Salt, Difference.DE_Combo))
pmod4 + geom_point(aes(fill=Category, color=Category, shape = Salt_ComboDiff), size = 3, stroke = 1) +
  geom_abline(intercept=0, slope= 1, linetype="dotted") +
  scale_color_manual(values=c("red", "grey", "green", "grey", "grey", "blue", "orange", "grey"))

########

#pmod + geom_point(aes(fill=Category, shape=Category, color = Category), size = 3, stroke = 1) +
#  scale_shape_manual(values=c(19, 19, 24,
#                              19, 19, 19,
#                              24, 24,
#                              19, 19))

