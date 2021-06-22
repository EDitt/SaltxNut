# Graphing Module ID/Tree Cut parameters

library(ggplot2)
library(ggthemes)
library(cowplot)

### Graph Tree Cut Parameters

TreeCutParameter <- read.csv("Network_Analyses/PreProcess/TreeCutVary.csv", header=T)

p1 <- ggplot(data = TreeCutParameter, aes(NumMods, MeanVarExp_noMin))
p1 + geom_point(aes(color=pam)) + geom_line(aes(color=pam)) 
p1g <- p1 + geom_point(aes(color=pam)) + geom_line(aes(color=pam)) +
  ylab("Mean %Var Exp") + xlab("") + 
  scale_color_manual(values=c("#1B9E77", "#D95F02")) +
  theme_cowplot(font_size = 12) + theme(legend.position = c(0.7, 0.1),
                                        legend.direction = "horizontal", complete = FALSE)

viridis(2)
wes_palette("Darjeeling1", 2, type = c("discrete"))
brewer.pal(3, "Dark2")
display.brewer.pal(3, "Dark2")

p2 <- ggplot(data = TreeCutParameter, aes(NumMods, MeanGeneNum))
p2 + geom_point(aes(color=pam)) + geom_line(aes(color=pam)) +
  scale_color_manual(values=c("#1B9E77", "#D95F02"))
p2g <- p2 + geom_point(aes(color=pam)) + geom_line(aes(color=pam)) +
  ylab("Mean # Genes/Mod") + xlab("# of Modules") +
  scale_color_manual(values=c("#1B9E77", "#D95F02")) + theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.background=element_rect(fill=NULL)) 

plot_grid(p1g, p2g, ncol=1, nrow=2)

p3 <- ggplot(data = TreeCutParameter[which(TreeCutParameter$pam=="FALSE"),], aes(NumMods, NumUnAssigned))
p3 + geom_point(color="#1B9E77") + geom_line(color="#1B9E77")
p3g <- p3 + geom_point(color="#1B9E77") + geom_line(color="#1B9E77") + theme_cowplot(font_size = 12) +
  ylab("# Un-assigned") + xlab("")

TreeCutParameter$transform <- TreeCutParameter$NumUnAssigned / 25

# to transform data
coeff <- 250

p3g <- p2g +  geom_line(data = TreeCutParameter[which(TreeCutParameter$pam=="FALSE"),], aes(y=transform), color="#1B9E77", lty=2) + 
  geom_point(data = TreeCutParameter[which(TreeCutParameter$pam=="FALSE"),], aes(y=transform), color="#1B9E77") +
  scale_y_continuous(sec.axis = sec_axis(~ . * 25, name="# Un-assigned")) 

plot_grid(p1g, p3g, align="hv", ncol=1, nrow=2)
ggsave(filename = "ResultsFiles/Coexpression/ModuleIDparameters.png")

### Graph Cut Height for Merge

ModuleClustParameter <- read.csv("Network_Analyses/PreProcess/VariableCutThreshold_PAMFALSE_split2.csv", header=T)

head(ModuleClustParameter)

p4 <- ggplot(data=ModuleClustParameter, aes(x=Threshold, y=MeanVarExp))

p4g <- p4 + geom_point() + geom_line() + theme_cowplot(font_size = 12) + ylim(0,1.1) + 
  ylab("Mean %Var Exp") + xlab("Cut Height")

ModuleClustParameter$ConvertModNum <- ModuleClustParameter$ModuleNumber / 100

p4g + geom_line(data = ModuleClustParameter, aes(y=ConvertModNum), lty=2) +
  geom_point(data = ModuleClustParameter, aes(y=ConvertModNum), pch=21) +
  scale_y_continuous(sec.axis = sec_axis(~ .*100, name= "Number of Modules"))

ggsave(filename = "ResultsFiles/Coexpression/ModuleClustCut.png")


