#### Experiment with parameters

library(ggplot2)
library(ggthemes)
library(cowplot)

viridis(2)
wes_palette("Darjeeling1", 2, type = c("discrete"))
brewer.pal(3, "Dark2")
display.brewer.pal(3, "Dark2")

#################################
### Graph Tree Cut Parameters ###
#################################

TreeCutParameter <- read.csv("Network_Analyses/PreProcess/TreeCutVary.csv", header=T)

p1 <- ggplot(data = TreeCutParameter, aes(NumMods, MeanVarExp_noMin))
p1 + geom_point(aes(color=pam)) + geom_line(aes(color=pam)) 
p1g <- p1 + geom_point(aes(color=pam)) + geom_line(aes(color=pam)) +
  ylab("Mean %Var Exp") + xlab("Number of Modules") + 
  scale_color_manual(values=c("#1B9E77", "#D95F02")) +
  theme_cowplot(font_size = 12) + theme(legend.position = c(0.5, 0.1),
                                        legend.direction = "horizontal", complete = FALSE) 
#+ geom_errorbar(aes(ymin=MinVarExp, ymax=MaxVarExp))


p2 <- ggplot(data = TreeCutParameter, aes(NumMods, MeanGeneNum))
p2 + geom_point(aes(color=pam)) + geom_line(aes(color=pam)) +
  scale_color_manual(values=c("#1B9E77", "#D95F02"))
p2g <- p2 + geom_point(aes(color=pam)) + geom_line(aes(color=pam)) +
  ylab("Mean # Genes/Mod") + xlab("# of Modules") +
  scale_color_manual(values=c("#1B9E77", "#D95F02")) + theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.background=element_rect(fill=NULL)) 

plot_grid(p1g, p2g, ncol=1, nrow=2)

# add in # un-assigned genes for pam=FALSE
p3 <- ggplot(data = TreeCutParameter[which(TreeCutParameter$pam=="FALSE"),], aes(NumMods, NumUnAssigned))
p3 + geom_point(color="#1B9E77") + geom_line(color="#1B9E77")
p3g <- p3 + geom_point(color="#1B9E77") + geom_line(color="#1B9E77") + theme_cowplot(font_size = 12) +
  ylab("# Un-assigned") + xlab("")

TreeCutParameter$transform <- TreeCutParameter$NumUnAssigned / 15

# to transform data
coeff <- 200
p3 + geom_point(color="#1B9E77") + geom_line(color="#1B9E77") 

p3g <- p2g +  geom_line(data = TreeCutParameter[which(TreeCutParameter$pam=="FALSE"),], aes(y=transform), color="#1B9E77", lty=2) + 
  geom_point(data = TreeCutParameter[which(TreeCutParameter$pam=="FALSE"),], aes(y=transform), color="#1B9E77") +
  scale_y_continuous(sec.axis = sec_axis(~ . *15, name="# Un-assigned")) 

plot_grid(p1g, p3g, align="hv", ncol=1, nrow=2)
ggsave(filename = "Network_Analyses/PreProcess/ModuleIDparameters.png")


###################################
# Graph Module Cluster Parameters #
###################################

### For pam=TRUE module
ModuleClustParameter <- read.csv("Network_Analyses/PreProcess/VariableCutThreshold_PAMTRUE.csv", header=T)

head(ModuleClustParameter)

p4 <- ggplot(data=ModuleClustParameter, aes(x=Threshold, y=MeanVarExp))

p4g <- p4 + geom_point() + geom_line() + theme_cowplot(font_size = 12) + ylim(0,1.1) + 
  ylab("Mean %Var Exp") + xlab("Cut Height")

ModuleClustParameter$ConvertGeneNum <- ModuleClustParameter$MeanNumGenes / 1000

p4g + geom_line(data = ModuleClustParameter[which(ModuleClustParameter$ConvertGeneNum < 1),], aes(y=ConvertGeneNum), lty=2) +
  geom_point(data = ModuleClustParameter[which(ModuleClustParameter$ConvertGeneNum < 1),], aes(y=ConvertGeneNum), pch=21) +
  scale_y_continuous(sec.axis = sec_axis(~ .*1000, name= "Mean # Genes/ Module"))

ModuleClustParameter$ConvertModNum <- ModuleClustParameter$ModuleNumber / 100

p4g + geom_line(data = ModuleClustParameter, aes(y=ConvertModNum), lty=2) +
  geom_point(data = ModuleClustParameter, aes(y=ConvertModNum), pch=21) +
  scale_y_continuous(sec.axis = sec_axis(~ .*100, name= "Num Modules"))


ggsave(filename = "Network_Analyses/PreProcess/ModuleClustCut.png")
