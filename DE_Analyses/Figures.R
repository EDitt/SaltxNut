
#########################
####### BAR CHART #######
#########################

viridis(5)
viridis(8)
brewer.pal(8, "Dark2")
display.brewer.pal(8, "Dark2")
show_col(pal_jco("default")(10))
show_col(viridis(n = 25))
show_col(viridis(n = 25, option = "E"))

#First attempt
#Opposite_col <- c("#8F7700FF")
#Opposite_col <- c("#AADC32FF") #from virids
#Opposite_col <- c("#E6AB02") #dark2
Opposite_col <- c("#3B3B3BFF") #2nd option from jco
Increased_col <- c("#A73030FF")
Decreased_col <- c("#EFC000FF")
NotDE_col <- c("#868686FF")
#NotDE_col <- c("#3B3B3BFF") # lighter grey so better contrast with army green color
Unchanged_col <- c("#3B528BFF") # used a darker blue from viridis
#Unchanged_col <- c("#003C67FF")

#Second (dark2)
Opposite_col <- c("#E6AB02") 
Increased_col <- c("#1B9E77")
Decreased_col <- c("#D95F02")
NotDE_col <- c("#666666")
Unchanged_col <- c("#3B528BFF") # used a darker blue from viridis

#Opposite_col <- c("#8F7700FF") 
#Opposite_col <- c("#3B3B3BFF") #dark brown
#Opposite_col <- c("darkkhaki") 

Opposite_col <- c("olivedrab4") 

Increased_col <- c("#A73030FF")
#Increased_col <- c("#CD534CFF") #lighter red
Decreased_col <- c("#EFC000FF")
#NotDE_col <- c("#868686FF")
NotDE_col <- c("grey55")
#Unchanged_col <- c("#3B528BFF") # used a darker blue from viridis
#Unchanged_col <- c("#0073C2FF") #bright blue
Unchanged_col <- c("#4A6990FF") #periwinkle blue

plot_ComboCats <- ggplot(df_all, aes(fill=labels, y=Number, x=Stress))

plot_ComboCats + geom_bar(position="stack", stat="identity", color="black", size=0.3) +
  ylab("Number of DE Genes") + xlab("Single Stress") +
  scale_fill_manual(values=c(Opposite_col, Increased_col, Decreased_col, NotDE_col, Unchanged_col,
                             Decreased_col, Unchanged_col,
                             Increased_col, Decreased_col, Unchanged_col,
                             Opposite_col, Increased_col, Decreased_col, Unchanged_col,
                             NotDE_col,
                             Opposite_col, Increased_col, Decreased_col, Unchanged_col),
                    #guide=guide_legend(reverse = TRUE),
                    name="Addition of Second Stress",
                    breaks=c("Unchanged", "Decreased_mag",
                             "Increased_mag", "Opposite_dir", "Not_DE"),
                    labels=c("Un-changed",
                             "Decreased Magnitude", "Increased Magnitude",
                             "Opposite Direction", "No Longer Differentially Expressed")) +
  theme_bw(base_size = 14)



x.seq <- c(1,2,4,5)
ggplot(data=transform(df_all, x=x.seq), aes(y=Number, x=Stress)) +
  geom_bar(position="stack", stat="identity", color="black", size=0.4, aes(fill=labels)) + 
  labs(x="", y="") + 
  scale_x_discrete(breaks = NA) + 
  geom_text(aes(x=c(sum(x.seq[1:2])/2, sum(x.seq[3:4])/2), y=0, 
                label=c("X","Y")), vjust=1.2, size=8)

Decreased_col <- c("#EFC000FF")
Decreased_col <- c("gold1")
Decreased_col <- c("#F8DF4BFF")
#Opposite_col <- c("khaki4")
#NotDE_col <- c("grey62")

test <- ggplot(df_all, aes(fill=labels, y=Number, x=Stress, dodge=Stress))
test + geom_col(position = "dodge")
, position_dodge(width = 0.5))

scale_x_discrete(breaks=1:4)
#scale_x_discrete(expand = c(0.5, 0.5)) +
scale_x_discrete(expand = waiver)
scale_x_discrete(expand_scale(mult = 0, add = 5))
scale_x_discrete("Stress", 
                 #breaks = c(1,3,5,6),
                 labels=c("Nutrient only", "Nutrient", "Salt", "Salt only"), position="bottom")
