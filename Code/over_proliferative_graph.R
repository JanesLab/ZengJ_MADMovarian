#install.packages("ggpubr")
#install.packages("plotrix") # Install plotrix R package
library("plotrix")   
library(ggpubr)
library(dplyr)
library(rstatix) 

#Read the data
setwd("/data/Clonal_distribution/")
over_prol <- read.csv(file = 'over_proliferative_WT.csv')
head(over_prol)

## Run the test and plot the resuls
stat.test <- aov(opc ~ month, data = over_prol) %>%
  tukey_hsd()

p <- ggerrorplot(over_prol, x = "month", y = "opc", add = "jitter", add.params = list(color = "darkgray", size=3), error.plot = "errorbar", xlab = "Age",
                 ylab = "Overproliferative clones (%)") +
                  scale_y_continuous(breaks = seq(0, 45, 5), 
                     limits = c(0,45), 
                     expand = c(0,0))
p + stat_summary(
  geom = "point",
  shape = 95,
  size = 10,
  col = "red",
  fun = "mean") + 
  stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  y.position = c(35, 38, 42)) 
ggsave('Overproliferative_clones.jpg', dpi=300, height=10, width = 12, units='cm')




