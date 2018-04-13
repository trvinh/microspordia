setwd("/home/vinh/Desktop/data/project/iniAnalysis")
data <- read.table("/home/vinh/Desktop/data/project/iniAnalysis/length_mean.stat", sep='\t',header=T)

library(ggplot2)
library(RColorBrewer)

data$Taxon <- ordered(data$Taxon, levels = c("E.hellem","E.intestinallis","E.cuniculi","N.ceranae","E.bieneusi","V.corneae","A.algerae","A.locustae","E.aedis","V.culicis","N.parisii"))

png("/home/vinh/Desktop/data/project/iniAnalysis/length_mean.stat.png", width=1600, height=800)
p = ggplot(data, aes(x=Taxon, y=meanLength, group=Type, colour=Type)) +
  geom_line(colour = "gray90") +
  geom_point(size=3)
p = p+labs(x="", y="Sequence length (mean)")
p = p+theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),legend.text=element_text(size=15)) + 
    theme(legend.title=element_blank())
p = p + scale_fill_brewer(palette="Set2")
p
dev.off()
