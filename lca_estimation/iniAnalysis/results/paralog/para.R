setwd("/home/vinh/Desktop/data/project/iniAnalysis")
data <- read.table("/home/vinh/Desktop/data/project/iniAnalysis/orthomcl_micros_orig.list.para", sep='\t',header=T)

library(ggplot2)
library(RColorBrewer)

data$Paralogs <- ordered(data$Paralogs, levels = c("1","2","3-5","6-10","11-20",">20"))

png("/home/vinh/Desktop/data/project/iniAnalysis/orthomcl_micros_orig.list.para.png", width=1600, height=800)
p = ggplot(data, aes(x=Paralogs, y=Groups, group=Taxon, colour=Taxon)) +
  geom_line() +
  geom_point()
p = p+labs(x="Number of paralogs", y="Number of homologous groups")
p = p+theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),legend.text=element_text(size=15))
p = p + scale_fill_brewer(palette="Set2")
p
dev.off()
