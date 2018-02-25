setwd("/home/vinh/Desktop/data/project/iniAnalysis/results/lengthTest")
data <- read.table("/home/vinh/Desktop/data/project/iniAnalysis/results/lengthTest/Pfam_homologous_vs_orphan.LEN.list", sep='\t',header=T)

library(ggplot2)
library(reshape2) # for melt
library(RColorBrewer)

data$taxon <- ordered(data$taxon, levels = c("E.hellem","E.intestinallis","E.cuniculi","N.ceranae","E.bieneusi","V.corneae","A.algerae","A.locustae","E.aedis","V.culicis","N.parisii"))

png("/home/vinh/Desktop/data/project/iniAnalysis/results/lengthTest/Pfam_homologous_vs_orphan.LEN.png", width=800, height=800)

p = ggplot(data, aes(x=type, y=length,fill=type))+
    geom_boxplot()+
    facet_wrap(~taxon,ncol=6,nrow=2)+
    theme_set(theme_gray(base_size = 18)
)

p = p+labs(x="", y="Length")
p = p+theme(axis.text.x = element_blank(),axis.text.y = element_text(size=15),strip.text.x = element_text(size = 18),legend.text=element_text(size=15))
p = p + scale_fill_brewer(palette="Set2")
p
dev.off()
