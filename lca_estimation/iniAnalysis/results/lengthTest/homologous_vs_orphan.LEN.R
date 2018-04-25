setwd("/Users/trvinh/work/thesis/microsporidia/lca_estimation/iniAnalysis/results/lengthTest")
data <- read.table("/Users/trvinh/work/thesis/microsporidia/lca_estimation/iniAnalysis/results/lengthTest/homologous_vs_orphan.LEN.list", sep='\t',header=T)

library(ggplot2)
library(reshape2) # for melt
library(RColorBrewer)

data$taxon <- ordered(data$taxon, levels = c("E.hellem","E.intestinallis","E.cuniculi","N.ceranae","E.bieneusi","V.corneae","A.algerae","A.locustae","E.aedis","V.culicis","N.parisii"))
data$length <- log(data$length)

head(data)
data$type <- as.character(data$type)
data$type

data[data$type == "homologous",]$type <- "orthologous protein"

png("/Users/trvinh/work/thesis/microsporidia/lca_estimation/iniAnalysis/results/lengthTest/final_homologous_vs_orphan.LEN.png", width=800, height=500)

p = ggplot(data, aes(x=type, y=length,fill=type))+
    geom_boxplot()+
    facet_wrap(~taxon,ncol=6,nrow=2)+
    theme_set(theme_gray(base_size = 18)
    ) +
  labs(x="", y="log(length)") + 
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text.x = element_text(size = 17),
        legend.text=element_text(size=15),
        legend.title = element_blank(),
        legend.position = "top") + 
  scale_fill_brewer(palette="Set2")
p
dev.off()
