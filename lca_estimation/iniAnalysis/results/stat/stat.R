setwd("/home/vinh/Desktop/data/project/iniAnalysis")
data <- read.table("/home/vinh/Desktop/data/project/iniAnalysis/orthomcl_micros_orig.list.stat", sep='\t',header=T)

library(ggplot2)
library(reshape2) # for melt
library(RColorBrewer)

data$Taxon <- ordered(data$Taxon, levels = c("E.hellem","E.intestinallis","E.cuniculi","N.ceranae","E.bieneusi","V.corneae","A.algerae","A.locustae","E.aedis","V.culicis","N.parisii"))
melted <- melt(data, "Taxon")

png("/home/vinh/Desktop/data/project/iniAnalysis/orthomcl_micros_orig.list.stat.png", width=1600, height=800)
melted$cat <- ''
melted[melted$variable == 'LCA_homologous_prots',]$cat <- "#Genes"
melted[melted$variable == 'Non.LCA_homologous_prots',]$cat <- "#Genes"
melted[melted$variable == 'Orphan_prots',]$cat <- "#Genes"
melted[melted$variable == 'LCA.OG',]$cat <- "#Groups"
melted[melted$variable == 'Non.LCA.OG',]$cat <- "#Groups"

p = ggplot(melted, aes(x = cat, y = value, fill = variable)) + 
    geom_bar(stat = 'identity', position = 'stack') + 
    facet_grid(~ Taxon) +
    theme_set(theme_gray(base_size = 15))
p = p+labs(x="", y="")
p = p+theme(axis.text.x = element_blank(),axis.text.y = element_text(size=15),strip.text.x = element_text(size = 18),legend.text=element_text(size=15)) + 
    theme(legend.title=element_blank())
p = p + scale_fill_brewer(palette="Set2")
p
dev.off()
