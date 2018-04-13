setwd("/Users/trvinh/work/thesis/microsporidia/lca_estimation/iniAnalysis/results/paralog")
data <- read.table("orthomcl_micros_orig.list.para", sep='\t',header=T)

library(ggplot2)
library(RColorBrewer)

data$Paralogs <- ordered(data$Paralogs, levels = c("1","2","3-5","6-10","11-20",">20"))
data$Taxon <- ordered(data$Taxon, 
                      # levels = c("E.intestinalis","E.hellem","E.cuniculi","N.ceranae","E.bieneusi","V.corneae","A.locustae","A.algerae","V.culicis","E.aedis","N.parisii"))
                      levels = c("E.hellem","E.intestinalis","E.cuniculi","N.ceranae","E.bieneusi","V.corneae","A.algerae","A.locustae","E.aedis","V.culicis","N.parisii"))

# png("/home/vinh/Desktop/data/project/iniAnalysis/orthomcl_micros_orig.list.para.png", width=1600, height=800)
p = ggplot(data, aes(x=Paralogs, y=Groups, group=Taxon, colour=Taxon)) +
  geom_line() +
  geom_point()
p = p+labs(x="Number of paralogs", y="Number of homologous groups")
p = p+theme_minimal()+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.text=element_text(size=12))
p = p + scale_fill_brewer(palette="Set2")
p
# dev.off()
