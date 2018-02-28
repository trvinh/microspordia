setwd("/Users/trvinh/work/thesis/microsporidia/lca_estimation/iniAnalysis/results/pfam")
data <- read.table("orthomcl_micros_orig.list.pfam", sep='\t',header=T)

library(ggplot2)
library(reshape2) # for melt
library(RColorBrewer)

data$Taxon <- ordered(data$Taxon, 
#levels = c("E.hellem","E.intestinallis","E.cuniculi","N.ceranae","E.bieneusi","V.corneae","A.algerae","A.locustae","E.aedis","V.culicis","N.parisii"))
levels = c("E.hellem","E.intestinalis","E.cuniculi","N.ceranae","E.bieneusi","V.corneae","A.algerae","A.locustae","E.aedis","V.culicis","N.parisii"))

melted <- melt(data, "Taxon")

# png("/home/vinh/Desktop/data/project/iniAnalysis/orthomcl_micros_orig.list.pfam.png", width=1600, height=800)
melted$cat <- ''
melted[melted$variable == 'homologous_Pfam',]$cat <- "homologous_prots"
melted[melted$variable == 'homologous_noPfam',]$cat <- "homologous_prots"
melted[melted$variable == 'orphan_homologPfam',]$cat <- "orphans"
melted[melted$variable == 'orphan_newPfam',]$cat <- "orphans"
melted[melted$variable == 'orphan_noPfam',]$cat <- "orphans"

melted$variable2[melted$variable == 'homologous_Pfam'] <- "orthologs have PFAM"
melted$variable2[melted$variable == 'homologous_noPfam'] <- "orthologs have no PFAM"
melted$variable2[melted$variable == 'orphan_homologPfam'] <- "orphans have orthologous PFAM"
melted$variable2[melted$variable == 'orphan_newPfam'] <- "orphans have new PFAM"
melted$variable2[melted$variable == 'orphan_noPfam'] <- "orphans have no PFAM"


p = ggplot(melted, aes(x = cat, y = value, fill = variable2)) + 
    geom_bar(stat = 'identity', position = 'stack') + 
    facet_grid(~ Taxon) +
    theme_set(theme_gray(base_size = 18)
)

p = p+labs(x="", y="")

p = p+ theme_minimal() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        strip.text.x = element_text(size = 10)
        ) + 
    theme(legend.title=element_blank(),
          legend.text = element_text(size=10),
          legend.position="bottom")+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

p = p + scale_fill_brewer(palette="Set2")
p
# dev.off()
