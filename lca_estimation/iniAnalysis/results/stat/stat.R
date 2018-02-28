setwd("/Users/trvinh/work/thesis/microsporidia/lca_estimation/iniAnalysis/results/stat")
data <- read.table("orthomcl_micros_orig.list.stat", sep='\t',header=T)

library(ggplot2)
library(reshape2) # for melt
library(RColorBrewer)
library(ggthemes)

data$Taxon <- ordered(data$Taxon, 
# levels = c("E.intestinalis","E.hellem","E.cuniculi","N.ceranae","E.bieneusi","V.corneae","A.locustae","A.algerae","V.culicis","E.aedis","N.parisii"))
levels = c("E.hellem","E.intestinalis","E.cuniculi","N.ceranae","E.bieneusi","V.corneae","A.algerae","A.locustae","E.aedis","V.culicis","N.parisii"))
head(data)
melted <- melt(data, "Taxon")

# png("/home/vinh/Desktop/data/project/iniAnalysis/orthomcl_micros_orig.list.stat.png", width=1600, height=800)
melted$cat <- ''
melted[melted$variable == 'LCA_homologous_prots',]$cat <- "#Genes"
melted[melted$variable == 'Non.LCA_homologous_prots',]$cat <- "#Genes"
melted[melted$variable == 'Orphan_prots',]$cat <- "#Genes"
melted[melted$variable == 'LCA.OG',]$cat <- "#Groups"
melted[melted$variable == 'Non.LCA.OG',]$cat <- "#Groups"
melted[melted$variable == 'Total_homologous',]$cat <- "#Genes"
melted[melted$variable == 'Orthologous_groups',]$cat <- "#Groups"


p = ggplot(melted, aes(x = cat, y = value, fill = variable)) + 
    geom_bar(stat = 'identity', position = 'stack') + 
    facet_grid(~ Taxon) +
    theme_set(theme_gray(base_size = 15))
p = p+labs(x="", y="")
p = p+theme(axis.text.x = element_blank(),axis.text.y = element_text(size=15),strip.text.x = element_text(size = 18),legend.text=element_text(size=15)) + 
    theme(legend.title=element_blank())
p = p + scale_fill_brewer(palette="Set2")
p
# dev.off()

melted
meltedSub <- melted[melted$variable == "Total_homologous" | melted$variable == "Orphan_prots",]
meltedSub$value2[meltedSub$variable == "Orphan_prots"] <- 0-meltedSub$value[meltedSub$variable == "Orphan_prots"]
meltedSub$value2[meltedSub$variable == "Total_homologous"] <- 0+meltedSub$value[meltedSub$variable == "Total_homologous"]
meltedSub$variable2[meltedSub$variable == "Orphan_prots"] <- "non-orthologous protein"
meltedSub$variable2[meltedSub$variable == "Total_homologous"] <- "orthologous protein"
meltedSub

# X Axis Breaks and Labels 
brks <- seq(-2500, 2500, 500)
lbls = paste0(as.character(c(seq(2500, 0, -500), seq(500, 2500, 500))))

meltedSub$Taxon <- factor(meltedSub$Taxon, 
    levels = rev(c("E.intestinalis","E.hellem","E.cuniculi","N.ceranae","E.bieneusi","V.corneae","A.locustae","A.algerae","V.culicis","E.aedis","N.parisii")))

ggplot(meltedSub, aes(x = Taxon, y = value2, fill = variable2)) +   # Fill column
  geom_bar(stat = "identity", width = .6) +   # draw the bars
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls) + # Labels
  coord_flip() +  # Flip axes
  labs(#title="Fraction of homologous and non-homologous proteins",
       y = "Number of proteins",
       fill = "") +
  # theme_tufte() +  # Tufte theme from ggfortify
  theme_minimal() +
  theme(plot.title = element_text(hjust = .5), 
        axis.ticks = element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=8)) +   # Centre plot title
  scale_fill_brewer(palette = "Set2")  # Color palette

meltedSub2 <- melted[melted$variable == "Total_homologous" | melted$variable == "Orthologous_groups",]
meltedSub2$variable2[meltedSub2$variable == "Total_homologous"] <- "Homologous proteins"
meltedSub2$variable2[meltedSub2$variable == "Orthologous_groups"] <- "Homologous groups"


p = ggplot(meltedSub2, aes(x = cat, y = value, fill = variable2)) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  facet_grid(~ Taxon) +
  theme_set(theme_gray(base_size = 15))
p = p+labs(x="", y="")
p = p+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        strip.text.x = element_text(size = 10),
        legend.text=element_text(size=10)) + 
  theme(legend.title=element_blank(),
        legend.position = "bottom")
p = p + scale_fill_brewer(palette="Set2")
p