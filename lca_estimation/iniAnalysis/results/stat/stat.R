setwd("/Users/trvinh/work/thesis/microsporidia/lca_estimation/iniAnalysis/results/stat")
data <- read.table("orthomcl_micros_orig.list.stat", sep='\t',header=T)

library(ggplot2)
library(reshape2) # for melt
library(RColorBrewer)
library(ggthemes)
library(grid)
library(gridExtra)
library(cowplot)

data$Taxon <- ordered(data$Taxon, 
# levels = c("E.intestinalis","E.hellem","E.cuniculi","N.ceranae","E.bieneusi","V.corneae","A.locustae","A.algerae","V.culicis","E.aedis","N.parisii"))
levels = c("E.hellem","E.intestinalis","E.cuniculi","N.ceranae","E.bieneusi","V.corneae","A.algerae","A.locustae","E.aedis","V.culicis","N.parisii"))
head(data)
melted <- melt(data, "Taxon")
head(melted)

##### Pyramid plot
dataTMP <- data[order(data[,"Genome_Size"]),]
dataTMP$Taxon
dataTMP$Taxon <- factor(dataTMP$Taxon, levels = dataTMP$Taxon)


pyramidDf <- dataTMP[,c("Taxon","Orphan_prots","Total_homologous")]
pyramidDf$Orphan_prots <- 0-pyramidDf$Orphan_prots
colnames(pyramidDf) <- c("Taxon","Orphan","Orthologous protein")
meltedPyramid <- melt(pyramidDf, "Taxon")
head(meltedPyramid)

# X Axis Breaks and Labels 
brks <- seq(-2000, 2500, 500)
lbls = paste0(as.character(c(seq(2000, 0, -500), seq(500, 2500, 500))))

# Plot
pyra <- ggplot(meltedPyramid, aes(x = Taxon, y = value, fill = variable)) +   # Fill column
  geom_bar(stat = "identity", width = .6) +   # draw the bars
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls) + # Labels
  scale_x_discrete(limits = rev(levels(meltedPyramid$Taxon))) + # reverse x-axis order
  coord_flip() +  # Flip axes
  theme_minimal() +  # Tufte theme from ggfortify
  labs(y = "Number of proteins") +
  theme(plot.title = element_text(hjust = .5), 
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_blank(),
        legend.position = "right") +   # Centre plot title
  scale_fill_brewer(palette = "Set2")  # Color palette
pyra

sizeDf <- dataTMP[,c("Taxon","Genome_Size","Total_Protein")]
# colnames(sizeDf) <- c("Taxon","Genome size")
# meltedSizeDf <- melt(sizeDf, "Taxon")
# head(meltedSizeDf)
head(sizeDf)
genome <- ggplot(sizeDf, aes(x = Taxon, y = Genome_Size)) +
  geom_bar(stat = "identity", width = .6, fill="orange") +
  scale_x_discrete(limits = rev(levels(sizeDf$Taxon))) + # reverse x-axis order
  coord_flip() +  # Flip axes
  theme_minimal() +  # Tufte theme from ggfortify +
  labs(y = "Genome size (Mb)") +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.position = "top")   # Centre plot title
genome

prot <- ggplot(sizeDf, aes(x = Taxon, y = Total_Protein)) +
  geom_bar(stat = "identity", width = .6, fill="steelblue") +
  scale_x_discrete(limits = rev(levels(sizeDf$Taxon))) + # reverse x-axis order
  coord_flip() +  # Flip axes
  scale_y_reverse() +
  theme_minimal() +  # Tufte theme from ggfortify +
  labs(y = "Number of predicted gene") +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.position = "top")   # Centre plot title
prot

# grid.arrange(pyra, size, nrow = 1, rel_widths = c(2/3, 1/3))
plot_grid(prot,pyra,genome, nrow = 1, labels = "", align = 'h', rel_widths = c(0.3,1,0.3))

#####


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
