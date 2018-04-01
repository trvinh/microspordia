library(ggplot2)
library(grid)
library(gridExtra)
library(treemapify)
library(plyr)

setwd("/Users/trvinh/work/thesis/microsporidia/koAnnotation")

df <- as.data.frame(read.table("pathways_compare.stat", sep='\t',header=TRUE,check.names=FALSE,comment.char="",stringsAsFactors = FALSE))
koDf <- df[df$type == "KO",]
rnDf <- df[df$type == "RN",]
rnDf <- rnDf[rnDf$count > 0,]
textSize <- 15

### treemap
head(koDf)

ko_plot <- ggplot(koDf, aes(area = count, fill = category, label = pathway)) +
  geom_treemap() +
  geom_treemap_text(grow = T, reflow = T, colour = "black") +
  facet_wrap( ~ source) +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "bottom", legend.text = element_text(size=textSize), legend.title = element_text(size=textSize),
        strip.text.x = element_text(size = textSize),
        plot.title = element_text(size=textSize)
  ) +
  labs(
    # title = "Pathway enrichment",
    fill = "Category"
  )
# ko_plot
ggsave("pathway_enrichment_ko.pdf", width = 30, height = 20, units = "cm")


ko_plot <- ggplot(koDf[koDf$source == "01_LCA_Microsporidia",], aes(area = count, fill = category, label = pathway)) +
  geom_treemap() +
  geom_treemap_text(grow = T, reflow = T, colour = "black") +
  facet_wrap( ~ source) +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "bottom", legend.text = element_text(size=textSize), legend.title = element_text(size=textSize),
        strip.text.x = element_text(size = textSize),
        plot.title = element_text(size=textSize)
  ) +
  labs(
    # title = "Pathway enrichment",
    fill = "Category"
  )
# ko_plot
ggsave("pathway_enrichment_ko_LCA.pdf", width = 30, height = 20, units = "cm")

rn_plot <- ggplot(rnDf, aes(area = count, fill = category, label = pathway)) +
  geom_treemap() +
  geom_treemap_text(grow = T, reflow = T, colour = "black") +
  facet_wrap( ~ source) +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "bottom", legend.text = element_text(size=textSize), legend.title = element_text(size=textSize),
        strip.text.x = element_text(size = textSize),
        plot.title = element_text(size=textSize)
        ) +
  labs(
    title = "Pathway enrichment",
    fill = "Category"
  )
# rn_plot
ggsave("pathway_enrichment_rn.pdf", width = 30, height = 20, units = "cm")
