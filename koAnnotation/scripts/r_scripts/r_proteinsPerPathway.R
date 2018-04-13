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

### barplot
# warp text of category
swr = function(string, nwrap=20) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)

koDf$category = swr(koDf$category)
# create plot
p <- ggplot(koDf, aes(x = pathway, y = count, fill = source)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  # coord_flip() +
  facet_grid(. ~ category, scales = "free", space = "free") +
  labs(y="Number of proteins")
p <- p  + scale_fill_brewer(palette = "Set2") +
  # theme_minimal() +
  theme(legend.title = element_blank(), legend.text = element_text(size=textSize*1.3),legend.position="top",
        axis.text.x = element_text(angle=70,hjust=1,size=textSize*1.5),
        axis.title.x = element_blank(), axis.title.y = element_text(size=textSize*1.5),
        strip.text.x = element_text(size = textSize*1.5)
  )
p
ggsave("pathway_comparison_ko.pdf", width = 50, height = 30, units = "cm")

rnDf$category = swr(rnDf$category)
# create plot
p2 <- ggplot(rnDf, aes(x = pathway, y = count, fill = source)) +
  geom_bar(stat = 'identity', position = 'dodge', width = 0.8) +
  # coord_flip() +
  facet_grid(. ~ category, scales = "free", space = "free") +
  labs(y="Number of reactions")
p2 <- p2  + scale_fill_brewer(palette = "Set2") +
  # theme_minimal() +
  theme(legend.title = element_blank(), legend.text = element_text(size=textSize),legend.position="top",
        axis.text.x = element_text(angle=60,hjust=1,size=textSize),
        axis.title.x = element_blank(), axis.title.y = element_text(size=textSize),
        strip.text.x = element_text(size = textSize)
  )
# p2
ggsave("pathway_comparison_rn.pdf", width = 50, height = 20, units = "cm")
