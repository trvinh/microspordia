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
sumKoDf <- aggregate(koDf$count, by=list(category=koDf$category,source=koDf$source), FUN=sum)
colnames(sumKoDf) <- c("category","source","count")
# sumKoDf <- ddply(koDf, .(category, source), summarize, count=sum(count)) # same as aggregate
sumKoDf <- ddply(sumKoDf, .(source), mutate, percentage = round(count/sum(count)*100,2))
sumKoDf


ndeep = 10

createWafflesData <- function(ndeep,sourceName, df){
  dfSub <- df[df$source == sourceName,]
  tb4waffles <- expand.grid(y = 1:ndeep,
      x = seq_len(ceiling(sum(dfSub$percentage) / ndeep)))
  catvec <- rep(dfSub$category,dfSub$percentage)
  tb4waffles$category <- c(catvec, rep(NA, nrow(tb4waffles) - length(catvec)))
  tb4waffles$source <- sourceName
  return(head(tb4waffles,n=100))
}

tb4wafflesLca <- createWafflesData(ndeep,"01_LCA_Microsporidia",sumKoDf)
tb4wafflesEnccu <- createWafflesData(ndeep,"02_E.cuniculi",sumKoDf)
tb4wafflesEnche <- createWafflesData(ndeep,"03_E.hellem",sumKoDf)
tb4wafflesEncin <- createWafflesData(ndeep,"04_E.intestinalis",sumKoDf)
tb4wafflesNosce <- createWafflesData(ndeep,"05_N.ceranae",sumKoDf)
tb4wafflesSacce <- createWafflesData(ndeep,"06_S.cerevisiae",sumKoDf)

tb4waffles <- rbind(tb4wafflesLca,tb4wafflesEnccu)
tb4waffles <- rbind(tb4waffles,tb4wafflesEnche)
tb4waffles <- rbind(tb4waffles,tb4wafflesEncin)
tb4waffles <- rbind(tb4waffles,tb4wafflesNosce)
tb4waffles <- rbind(tb4waffles,tb4wafflesSacce)
tb4waffles$category[is.na(tb4waffles$category)] <- "Metabolism"

wafflePlot <- ggplot(tb4waffles, aes(x = x, y = y, fill = category)) + 
  geom_tile(color = "white") + # The color of the lines between tiles +
  facet_wrap( ~ source) +
  scale_fill_brewer(name="Pathway category",palette = "Set2")
wafflePlot <- wafflePlot +
  theme(#panel.border = element_rect(size = 2),
        #plot.title = element_text(size = rel(1.2)),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",legend.text = element_text(size=textSize),legend.title = element_blank(),
        strip.text.x = element_text(size = textSize))
wafflePlot

# p <- ggplot(sumKoDf, aes(x = category, y = count, fill = source)) + 
#   geom_bar(stat = 'identity', position = 'dodge') + 
#   # coord_flip() +
#   facet_grid(. ~ category, scales = "free", space = "free") +
#   labs(y="Number of proteins")
# p <- p  + scale_fill_brewer(palette = "Set2") + 
#   # theme_minimal() +
#   theme(legend.title = element_blank(), legend.text = element_text(size=textSize),legend.position="top",
#         axis.text.x = element_blank(),
#         axis.title.x = element_blank(), axis.title.y = element_text(size=textSize),
#         strip.text.x = element_text(size = textSize)
#   )
# p
ggsave("categories_comparison_ko.pdf", width = 30, height = 20, units = "cm")

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
    title = "Pathway enrichment",
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
    title = "Pathway enrichment",
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
