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

### waffle plot
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
