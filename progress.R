### thesis writing progess
library(ggplot2)
chapter <- rep(c("01-Introduction","02-Estimation of LCA protein set",
                 "03-PhyloProfile",
                 "04-Distribution analysis","05-HamFAS",
                 "06-Metabolic pathway analysis",
                 "07-Discussion & Outlook"), each = 3)
type <- rep(c("draft","refined","illustration"),5)
progress <- c(100,100,100, # chapter 01 Introduction
              100,100,95, # chapter 02 LCA estimation
              50,1,80, # chapter 03 PhyloProfile
              100,100,90, # chapter 04 Distribution analysis
              100,95,90, # chapter 05 HamFAS
              80,1,50, # chapter 06 pathway analysis
              1,1,1 # chapter 07 Discussion & Outlook
              )

df <- as.data.frame(cbind(chapter,type,progress), stringsAsFactors = FALSE)
df$progress <- as.numeric(df$progress)
df$type <- factor(df$type, levels = c("draft","refined","illustration"))

ggplot(df, aes(x = chapter, y = progress)) +
  geom_bar(aes(fill = type), stat = "identity", position = position_dodge(width = 0.7), width = 0.5) +
  geom_text(aes(label = paste0(progress,"%"), group = type), color = "white", size = 2, position = position_dodge(width = 0.7), hjust = 1) +
  coord_flip() +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "top", legend.title = element_blank())

ggsave("/Users/trvinh/work/thesis/microsporidia/progress.pdf", width = 20, height = 10, units = "cm")

