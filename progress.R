### thesis writing progess
library(ggplot2)
chapter <- rep(c("01-Introduction","02-Estimation of LCA protein set",
                 "03-PhyloProfile",
                 "04-Distribution analysis",
                 "05-Metabolic pathway analysis"), each = 3)
type <- rep(c("draft","refining","illustration"),5)
progress <- c(100,100,100, # chapter 01
              100,98,80, # chapter 02
              1,1,10, # chapter 03
              70,1,80, # chapter 04
              80,1,70 # chapter 05
              )

df <- as.data.frame(cbind(chapter,type,progress), stringsAsFactors = FALSE)
df$progress <- as.numeric(df$progress)
df$type <- factor(df$type, levels = c("draft","refining","illustration"))

ggplot(df, aes(x = chapter, y = progress)) +
  geom_bar(aes(fill = type), stat = "identity", position = position_dodge(width = 0.7), width = 0.5) +
  geom_text(aes(label = paste0(progress,"%"), group = type), color = "white", position = position_dodge(width = 0.7), hjust = 1) +
  coord_flip() +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "top")
